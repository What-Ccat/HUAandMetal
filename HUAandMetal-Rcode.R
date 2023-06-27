##################################################
##                                              ##
##    变量选择——常用的Lasso和Elastic Net 分析   ##
##                                              ##
##    原作者：Ahlam Abuawad    编译：What！Cat  ##
##                                              ##
##################################################

#install_github("gabrielrvsc/HDeconometrics")

#### 加载需要的package####
library(tidyverse)
library(janitor)
library(caret)
library(Hmisc)
library(glmnet)
library(ggplot2)
library(HDeconometrics)
library(plotmo)

#### 数据的导入与清洗 ####
alldata <- read.csv("读取文件")
names(alldata) ##看看有几个数据
options(scipen = 999) ## 关闭科学计数法

#### 根据需要清理相应的数据
cleandata <- na.omit(alldata)

## 建立Lasso回归需要的自变量和因变量（需要完整的数据集）
## 创建一个预测变量矩阵作为 x
x = model.matrix(结局变量 ~ ., 数据集)[,-1] 
y = cleandata$结局变量 ## 提取结果向量
colnames(x)

#####快速描述清洗好的数据#####
describe(cleandata) ##查看描述性数据
table(cleandata$结局变量) ##查看数据，看看是否符合预期

#####数据可视化######
## 先绘制经过Log转换的暴露变量x
names(x)
featurePlot(x = x[,14:20],
            y = y,
            between = list(x = 1, y = 1), 
            type = c("g", "p", "smooth"))

## 再绘制协变量的分布图
featurePlot(x = x[,1:13],
            y = y,
            between = list(x = 1, y = 1),
            type = c("g", "p", "smooth"))
# 警告！这是关于绘图参数的，可以忽略

######利用交叉验证Lasso分析进行选择#######
## 我们将首先查看应用于完整设计矩阵的Lasso。下面，我们将使用交叉验证法来识别调整参数；
## 这是一个随机过程，因此我们将set.seed以确保可重复性。
set.seed(2022)

## 首先，我们将使用指定的调整参数值网格，并为每个值拟合套索模型。
lam_grid <- .5 ^ (-20:20)

lasso.mod = glmnet(x, y, alpha = 1,family = c( "binomial"), lambda = lam_grid)

## 绘制系数路径图
plot(lasso.mod)

## 查看数据的系数
coef(lasso.mod)[,10]

## 一些内置函数将进行交叉验证分析并确定“最佳”调整参数。
set.seed(2022)

## cv.glmnet 的 n-folds 设置为默认值 10，即10折交叉验证
cv.out = cv.glmnet(x, y,family = c( "binomial"), alpha = 1,nfolds = 10)

## 绘制均方误差和log（λ）关系图
plot(cv.out)

## 查看经过交叉验证法调整后的参数
coef(cv.out)

## 找到导致最小CV误差的λ
best_lambda = cv.out$lambda.1se
best_lambda
log(best_lambda)

## 建立ln(λ)绘图方程——改自CSDN的hetallian
get_plot<- function(the_fit,the_fit_cv,the_lamb,toplot = seq(1,50,2),ylim,col,legend){
  Coefficients <- coef(the_fit, s = the_lamb)
  coeall <- coef(the_fit, s = the_fit_cv$lambda[toplot])
  coe <- coeall
  sp <- spline(log(the_fit_cv$lambda[toplot]),coe[1,],n=100)
  plot(sp,type='l',col =1,lty=1,ylim=ylim,
       ylab = 'Coefficients', xlab = 'log Lambda')
  abline(v=log(best_lambda),col="red",lwd=1.5)
  abline(h=0,lty=2) 
  if (col==T) {
    for(i in c(2:nrow(coe))){
      lines(spline(log(the_fit_cv$lambda[toplot]),coe[i,],n=1000),
            col =i)
    }} 
  else if (col==F){
    for(i in c(2:nrow(coe))){
      lines(spline(log(the_fit_cv$lambda[toplot]),coe[i,],n=1000),
            col ="#797979")
    }}
  
  if (legend==T){
    legend("right",legend=rownames(coe),col=c(1:nrow(coe)),
           lty=1,xpd=T,x.intersp=0.5,
           cex=0.5)}
}


## 引入最优λ（best_lambda)，变量筛选图
get_plot(lasso.mod,cv.out,best_lambda,ylim=c(-0.5,1),col=F,legend=F)

##或者用HDeconometrics包绘制Lasso图,画出来好看些的其实##
## 重新拟合lasso，使用HDeconometrics的自带语句
lasso=ic.glmnet(x, y,family = c( "binomial"),crit = "bic")

## 绘制图形
plot(lasso$glmnet,"lambda",ylim=c(-0.5,1))
abline(v=log(best_lambda),col="red",lwd=1.5)
abline(h=0,lty=2)

## 利用plomo包的语句画图
plot_glmnet(lasso.mod,xvar="lambda",label=10,nresponse = 3) 

## 输出参数的函数
coeff2dt <- function(fitobject, s) {
  coeffs <- coef(fitobject, s) 
  coeffs.dt <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x) 
  
  # reorder the variables in term of coefficients
  return(coeffs.dt)
}

coeffs.table <- coeff2dt(fitobject = cv.out, s = "lambda.1se")
coeffs.table = coeffs.table[-1,]


## 画图
ggplot(data = coeffs.table) +
  geom_col(aes(x = name, y = coefficient, fill = {coefficient > 0})) +
  xlab(label = "") + ylim(-0.1,0.7)+
  ggtitle(expression(paste("Lasso Coefficients with ", lambda, " = 0.0026"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


########################################
##                                    ##
##       贝叶斯核机器回归（BKMR）     ##
##                                    ##
########################################

########加载包########
library(bkmr)
library(ggplot2)

### 指定路径
alldata <- read.csv("文件路径")

### 清除缺失
cleandata <- na.omit(subdata)

### 定义暴露结局协变量

covar <- data.matrix(cleandata[,c("协变量1","协变量2")])

expos <- data.matrix(cleandata[,c("暴露1","暴露2")])

Y <- cleandata$结局

### 拟合模型，基本不用改
fitkm <- kmbayes(Y, Z = expos, X = covar, iter = 20000,family = "binomial",verbose = FALSE, varsel = TRUE)

pred.resp.univar <- PredictorResponseUnivar(fit = fitkm)

expos.pairs <- subset(data.frame(expand.grid(expos1 = c(1,2,3,4,5),expos2 = c(1,2,3,4,5))), expos1 < expos2)

pred.resp.bivar <- PredictorResponseBivar(fit = fitkm,min.plot.dist = 0.5,z.pairs = expos.pairs)

pred.resp.bivar.levels <- PredictorResponseBivarLevels(pred.resp.bivar,expos,qs = c(0.10, 0.5, 0.90))

risks.overall <- OverallRiskSummaries(fit = fitkm, qs = seq(0.10, 0.90, by = 0.10), q.fixed = 0.50)

risk.singvar <- SingVarRiskSummaries(fit = fitkm,qs.diff = c(0.10,0.90),q.fixed = c(0.10,0.5,0.90))

risk.int <- SingVarIntSummaries(fit=fitkm, qs.diff = c(0.10,0.90),qs.fixed= c(0.10,0.90))

######调查模型的收敛情况######
TracePlot(fit = fitkm, par = "beta") ##查看不同的参数随着随机数的变化
TracePlot(fit = fitkm, par = "sigsq.eps")
TracePlot(fit = fitkm, par = "r", comp = 5)

#####提取模型参数,估计后验证包含概率Pip######
ExtractPIPs(fitkm)  ##pip值越高，其越在模型中就越重要


#####绘制单变量截面图######
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, 
                             ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~variable, ncol = 5) +
  xlab("expos") +
  ylab("h(expos)")

#####根据单变量截面图绘制双变量截面图######
#######输出变量轮廓图#####
ggplot(pred.resp.bivar, aes(z1, z2, fill = est)) + 
  geom_raster() + 
  facet_grid(variable2 ~ variable1) +
  scale_fill_gradientn(colours = c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  xlab("expos1") +
  ylab("expos2") +
  ggtitle("h(expos1, expos2)")

#########对变量两两组合的暴露-反应曲线进行绘制######
#####绘制暴露-反应曲线######
ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")

#######在双变量的基础上进行三变量绘制######
######累积效应分析######
risks.overall

######绘制累计效应图######
ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd)) + 
  scale_x_continuous(breaks = c(0.1,0.3,0.5,0.7,0.9))+
  geom_hline(yintercept = 0, lty = 2, col = "brown") +
  geom_pointrange()

######绘制单因素暴露的效果#####
#####绘制各百分位数时的贡献####
ggplot(risk.singvar,aes(variable,est,ymin=est-1.96*sd,ymax=est +1.96*sd,col=q.fixed))+
  geom_pointrange(position=position_dodge(width=0.75))+geom_hline(yintercept = 0, lty = 2, col="red")+
  coord_flip()

######变量的交互作用#####
#####绘制效应图####
ggplot(risk.int,aes(variable,est,ymin =est - 1.96*sd,ymax =est+1.96*sd))+
  geom_pointrange(position=position_dodge(width=0.75))+
  geom_hline(yintercept = 0, lty = 2, col="brown")+
  coord_flip()

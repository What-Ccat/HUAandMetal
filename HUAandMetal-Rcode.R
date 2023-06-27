##################################################
##                                              ##
##    ����ѡ�񡪡����õ�Lasso��Elastic Net ����   ##
##                                              ##
##    ԭ���ߣ�Ahlam Abuawad    ���룺What��Cat  ##
##                                              ##
##################################################

#install_github("gabrielrvsc/HDeconometrics")

#### ������Ҫ��package####
library(tidyverse)
library(janitor)
library(caret)
library(Hmisc)
library(glmnet)
library(ggplot2)
library(HDeconometrics)
library(plotmo)

#### ���ݵĵ�������ϴ ####
alldata <- read.csv("��ȡ�ļ�")
names(alldata) ##�����м�������
options(scipen = 999) ## �رտ�ѧ������

#### ������Ҫ������Ӧ������
cleandata <- na.omit(alldata)

## ����Lasso�ع���Ҫ���Ա��������������Ҫ���������ݼ���
## ����һ��Ԥ�����������Ϊ x
x = model.matrix(��ֱ��� ~ ., ���ݼ�)[,-1] 
y = cleandata$��ֱ��� ## ��ȡ�������
colnames(x)

#####����������ϴ�õ�����#####
describe(cleandata) ##�鿴����������
table(cleandata$��ֱ���) ##�鿴���ݣ������Ƿ����Ԥ��

#####���ݿ��ӻ�######
## �Ȼ��ƾ���Logת���ı�¶����x
names(x)
featurePlot(x = x[,14:20],
            y = y,
            between = list(x = 1, y = 1), 
            type = c("g", "p", "smooth"))

## �ٻ���Э�����ķֲ�ͼ
featurePlot(x = x[,1:13],
            y = y,
            between = list(x = 1, y = 1),
            type = c("g", "p", "smooth"))
# ���棡���ǹ��ڻ�ͼ�����ģ����Ժ���

######���ý�����֤Lasso��������ѡ��#######
## ���ǽ����Ȳ鿴Ӧ����������ƾ����Lasso�����棬���ǽ�ʹ�ý�����֤����ʶ�����������
## ����һ��������̣�������ǽ�set.seed��ȷ�����ظ��ԡ�
set.seed(2022)

## ���ȣ����ǽ�ʹ��ָ���ĵ�������ֵ���񣬲�Ϊÿ��ֵ�������ģ�͡�
lam_grid <- .5 ^ (-20:20)

lasso.mod = glmnet(x, y, alpha = 1,family = c( "binomial"), lambda = lam_grid)

## ����ϵ��·��ͼ
plot(lasso.mod)

## �鿴���ݵ�ϵ��
coef(lasso.mod)[,10]

## һЩ���ú��������н�����֤������ȷ������ѡ�����������
set.seed(2022)

## cv.glmnet �� n-folds ����ΪĬ��ֵ 10����10�۽�����֤
cv.out = cv.glmnet(x, y,family = c( "binomial"), alpha = 1,nfolds = 10)

## ���ƾ�������log���ˣ���ϵͼ
plot(cv.out)

## �鿴����������֤��������Ĳ���
coef(cv.out)

## �ҵ�������СCV���Ħ�
best_lambda = cv.out$lambda.1se
best_lambda
log(best_lambda)

## ����ln(��)��ͼ���̡�������CSDN��hetallian
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


## �������Ŧˣ�best_lambda)������ɸѡͼ
get_plot(lasso.mod,cv.out,best_lambda,ylim=c(-0.5,1),col=F,legend=F)

##������HDeconometrics������Lassoͼ,�������ÿ�Щ����ʵ##
## �������lasso��ʹ��HDeconometrics���Դ����
lasso=ic.glmnet(x, y,family = c( "binomial"),crit = "bic")

## ����ͼ��
plot(lasso$glmnet,"lambda",ylim=c(-0.5,1))
abline(v=log(best_lambda),col="red",lwd=1.5)
abline(h=0,lty=2)

## ����plomo������仭ͼ
plot_glmnet(lasso.mod,xvar="lambda",label=10,nresponse = 3) 

## ��������ĺ���
coeff2dt <- function(fitobject, s) {
  coeffs <- coef(fitobject, s) 
  coeffs.dt <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x) 
  
  # reorder the variables in term of coefficients
  return(coeffs.dt)
}

coeffs.table <- coeff2dt(fitobject = cv.out, s = "lambda.1se")
coeffs.table = coeffs.table[-1,]


## ��ͼ
ggplot(data = coeffs.table) +
  geom_col(aes(x = name, y = coefficient, fill = {coefficient > 0})) +
  xlab(label = "") + ylim(-0.1,0.7)+
  ggtitle(expression(paste("Lasso Coefficients with ", lambda, " = 0.0026"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


########################################
##                                    ##
##       ��Ҷ˹�˻����ع飨BKMR��     ##
##                                    ##
########################################

########���ذ�########
library(bkmr)
library(ggplot2)

### ָ��·��
alldata <- read.csv("�ļ�·��")

### ���ȱʧ
cleandata <- na.omit(subdata)

### ���屩¶���Э����

covar <- data.matrix(cleandata[,c("Э����1","Э����2")])

expos <- data.matrix(cleandata[,c("��¶1","��¶2")])

Y <- cleandata$���

### ���ģ�ͣ��������ø�
fitkm <- kmbayes(Y, Z = expos, X = covar, iter = 20000,family = "binomial",verbose = FALSE, varsel = TRUE)

pred.resp.univar <- PredictorResponseUnivar(fit = fitkm)

expos.pairs <- subset(data.frame(expand.grid(expos1 = c(1,2,3,4,5),expos2 = c(1,2,3,4,5))), expos1 < expos2)

pred.resp.bivar <- PredictorResponseBivar(fit = fitkm,min.plot.dist = 0.5,z.pairs = expos.pairs)

pred.resp.bivar.levels <- PredictorResponseBivarLevels(pred.resp.bivar,expos,qs = c(0.10, 0.5, 0.90))

risks.overall <- OverallRiskSummaries(fit = fitkm, qs = seq(0.10, 0.90, by = 0.10), q.fixed = 0.50)

risk.singvar <- SingVarRiskSummaries(fit = fitkm,qs.diff = c(0.10,0.90),q.fixed = c(0.10,0.5,0.90))

risk.int <- SingVarIntSummaries(fit=fitkm, qs.diff = c(0.10,0.90),qs.fixed= c(0.10,0.90))

######����ģ�͵��������######
TracePlot(fit = fitkm, par = "beta") ##�鿴��ͬ�Ĳ�������������ı仯
TracePlot(fit = fitkm, par = "sigsq.eps")
TracePlot(fit = fitkm, par = "r", comp = 5)

#####��ȡģ�Ͳ���,���ƺ���֤��������Pip######
ExtractPIPs(fitkm)  ##pipֵԽ�ߣ���Խ��ģ���о�Խ��Ҫ


#####���Ƶ���������ͼ######
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, 
                             ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~variable, ncol = 5) +
  xlab("expos") +
  ylab("h(expos)")

#####���ݵ���������ͼ����˫��������ͼ######
#######�����������ͼ#####
ggplot(pred.resp.bivar, aes(z1, z2, fill = est)) + 
  geom_raster() + 
  facet_grid(variable2 ~ variable1) +
  scale_fill_gradientn(colours = c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  xlab("expos1") +
  ylab("expos2") +
  ggtitle("h(expos1, expos2)")

#########�Ա���������ϵı�¶-��Ӧ���߽��л���######
#####���Ʊ�¶-��Ӧ����######
ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")

#######��˫�����Ļ����Ͻ�������������######
######�ۻ�ЧӦ����######
risks.overall

######�����ۼ�ЧӦͼ######
ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd)) + 
  scale_x_continuous(breaks = c(0.1,0.3,0.5,0.7,0.9))+
  geom_hline(yintercept = 0, lty = 2, col = "brown") +
  geom_pointrange()

######���Ƶ����ر�¶��Ч��#####
#####���Ƹ��ٷ�λ��ʱ�Ĺ���####
ggplot(risk.singvar,aes(variable,est,ymin=est-1.96*sd,ymax=est +1.96*sd,col=q.fixed))+
  geom_pointrange(position=position_dodge(width=0.75))+geom_hline(yintercept = 0, lty = 2, col="red")+
  coord_flip()

######�����Ľ�������#####
#####����ЧӦͼ####
ggplot(risk.int,aes(variable,est,ymin =est - 1.96*sd,ymax =est+1.96*sd))+
  geom_pointrange(position=position_dodge(width=0.75))+
  geom_hline(yintercept = 0, lty = 2, col="brown")+
  coord_flip()
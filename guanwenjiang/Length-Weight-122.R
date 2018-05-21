library(readr)
RData<- read_csv("~/R data/guanwenjiang/data/Length-Weight.csv")
####lm()函数与lsfit()函数采用最小二乘法进行线性回归
names(RData)
length<-log(RData[,1])
weight<-log(RData[,2])
length<- unlist(length)#使用length<- unlist(length)转换类型，然后进行回归
weight<- unlist(weight)#使用weight<- unlist(weight)转换类型，然后进行回归
lmF<- lm(weight~length)
#lsF<-lsfit(length,weight)
#coef(lsF)
coefficients(lmF)
summary(lmF)
############################################
#正太分布最大似然估计
#读取数据
RData<- read_csv("~/R data/guanwenjiang/data/Length-Weight.csv")
length<-log(RData[,1])
weight<-log(RData[,2])
length<- unlist(length)#使用length<- unlist(length)转换类型，然后进行回归
weight<- unlist(weight)#使用weight<- unlist(weight)转换类型，然后进行回归
#编写负对数似然函数
MinusLogL<- function(Sigma=log(0.2),alpha=2.0,beta=2.0)
{
  n<-length(weight)
  Sigma2=2*exp(Sigma)^2
  n*Sigma+0.5*n*log(2*pi)+sum((weight-(alpha+beta*length))*weight-(alpha+beta*length))/Sigma2
}
#利用其密度函数，可以简化MinusLogL()函数
MinusLogL<- function(Sigma=log(0.2),alpha=-5.0,beta=3.0)
{
  -sum(stats::dnorm(weight,mean=alpha+beta*length,sd=exp(Sigma),log = TRUE) )
}
#利用stats4包中mle()函数，进行最大似然法拟合
library(stats4)
library(stats)
MLEFit<- mle(MinusLogL)
summary(MLEFit)
#利用logLik()获取对数似然值
logLik(MLEFit)
#利用vcov()函数获得参数的协方差矩阵
vcov(MLEFit)
#profile()函数提供参数概率剖面图
oldp<-par(mfcol=c(3,1))
plot(profile(MLEFit),adsVal=FALSE)
par(oldp)
#获取最大似然参数估计值
sigma<-coef(MLEFit)[1]
alpha<-coef(MLEFit)[2]
beta<-coef(MLEFit)[3]
#beta参数的取值范围
beta2<-seq(2.7,3.6,by=0.01)
N<-length(beta2)
LL<-numeric(N)
#固定beta参数，对其他参数进行最大似然估计，并获得相应的似然值
for(i in 1:N)
{
  MLEFit1<-mle(MinusLogL,fix=list(beta=beta2[i]))
  LL[i]<-logLik(MLEFit1)
}
#获取最大似然估计值
MaxLL<-MinusLogL(Sigma,alpha,beta)
#绘制似然值随beta的变化曲线
plot(beta2,LL)
#获取最小似然估计值
MinLL<-logLik(MLEFit)-qchisq(1-0.05,1)/2
#beta的置信区间
abline(h=MinLL,col="green")
#参数的置信区间也可以利用confint()函数直接获取
confint(MLEFit)
betaLU<-confint(MLEFit)[3,]
abline(v=betaLU[1],col="red")
abline(v=betaLU[2],col="red")



#对数正态分布最大似然估计
length<-log(RData[,1])
weight<-RData[,2]
length<- unlist(length)#使用length<- unlist(length)转换类型，然后进行回归
weight<- unlist(weight)
MinusLogL<- function(Sigma=log(0.2),alpha=-5.0,beta=3.0)
{
  if(Sigma>-10 & Sigma<10)
  -sum(stats::dlnorm(weight,meanlog=alpha+beta*length-0.5*exp(Sigma)*exp(Sigma),sdlog=exp(Sigma),log = TRUE) )
  else
    NA
}
MLEFit<-mle(MinusLogL)
coef(MLEFit)
summary(MLEFit)
confint(MLEFit)
logLik(MLEFit)
vcov(MLEFit)
oldp<-par(mfcol=c(3,1))
plot(profile(MLEFit),adsVal=FALSE)
par(oldp)

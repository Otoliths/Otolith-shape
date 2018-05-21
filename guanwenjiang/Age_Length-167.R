#read data
library(readr)
RData <- read_csv("~/R data/guanwenjiang/data/Age_Length.csv")
View(RData)
ObsL<- RData[,2]
Age<- RData[,1]

#似然函数
MinusLongL<- function(Linf=log(100),K=0.25,t0=-0.15,sigma=0.2)
{
  if(K>0 & K<1)
  {
    L<- exp(Lif)*(1.0-exp(-K*(Age-t0)))
    -sum(stats::dnorm(ObsL,mean=L,sd=exp(Sigma),log=TRUE))
  }
  else
    NA
}

#最大似然法估计参数
MLEFit<- mle(MinusLongL)
summary(MLEFit)
AIC(MLEFit)


#################
Age1<- 1:6
L1<- c(166.4,218.2,288.1,346.1,1402.1,440.2)
L1L<-log(L1)
#逻辑斯谛生长方程，正太分布
LogM<- nls(L1~Linf/(1+exp(alpha-r*Age1)),start = list(Linf=445.9,alpha=1.2,r=0.5))
#逻辑斯谛生长方程，对数正太分布
LogML<- nls(L1L~log(Linf)-log(1+exp(alpha-r*Age1)),start = list(Linf=445.9,alpha=1.2,r=0.5))

#Gompertz生长方程，正态分布
LogMG<- nls(L1~Linf*exp(-g*exp(-r*Age1)),start = list(Linf=445.9,g=1.2,r=0.5))
#Gompertz生长方程，对数正态分布
LogMGL<- nls(L1L~log(Linf)-g*exp(-r*Age1),start = list(Linf=445.9,g=1.2,r=0.5))

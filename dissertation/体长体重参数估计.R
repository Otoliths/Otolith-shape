data <- read.csv("./data/Schizothorax nukiangensis.csv")
head(data)
str(data)


##最小二乘法估计a，b参数
fit1 <- lm(log(data$BW)~log(data$SL),
           data = subset(data,Abbreviation == "MKZ"))
summary(fit1)
#估计参数
coef(fit1)
#估计参数的置信区间，可调用confit()函数
confint(fit1)
#对模型进行方差分析，可调用anova()函数
anova(fit1)
#采用plot()函数进行回归诊断分析
mkz <- par(mfcol = c(2,3))
plot(fit1,which = 6, col = "red")#which=1:6
par(mkz)

#回归诊断的总括
influence.measures(fit1)
library(car)
influencePlot(fit1)

#

#
##------------------------------------
##(1)残差与拟合值图，用于分析残差分布是否存在趋势，以检验方差是否为常数或因变量与自变量的关系假设是否正确；
##(2)Normal Q-Q 绘制，用于分析残差是否服从正太分布，是否存在重尾，轻尾，正偏，负篇等情况；
##(3)Scale-Location 为标准化残差平方根的分布用于分析残差的分布趋势及数据对回归的影响，是否存在异常值。一般把标准化残差的绝对值大于或等于2的观测点认为是可疑点，大于或等于3的观测点认为是异常值。对可疑点和异常值可用影响分析influence()函数进一步分析；
##(4)Cook 距离(cook.distance())用于显示对回归有重要影响的点，一般|Di|>4/n时，该点被认为是强影响点，其中n为样本容量
##-------------------------------------
weight <- predict(fit1, se.fit = T,
                  interval = "confidence",level = 0.95)
weight
#R中分别用residuals()，rstandard()h函数与rstudent()函数来计算残差，标准残差和学生残差
plot(residuals(fit1)~weight$fit[,1])
rstandard(fit1)
rstudent(fit1)
plot(cooks.distance(fit1))

#------------------------------------------------------
#------------------------------------------------------
#最大似然法估计a，b参数
Rdata <- read.csv("./data/Schizothorax nukiangensis.csv")
length <- log(Rdata$SL)
weight <- log(Rdata$BW)
#编写负对数似然函数，并赋初始值(可通过默认参数获得)
MinusLogL <- function(Sigma = log(0.2), alpha = 2.0, beta = 2.0)
{
  n <- length(weight)
  Sigma2 = 2*exp(Sigma)^2
  n*Sigma + 0.5*n*log(2*pi) +
    sum((weight - (alpha + beta*length)) *
          (weight - (alpha + beta*length)))/Sigma2
}
#利用其密度函数，可以简化MinusLogL()函数
MinusLogL <- function(Sigma = log(0.2), alpha = -5.0, beta = 3.0)
{
  -sum(stats::dnorm(weight,mean = alpha + beta*length,
                    sd = exp(Sigma), log = TRUE))
}
#利用stats4软件包中的mle()函数，进行最大似然法拟合
library("stats4")
MLEFit <- mle(MinusLogL)
#获取估计参数及标准误，最大似然法与最小二乘法估计值一致，但标准误存在较小误差
summary(MLEFit)
#利用logLik()函数获取对数似然值
logLik(MLEFit)
#利用vcov()函数获取参数方程的协方差矩阵
vcov(MLEFit)
#profile()函数，能提供参数概率剖面图
oldp <- par(mfcol = c(3, 1))
plot(profile(MLEFit), absVal = FALSE)
par(oldp)
#获取最大似然估计参数值
sigma <- coef(MLEFit)[1]
alpha <- coef(MLEFit)[2]
beta <- coef(MLEFit)[3]
#beta参数的取值范围
beta2 <- seq(2.7, 3.6, by = 0.01)
N <- length(beta2)
LL <- numeric(N)
#固定beta参数，对其他参数进行最大似然估计，并获得相应的似然值
for (i in 1:N)
{
  MLEFit1 <- mle(MinusLogL, fix = list(beta = beta2[i]))
  #固定beta参数，对其他参数进行最大似然估计
  LL[i] <- logLik(MLEFit1)
}
#获取最大似然估计时的似然值
MaxLL <- MinusLogL(sigma, alpha, beta)
#绘制似然值随beta变化的曲线
plot(beta2, LL)
#获取最小似然值
MinLL <- logLik(MLEFit) - qchisq(1-0.05, 1)/2
#获得beta的置信区间，直线与曲线的交点为其置信区间的上下边界
abline(h = MinLL, col = "green")
#参数的置信区间也可以利用confint()函数直接获取
confint(MLEFit)
betaLU <- confint(MLEFit)[3,]
abline(v =  betaLU[1], col = "red")
abline(v =  betaLU[2], col = "red")
######---------------------------------------------

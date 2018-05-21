
x<-runif(50)
y<-3+2*rnorm(50,mean = 0,sd=2)
data<-data.frame(y,x)
library(R2WinBUGS)
library(arm)
library(BRugs)
N<-length(x)
data<-list("N","y","x")
inits<-function()
{
  list(alpha=0,beta=0,tau.y=1)
}
parameters<- c("alpha","beta","tau.y")
regression.sim<-bugs(data,inits,parameters,model.file = "D:/Documents/R data/guanwenjiang/example.bug",n.chains=1,n.iter=5000,n.burnin=1000,n.thin=1,debug=T,program="WinBUGS")

#geweke.diag(regression.sim)
#heidel.diag(regression.sim)
#raftery.diag(as.mcmc(regression.sim))
plot(regression.sim)
print(regression.sim)
summary(regression.sim)
library(mcmcplots)
mcmcplot(regression.sim,dir = getwd())
attach.bugs(regression.sim)
plot(alpha,type = "l")
plot(beta,type = "l")
plot(tau.y,type = "l")


plot(density(alpha))
abline(v=3,lty=2)
quantile(alpha,probs = c(0.025,0.975))

plot(density(beta))
abline(v=2,lty=2)
quantile(beta,probs = c(0.025,0.975))

plot(density(sqrt(1/tau.y)))
abline(v=2,lty=2)
quantile(sqrt(1/tau.y),probs = c(0.025,0.975))

plot(alpha~beta)

detach.bugs()

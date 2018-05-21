###Setting your working directory
################# setwd("………..")
##Loading available K and Loo estimates from studies included in Fishbase
FishBaseData <- read.table("./data/FishBaseMullusBarbatus.txt",header = T)
FishBaseData
##Entering mean length at age data from Sieli et al. 2011
Age <- c(1,2,3,4,5,6,7)
Length <- c(10.93,12.1,15.23,17.24,18.36,20.18,19.77)

###### METHOD 1##############
##Fitting the von Bertalanffy growth model with least square non linear regression
##### ##Setting the von Bertalanffy growth model equation
vonBertalanffy <- Length~Loo*(1-exp(-K*(Age-To)))
###Setting initial values for the 3 parameters of the model
InitLoo <- max(Length)
InitK <- mean(FishBaseData$K)
InitTo <- 0
##Running the non linear regression model
NLRBertalanffy <- nls(vonBertalanffy,
                      data = list(Age = Age,Length = Length),
                      start = list(Loo = InitLoo,K = InitK,To = InitTo))

######METHOD 2 ######
######Fitting the von Bertalanffy growth model with Bayesian inference using previous knowledge###
#### Fitting normal distributions to available FishBase data for using them as priors. Log-normal distributions could also be fitted, but let’s stick with the simpler normal distribution for better comprehension.
library(MASS)
PriorLoo <- fitdistr(FishBaseData$Loo,"normal")
PriorK <- fitdistr(FishBaseData$K,"normal")
###Printing the fitted distributions PriorLoo PriorK
library(R2jags)
runif(1)
#######Writing the Model in Jags notation
Model = "model { for( i in 1 :Nobs)
                  { Length[i] ~ dnorm(mu[i],tau)
                   mu[i]<-Loo*(1-exp(-K*(Age[i]-To)))
                  }
##Setting our informative priors for K and Loo as estimated above
    Loo ~ dnorm(24.9,0.055)
    K ~ dnorm (0.368,22.3)
########Using an uniformative prior for To and variability at the individual level
    tau ~ dgamma(0.001,0.001)
    To ~ dnorm (0,0.0001)
      }"
#Saving the model file in the working directory
cat(Model, file = "Model.bug")
####Defining data
DataJ <- list(Age = Age,Length = Length,Nobs = length(Age))
##Running the model in Jags
Jags <- jags(model.file = "Model.bug", working.directory = NULL,
             data = DataJ, parameters.to.save = c("K","Loo","To"),
             n.chains = 3,n.thin = 10, n.iter = 20000, n.burnin = 10000)
#####Getting results from the jags model
Jags$BUGSoutput$mean$K
#####Getting results from the non linear regression
summary(NLRBertalanffy)
#####Plotting our results
K_NL <- coef(NLRBertalanffy)[2]
K_Jags <- Jags$BUGSoutput$mean$K
Loo_NL <- coef(NLRBertalanffy)[1]
Loo_Jags <- Jags$BUGSoutput$mean$Loo
sd_K_NL <- summary(NLRBertalanffy)$coefficients[2,2]
sd_K_Jags <- Jags$BUGSoutput$sd$K
sd_Loo_NL <- summary(NLRBertalanffy)$coefficients[1,2]
sd_Loo_Jags <- Jags$BUGSoutput$sd$Loo
a1 <- Loo_NL - 3*sd_Loo_NL
a2 <- Loo_NL + 3*sd_Loo_NL
b1 <- K_NL - 3*sd_K_NL
b2 <- K_NL + 3*sd_K_NL
xK <- seq(b1,b2,length.out = 200)
xLoo <- seq(a1,a2,length.out = 200)
par(mfrow = c(1,2))
######Plotting posterior distribution of K from Jags
plot(xK,dnorm(xK,mean = K_Jags,sd = sd_K_Jags)/100,
     type = "l",lwd = 1,xlab = "K",ylab = "")
######Adding distribution of K from non linear regression
lines(xK,dnorm(xK,K_NL,sd_K_NL)/100,lty = 3, col = "red")
######Plotting posterior distribution of Loo from Jags
plot(xLoo,dnorm(xLoo,Loo_Jags,sd_Loo_Jags),type = "l",lwd = 1,xlab = "Loo",ylab = "")
######Adding distribution of K from non linear regression
lines(xLoo,dnorm(xLoo,Loo_NL,sd_Loo_NL),lty = 3, col = "red")
jpeg("PlotNLvsBayesianGrowth.jpeg")
dev.off()

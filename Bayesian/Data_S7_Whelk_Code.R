### West 1988 data

### clear workspace

rm(list = ls())

### load packages

library(tidyr)
library(dplyr)
library(rjags)
library(loo)

### set working directory

setwd("~/whelkdata")

### load data

data  <- read.csv("data/Whelk_Data_S6.csv")

### convert data to wide format

unique(data$Prey)

data2 <- data %>% group_by(Site, Size.Class, Snail.ID, Prey) %>% summarise(Number = n())


data2 <- spread(data2, Prey, Number, fill = 0)

### Prepare for analysis

### convert to data matrix for analysis in JAGS

data3 <- data.matrix(data2)

### load model specification

sink("multinom.dir.hier.txt")
cat("
    data{
    # calculate the number of prey that each individual was observed feeding on
    for (i in 1:N) {
    nisk[i]  <- sum(yisk[i,1:J])
    }
    }
    model {
    # model description for individuals nested within size class and site
    for (i in 1:N) {
    yisk[i, 1:J] ~ dmulti(pisk[i,1:J], nisk[i]) ### observations of individual diets are distributed multinomially with probability vector pisk and number of obs. nisk
    pisk[i, 1:J] ~ ddirch(alphask[sizeclass[i],site[i], 1:J]) ### pisk is distributed dirichlet with vector alphask
    }
    for (s in 1:S) {
    for (k in 1:K) {
    for (j in 1:J) {
    alphask[s,k,j]  <- qsk[s,k,j] * wsk[s,k] + 0.05 ### each component of alphask is determined by the mean diet for the size class and site, qsk, and the concentration parameter wsk
    }
    }
    }
    for (s in 1:S) {
    for (k in 1:K) {
    qsk[s,k, 1:J] ~ddirch(alphak[k, 1:J]) ### qsk is distributed dirichlet alphak
    }
    }
    for (k in 1:K) {
    for (j in 1:J) {
    alphak[k,j]  <- qk[k,j] * wk[k] + 0.05 ### each component of alphak is determined by the mean diet for the site, qk, and the concentration parameter wk
    }
    }
    for (k in 1:K) {
    qk[k, 1:J] ~ ddirch(alphapop[1:J]) ### qk is distributed dirichlet alphapop
    }
    for (j in 1:J) {
    alphapop[j]  <- q[j] * w + 0.05 ### each component of alphapop is determined by the overall mean population diet, q, and the concentration parameter w
    }
    q[1:J] ~ ddirch(alpha[]) ### alpha represents the prior for q. Here it is a vector of 1's or a uniform prior
    for (j in 1:J){
    alpha[j]  <- 1
    }
    for (s in 1:S) {
    for (k in 1:K) {
    wsk[s,k] ~ dunif(.1, 30) ### define the prior for wsk
    }
    }
    for (k in 1:K) {
    wk[k] ~ dunif(.1, 30) ### define the prior for wk
    }
    w ~ dunif(0.1, 30) ### define the prior for w
    # calculations to estimate individual PSi and the population mean PS (IS)
    for (i in 1:N){
    for(j in 1:J){
    diff.prop[i,j]  <- abs(pisk[i, j] - qsk[sizeclass[i],site[i],j])
    }
    PS[i]  <- 1 - 0.5 * sum(diff.prop[i,1:J])
    }
    mean.PS  <- mean(PS[])
    # calculate log likelihood for each data point at each iteration for calculating waic
    for (i in 1:N){
      log_lik[i] <- logdensity.multi(yisk[i,], pisk[i,1:J], nisk[i])
    }
    }
    ", fill = TRUE)
sink()

### set up the data to be entered into the JAGS model (define the data variables to be supplied to the model)

dataJags  <- list(yisk = data3[,4:22], N = as.numeric(length(data3[,1])), J = 19, K = 2, S = 3,
              sizeclass = data3[,2], site = data3[,1])

### calculations to give initial values for q and w, other parameters are given randomly generated initial values

RowTotal  <- apply(data3[,4:22], 1, sum)

Props  <- as.numeric(apply(data3[,4:22], 2, function(x){sum(x)/sum(RowTotal)}))
inits  <- function(){list(q = Props, w = (1))}

### initialize the jags model with three chains and a burn in period of 10,000 iterations

jags <- jags.model("multinom.dir.hier.txt", data = dataJags,
                   n.chains = 3, n.adapt = 10000, inits = inits )

### supply the parameters that the model should estimate

params  <- c("pisk", "wsk", "wk", "qsk", "qk", "w", "q", "PS", "mean.PS", "log_lik")

### sample 2,000 iterations for each chain with a thinning of 100 iterations

samples  <- jags.samples(jags, variable.names = params, n.iter = 2000,
                         n.thin = 100)

### create matrix of log-likelihoods to calculate waic
### matrix is SxN with S being the number of simulations and N being the number of observations in the model

log_lik_mat1 <- t(cbind(samples$log_lik[,,1], samples$log_lik[,,2], samples$log_lik[,,3]))

### use package 'loo' to calculate WAIC

SiteSizeWAIC <- waic(log_lik_mat1)

### Calculate DIC if preferred over WAIC

WhelkSiteSizeDIC <- dic.samples(jags, 1000, n.thin = 100, type = "pD")


### analyze data with site dropped
sink("multinom.dir.hier.txt")
cat("
    data{
    for (i in 1:N) {
    nis[i]  <- sum(yis[i,1:J]) # calculate sample size for each individual
    }
    }
    model {
    for (i in 1:N) {
    yis[i, 1:J] ~ dmulti(pis[i,1:J], nis[i]) # individual observations used to estimate preferences pis
    pis[i, 1:J] ~ ddirch(alphas[sizeclass[i], 1:J])  # each individual's preferences are drawn from a dirichlet distribution with vector alphas
    }
    for (s in 1:S) {
    for (j in 1:J) {
    alphas[s,j]  <- qs[s,j] * ws[s] + 0.05 # components of alphas are defined by size class level preference qs and concentration parameter ws
    }
    }
    for (s in 1:S) {
    qs[s, 1:J] ~ ddirch(alphapop[1:J]) # size class level preferences are distributed Dirichlet with vector alphapop
    }
    for (j in 1:J) {
    alphapop[j]  <- q[j] * w + 0.05 # components of alphapop are defined by the population preference q and the concentration parameter w
    }
    q[1:J] ~ ddirch(alpha[]) # define the prior for q
    for (j in 1:J){
    alpha[j]  <- 1 # a vector of 1's leads to a uniform distribution
    }
    for (s in 1:S) {
    ws[s] ~ dunif(0, 30) # set prior for ws's
    }
    w ~ dunif(0, 30) # set prior for w
    for (i in 1:N){
    for(j in 1:J){
    diff.prop[i,j]  <- abs(pis[i, j] - qs[sizeclass[i],j]) # calculate PSi
    }
    PS[i]  <- 1 - 0.5 * sum(diff.prop[i,1:J])
    }
    mean.PS  <- mean(PS[]) # calculate mean PSi for each iteration
    # calculate log-likelihood for loo calculation
    for (i in 1:N) {
      log_lik[i] <- logdensity.multi(yis[i, ], pis[i,1:J], nis[i])
    }
    }
    ", fill = TRUE)
sink()


### set up the data to be analyzed by JAGS

dataJags  <- list(yis = data3[,4:22], N = as.numeric(length(data3[,1])), J = 19, S = 3,
                  sizeclass = data3[,2])

### calculations to give initial values for q and w, other parameters are given randomly generated initial values

RowTotal  <- apply(data3[,4:22], 1, sum)

Props  <- as.numeric(apply(data3[,4:22], 2, function(x){sum(x)/sum(RowTotal)}))
inits  <- function(){list(q = Props, w = (1))}

### initialize the jags model with three chains and a burn in period of 10,000 iterations


jags2 <- jags.model("multinom.dir.hier.txt", data = dataJags,
                   n.chains = 3, n.adapt = 10000, inits = inits )

### supply the parameters that the model should estimate

params  <- c("pis", "ws", "w", "qs", "q", "PS", "mean.PS", "log_lik")

### sample 2,000 iterations for each chain with a thinning of 100 iterations

samples2 <- jags.samples(jags2, variable.names = params, n.iter = 2000,
                        n.thin = 100)

### create matrix of log-likelihoods to calculate waic
### matrix is SxN with S being the number of simulations and N being the number of observations in the model

log_lik_mat2 <- t(cbind(samples2$log_lik[,,1], samples2$log_lik[,,2], samples2$log_lik[,,3]))

### use package 'loo' to calculate WAIC

SizeWAIC <- waic(log_lik_mat2)

### Calculate DIC values if preferred over WAIC

WhelkSizeDIC <- dic.samples(jags, 1000, n.thin = 100, type = "pD")

### Model without size class or site

sink("multinom.dir.hier.txt")
cat("
    data{
    for (i in 1:N) {
    ni[i]  <- sum(yi[i,1:J]) ### calculate the number of prey observed for each individual
    }
    }
    model {
    for (i in 1:N) {
    yi[i, 1:J] ~ dmulti(pi[i,1:J], ni[i]) # individual observations used to estimate preferences pi
    pi[i, 1:J] ~ ddirch(alphapop[1:J]) # individual's preferences are drawn from a dirichlet distribution
    }
    for (j in 1:J) {
    alphapop[j]  <- q[j] * w + 0.05 # the components of alphapop are determined by the population mean diet, q, and the concentration parameter w
    }
    q[1:J] ~ ddirch(alpha[]) ### define the uniform prior for q
    for (j in 1:J){
    alpha[j]  <- 1
    }
    w ~ dunif(0.1, 30) ### define the prior for w
    ### calculate PSi
    for (i in 1:N){
    for(j in 1:J){
    diff.prop[i,j]  <- abs(pi[i, j] - q[j])
    }
    PS[i]  <- 1 - 0.5 * sum(diff.prop[i,1:J])
    }
    mean.PS  <- mean(PS[])
    # calculate log likelihood for each data point at each iteration for calculating waic
    for (i in 1:N){
      log_lik[i] <- logdensity.multi(yi[i, ], pi[i, ], ni[i])
    }
    }
    ", fill = TRUE)
sink()

### set up the data to be analyzed by JAGS

dataJags  <- list(yi = data3[,4:22], N = as.numeric(length(data3[,1])), J = 19)

### calculations to give initial values for q and w, other parameters are given randomly generated initial values

RowTotal  <- apply(data3[,4:22], 1, sum)

Props  <- as.numeric(apply(data3[,4:22], 2, function(x){sum(x)/sum(RowTotal)}))
inits  <- function(){list(q = Props, w = (1))}

### initialize the jags model with three chains and a burn in period of 10,000 iterations

jags3 <- jags.model("multinom.dir.hier.txt", data = dataJags,
                   n.chains = 3, n.adapt = 10000, inits = inits )

### supply the parameters that the model should estimate

params  <- c("pi", "q", "w", "PS", "mean.PS", "log_lik")

### sample 2,000 iterations for each chain with a thinning of 100 iterations

samples3  <- jags.samples(jags3, variable.names = params, n.iter = 2000,
                         n.thin = 100)

### create matrix of log-likelihoods to calculate waic
### matrix is SxN with S being the number of simulations and N being the number of observations in the model

log_lik_mat3 <- t(cbind(samples3$log_lik[,,1], samples3$log_lik[,,2], samples3$log_lik[,,3]))

### use package 'loo' to calculate WAIC

BaseWAIC <- waic(log_lik_mat3)

### get DIC values if preferred over WAIC

WhelkDIC <- dic.samples(jags3, 1000, n.thin = 100, type = "pD")


### Analysis of diet specialization with heterogeneity in sample sizes among individuals

### Load packages, set working directory, set BUGS directory

# Clear workspace
rm(list = ls())
### load required packages
library(R2OpenBUGS)
library(RInSp)
library(MCMCpack)
### set working directory
#setwd("~/IndSpecEst")

### set BUGS directory (location of the BUGS.exe file on the machine)
bd  <- "F:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

### write model and save to working directory
sink("multinom.dir.hier.txt")
cat("
    model {
    for (i in 1:N) {
    y [i, 1:J] ~ dmulti(pi[i, 1:J], n[i])   # multinomial distribution of response vectors
    n[i] <- sum(y[i, ]) # calculate sample size per individual
    pi[i, 1:J] ~ ddirch(alpha[]) # individual preferences are distributed Dirichlet(alpha)
    }
    for (j in 1:J){
    alpha[j] <- q[j] * w # estimate the parameters that make up alpha
    }
    q[1:J]  ~ ddirch(u[1:J]) # prior for q
    for (j in 1:J){
    u[j]  <- 1
    }
    w ~ dunif(0, 30) # prior for w
    for (i in 1:N){
    for(j in 1:J){
    diff.prop[i,j]  <- abs(pi[i, j] - q[j])
    }
    PS[i]  <- 1 - 0.5 * sum(diff.prop[i,1:J]) # calculate PSi
    }
    mean.PS  <- mean(PS[]) # calculate IS
    }
    ", fill = TRUE)
sink()

### heterogeneity in observations per individual
### Uniform
# Prepare array to save results
datahetunif  <- array(dim = c(100, 11, 500))
name  <- list(NULL, c("p1", "p2", "p3", "p4", "truePSi", "BayesPSi", "BayesVarPSi",
                      "Bayes95low", "Bayes95high", "RinSpPSi", "RinSpVarPSi"), NULL)
dimnames(datahetunif)  <- name
whet  <- vector()
PopProbhet  <- matrix(nrow = 500, ncol = 4)

### Set up loop
for (i in 1:500){
  try( {
    q = rdirichlet(1, c(1,1,1,1)) # overall population preference
    PopProbhet[i,]  <- q
    scale  <- runif(1, 1, 10)
    whet[i]  <- scale
    p_mat = rdirichlet(100, scale*q) # individual-level preferences
    datahetunif[,1:4,i]  <- p_mat
    X = t(apply(p_mat,1,function(x) rmultinom(1,ceiling(100*rbeta(1,1,1)),x))) # individual-level observations

    ### get true PS
    true.PS  <- c()
    prop.diff  <- matrix(nrow = 100, ncol = 4)
    for (t in 1:100) {
      for (j in 1:4) {
        prop.diff[t,j]  <- p_mat[t, j] - q[j]
      }
      true.PS[t]  <- 1 - 0.5 * sum(abs(prop.diff[t,]))
    }
    datahetunif[,5,i]  <- true.PS

    ### prepare data for BUGS analysis
    bug.data  <- list(y = X, N = as.numeric(length(X[,1])), J = as.numeric(length(X[1,])))

    ### parameters to estimate
    params  <- c("q", "pi", "w", "PS", "mean.PS")

    ### initial values function for the Markov chains
    inits  <- function(){list(q = c(.25, .25, .25, .25), w = (1)) }
    ### MCMC settings
    # number of Markov Chains
    nc  <- 3
    # number of samples from Markov Chains
    ni  <- 1000
    # number of iterations for Burn-in period
    nb  <- 200
    ### Thinning - Sample every 'nt' iterations
    nt  <- 5

    ### start Gibbs sampler
    out  <- bugs(bug.data, inits, params, model.file = "multinom.dir.hier.txt", n.thin = nt, n.chains = nc,
                 n.burnin = nb, n.iter = ni, working.directory = getwd(), OpenBUGS.pgm = bd, debug = FALSE)

    ### Run analysis using RInSp
    RInSpdata  <- import.RInSp(X)
    PSInSpdata  <- PSicalc(RInSpdata, pop.diet = "average", replicates = 5000)
    ### save data from analyses
    datahetunif[,6,i]  <- out$median$PS
    datahetunif[,7,i]  <- (out$sd$PS)^2
    datahetunif[,8,i]  <- apply(out$sims.list$PS, 2, function(x) quantile(x,probs = 0.025))
    datahetunif[,9,i]  <- apply(out$sims.list$PS, 2, function(x) quantile(x,probs = 0.975))
    datahetunif[,10,i]  <- PSInSpdata$PSi
    datahetunif[,11,i]  <- PSInSpdata$VarPSi

  }
  )
}

### heterogeneity in observations per individual
### weighed to low sample size - Beta(0.5, 1)
datahetlow  <- array(dim = c(100, 11, 500))
name  <- list(NULL, c("p1", "p2", "p3", "p4", "truePSi", "BayesPSi", "BayesVarPSi",
                      "Bayes95low", "Bayes95high", "RinSpPSi", "RinSpVarPSi"), NULL)
dimnames(datahetlow)  <- name
whet  <- vector()
PopProbhetlow  <- matrix(nrow = 500, ncol = 4)

### Set up loop
for (i in 1:500){
  try( {
    q = rdirichlet(1, c(1,1,1,1))
    PopProbhetlow[i,]  <- q
    scale  <- runif(1, 1, 10)
    whet[i]  <- scale
    p_mat = rdirichlet(100, scale*q)
    datahetlow[,1:4,i]  <- p_mat
    X = t(apply(p_mat,1,function(x) rmultinom(1,ceiling(100*rbeta(1,.5,1)),x)))

    ### get true PS
    true.PS  <- c()
    prop.diff  <- matrix(nrow = 100, ncol = 4)
    for (t in 1:100) {
      for (j in 1:4) {
        prop.diff[t,j]  <- p_mat[t, j] - q[j]
      }
      true.PS[t]  <- 1 - 0.5 * sum(abs(prop.diff[t,]))
    }
    datahetlow[,5,i]  <- true.PS

    # prepare data for BUGS analysis
    bug.data  <- list(y = X, N = as.numeric(length(X[,1])), J = as.numeric(length(X[1,])))

    ### parameters to estimate
    params  <- c("q", "pi", "w", "PS", "mean.PS")

    ### inits function
    inits  <- function(){list(q = c(.25, .25, .25, .25), w = (1)) }
    ### MCMC settings
    nc  <- 3
    ni  <- 1000
    nb  <- 200
    nt  <- 5

    ### start Gibbs sampler
    out  <- bugs(bug.data, inits, params, model.file = "multinom.dir.hier.txt", n.thin = nt, n.chains = nc,
                 n.burnin = nb, n.iter = ni, working.directory = getwd(), OpenBUGS.pgm = bd, debug = FALSE)

    ### Run analysis using RInSp
    RInSpdata  <- import.RInSp(X)
    PSInSpdata  <- PSicalc(RInSpdata, pop.diet = "average", replicates = 5000)
    ### save data from analyses
    datahetlow[,6,i]  <- out$median$PS
    datahetlow[,7,i]  <- (out$sd$PS)^2
    datahetlow[,8,i]  <- apply(out$sims.list$PS, 2, function(x) quantile(x,probs = 0.025))
    datahetlow[,9,i]  <- apply(out$sims.list$PS, 2, function(x) quantile(x,probs = 0.975))
    datahetlow[,10,i]  <- PSInSpdata$PSi
    datahetlow[,11,i]  <- PSInSpdata$VarPSi

  }
  )
}


### heterogeneity in observations per individual
### weighed to high sample size - Beta(1, 0.5)

datahethigh  <- array(dim = c(100, 11, 500))
name  <- list(NULL, c("p1", "p2", "p3", "p4", "truePSi", "BayesPSi", "BayesVarPSi",
                      "Bayes95low", "Bayes95high", "RinSpPSi", "RinSpVarPSi"), NULL)
dimnames(datahethigh)  <- name
whet  <- vector()
PopProbhethigh  <- matrix(nrow = 500, ncol = 4)

### Set up loop
for (i in 1:500){
  try( {
    q = rdirichlet(1, c(1,1,1,1))
    PopProbhethigh[i,]  <- q
    scale  <- runif(1, 1, 10)
    whet[i]  <- scale
    p_mat = rdirichlet(100, scale*q)
    datahethigh[,1:4,i]  <- p_mat
    X = t(apply(p_mat,1,function(x) rmultinom(1,ceiling(100*rbeta(1,1,.5)),x)))

    ### get true PS
    true.PS  <- c()
    prop.diff  <- matrix(nrow = 100, ncol = 4)
    for (t in 1:100) {
      for (j in 1:4) {
        prop.diff[t,j]  <- p_mat[t, j] - q[j]
      }
      true.PS[t]  <- 1 - 0.5 * sum(abs(prop.diff[t,]))
    }
    datahethigh[,5,i]  <- true.PS

    # prepare data for BUGS analysis
    bug.data  <- list(y = X, N = as.numeric(length(X[,1])), J = as.numeric(length(X[1,])))

    ### parameters to estimate
    params  <- c("q", "pi", "w", "PS", "mean.PS")

    ### inits function
    inits  <- function(){list(q = c(.25, .25, .25, .25), w = (1)) }
    ### MCMC settings
    nc  <- 3
    ni  <- 1000
    nb  <- 200
    nt  <- 5

    ### start Gibbs sampler
    out  <- bugs(bug.data, inits, params, model.file = "multinom.dir.hier.txt", n.thin = nt, n.chains = nc,
                 n.burnin = nb, n.iter = ni, working.directory = getwd(), OpenBUGS.pgm = bd, debug = FALSE)

    ### Run analysis using RInSp
    RInSpdata  <- import.RInSp(X)
    PSInSpdata  <- PSicalc(RInSpdata, pop.diet = "average", replicates = 5000)
    ### save data from analyses
    datahethigh[,6,i]  <- out$median$PS
    datahethigh[,7,i]  <- (out$sd$PS)^2
    datahethigh[,8,i]  <- apply(out$sims.list$PS, 2, function(x) quantile(x,probs = 0.025))
    datahethigh[,9,i]  <- apply(out$sims.list$PS, 2, function(x) quantile(x,probs = 0.975))
    datahethigh[,10,i]  <- PSInSpdata$PSi
    datahethigh[,11,i]  <- PSInSpdata$VarPSi

  }
  )
}

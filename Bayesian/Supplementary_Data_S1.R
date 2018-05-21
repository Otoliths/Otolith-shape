### Estimating PSi and its variance
### Load packages, set working directory, set BUGS directory
rm(list = ls())
### load required packages
library(R2OpenBUGS)
library(RInSp)
library(MCMCpack)
### set working directory
setwd("~/IndSpecEst")
### set BUGS directory (location of the BUGS.exe file on your machine)
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

### 10 observations per individual
### create array to store data
data10  <- array(dim = c(100, 11, 500))
name  <- list(NULL, c("p1", "p2", "p3", "p4", "truePSi", "BayesPSi", "BayesVarPSi",
                      "Bayes95low", "Bayes95high", "RinSpPSi", "RinSpVarPSi"), NULL)
dimnames(data10)  <- name
w10  <- vector()
PopProb10  <- matrix(nrow = 500, ncol = 4)

### Set up loop
for (i in 1:500){
  try( {
    q = rdirichlet(1, c(1,1,1,1)) # overall preferences
    PopProb10[i,]  <- q
    scale  <- runif(1, 1, 10) # concentration parameter for the population Dirichlet distribution
    w10[i]  <- scale
    p_mat = rdirichlet(100, scale*q) # individual-level preferences
    data10[,1:4,i]  <- p_mat
    X = t(apply(p_mat,1,function(x) rmultinom(1,10,x))) # individual-level observations

    ### get true PS
    true.PS  <- c()
    prop.diff  <- matrix(nrow = 100, ncol = 4)
    for (t in 1:100) {
      for (j in 1:4) {
        prop.diff[t,j]  <- p_mat[t, j] - q[j]
      }
      true.PS[t]  <- 1 - 0.5 * sum(abs(prop.diff[t,]))
    }
    data10[,5,i]  <- true.PS

    #prepare data for BUGS model
    bug.data  <- list(y = X, N = as.numeric(length(X[,1])), J = as.numeric(length(X[1,])))

    ### parameters to estimate
    params  <- c("q", "pi", "w", "PS", "mean.PS")

    ### initial values function for Markov Chains
    inits  <- function(){list(q = c(.25, .25, .25, .25), w = (1)) }
    ### MCMC settings
    ### Number of chains
    nc  <- 3
    # Number of sampling iterations
    ni  <- 1000
    # Number of iterations in the Burn-in period
    nb  <- 200
    ### Thinning - Sample every 'nt' iterations
    nt  <- 5

    ### start Gibbs sampler
    out  <- bugs(bug.data, inits, params, model.file = "multinom.dir.hier.txt", n.thin = nt, n.chains = nc,
                 n.burnin = nb, n.iter = ni, working.directory = getwd(), OpenBUGS.pgm = bd, debug = FALSE)

    ### Run analysis using RInSp
    # Create an RInSp data object
    RInSpdata  <- import.RInSp(X)
    # Calculate PSi
    PSInSpdata  <- PSicalc(RInSpdata, pop.diet = "average", replicates = 5000)
    ### save data from analyses
    data10[,6,i]  <- out$median$PS
    data10[,7,i]  <- (out$sd$PS)^2
    data10[,8,i]  <- apply(out$sims.list$PS, 2, function(x) quantile(x,probs = 0.025))
    data10[,9,i]  <- apply(out$sims.list$PS, 2, function(x) quantile(x,probs = 0.975))
    data10[,10,i]  <- PSInSpdata$PSi
    data10[,11,i]  <- PSInSpdata$VarPSi

  }
  )
}

### To change the number of observations per individual, the second term in the rmultinom function in the line,
# X = t(apply(p_mat,1,function(x) rmultinom(1,10,x))) # individual-level observations
# is changed to the number of desired samples per individual.
# For example, X = t(apply(p_mat,1,function(x) rmultinom(1,100,x)))
# will simulate data with 100 observations per individual. The rest of the code can be used
# unaltered.

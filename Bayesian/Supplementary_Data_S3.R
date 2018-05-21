### Population level estimation of individual diet specialization

### Load packages, set working directory, set BUGS directory

### clear workspace
rm(list = ls())
### load required packages
library(R2OpenBUGS)
library(RInSp)
library(MCMCpack)
### set working directory
setwd("~/IndSpecEst")
### set BUGS directory (the location of the BUGS.exe file on your machine)
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


### set up vectors to save simulation results
wactual  <- c()
west  <- c()
qactual  <- matrix(nrow = 1000, ncol = 4)
WICTNWactual  <- c()
ISactual  <- c()
WICTNWest <- c()
ISest <- c()

for (i in 1:1000){
  try( {
    q = rdirichlet(1, c(1,1,1,1)) # overall proportions
    qactual[i,]  <- q
    scale  <- runif(1, .5, 10)
    wactual[i]  <- scale
    p_mat = rdirichlet(100, scale*q) # individual-level proportions
    X = t(apply(p_mat,1,function(x) rmultinom(1,round(runif(1, 3, 50)),x)))

    ### Bayesian estimation
    bug.data  <- list(y = X, N = as.numeric(length(X[,1])), J = as.numeric(length(X[1,])))

    ### parameters to estimate
    params  <- c("w")

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

    ### Run PS analysis using RInSp
    RInSpdata  <- import.RInSp(X)
    PSInSpdata  <- PSicalc(RInSpdata, pop.diet = "average", replicates = 5000)

    ### Run WIC/TNW analysis using RInSp

    WTdMCdat  <- WTdMC(RInSpdata, pop.diet = "average", replicates = 5000)

    ### save data from analyses
    WICTNWest[i]  <- WTdMCdat$WonT
    west[i]  <- out$median$w
    ISest[i]  <- PSInSpdata$IS
    PSi <- c()
    for(j in 1:100){PSi[j] <- 1 - 0.5 * sum(abs(p_mat[j,] - q))}
    ISactual[i] <- mean(PSi)
    WICstep1 <- c()
    indtot <- apply(X, 1, sum)
    poptot <- sum(indtot)
    indprop <- indtot/poptot
    for (j in 1:100) {
      WICstep1[j] <- indprop[j] * -(sum(p_mat[j,] * log(p_mat[j,])))
    }
    WIC <- sum(WICstep1)
    TNW <- -sum(q * log(q))
    WICTNWactual[i] <- WIC/TNW
  })
}

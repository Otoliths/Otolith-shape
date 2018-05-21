### SIMULATING INDIVIDUALS THEN FITTING DATA to the SIMULATION

rm(list=ls())

palette("default")

####### WRITE FUNCTIONS TO CALCULATE THE MEAN AND VARIANCE OF THE POPULATION.

skew	<- function(x){
  m3	<- sum((x-mean(x))^3)/length(x)
  s3	<- sqrt(var(x))^3
  m3/s3}

EXPECT	<-	function(X0,ALPHA,BETA,PSI,G){

  if(max(G)>0){
    AGE	<-	max(which(G>0))

    ## X0.term
    FIRST.term	<- exp(log(X0) + ALPHA*log(BETA) - (ALPHA) * log(AGE+BETA))
    ## Gamma.terms
    term.1		<-	log(G[1:AGE]) +
      ALPHA * log(BETA) -
      (ALPHA + PSI -1) * log(AGE - seq(1:AGE) + BETA) +
      lgamma(ALPHA + PSI-1) -
      lgamma(ALPHA)
    term.2		<-	log(G[1:AGE]) +
      ALPHA * log(BETA) -
      (ALPHA + PSI -1) * log(AGE - seq(1:AGE) + 1 + BETA) +
      lgamma(ALPHA + PSI-1) -
      lgamma(ALPHA)

    sum(exp(term.1) - exp(term.2)) + FIRST.term
  }
}

VAR		<-	function(X0,DELTA,PSI,ALPHA,BETA,G){
  if(max(G)>0){
    AGE	<-	max(which(G>0))

    J	<-	seq(1,AGE,by=1)
    Z	<-	G[1:AGE]

    term.1	<-	exp( 2* log(X0) + ALPHA * log(BETA) - ALPHA* log(2*AGE + BETA) )

    term.2	<- sum(
      exp(
        log(2) + log(X0) + log(Z) + lgamma(ALPHA + PSI -1) - lgamma(ALPHA) +
          log( exp(ALPHA * log(BETA) - (ALPHA + PSI -1) * log(2*AGE -J + BETA)) -
                 exp(ALPHA * log(BETA)- (ALPHA + PSI -1) * log(1 + 2*AGE -J + BETA)) )
      ))

    Q.J	<-	Z
    Q.I	<-	Z

    temp.sum	<-	data.frame(matrix(0,length(Z),1))
    for(i in 1:length(Q.I)){
      I.numb =i
      fixer		<- rep(1,length(Q.I))
      fixer[i]	<- (1 + DELTA^-1)
      temp.sum[i,1]	<- sum(
        exp(
          log(Q.I[i]) + log(Q.J) + lgamma(ALPHA + 2*PSI -2) - lgamma(ALPHA) +
            log(
              exp(ALPHA * log(BETA) - (ALPHA + 2*PSI -2) * log(2*AGE-I.numb-J+BETA)) -
                exp(ALPHA * log(BETA) + log(2) - (ALPHA + 2*PSI -2) * log(1+2*AGE-I.numb-J+BETA)) +
                exp(ALPHA * log(BETA) - (ALPHA + 2*PSI -2) * log(2+ 2*AGE-I.numb-J+BETA))
            )
        ) * fixer
      )
    }

    term.3	<- sum(temp.sum)

    term.1 + term.2 + term.3
  }
}

#####################################################################################
#####################################################################################

#### READ IN  SCENARIO TO RUN FROM FILE
setwd("/Users/ole/Documents/Von B Growth/INPUT SCENARIOS")

input.dat <- read.csv(file="july-25-2011-scenarios.csv",header=T)

for(XXX in 10:length(input.dat$scenario)){

  NAME	<- paste(input.dat$scenario[XXX])

  #### INPUTS
  # Process
  ######################################################################################
  mean.k		<- input.dat$mean.k[XXX]
  #sd.k		<- A.true/B.true^2

  psi.true	<- input.dat$psi[XXX] #.00000001
  psi.assume	<- input.dat$assume.psi[XXX] #.00000001

  DELTA.true	<-	1/ input.dat$cv.env[XXX]^2
  ##### K ALPHA AND K BETA
  CV 		<- input.dat$cv.k[XXX]
  A.true  <-	1 / CV^2
  B.true  =	A.true / input.dat$mean.k[XXX]
  T  		=	1

  ##### STARTING SIZE
  x0.true	=	1
  nu.true	=	0.00001
  ##### AMONG VARIATION IN Q
  L.infin 	= input.dat$L.infinity[XXX]
  gam.true	= L.infin / mean.k^(psi.true-1)		 #	c(85,90,95,85,90,95,85,90,95)

  ############## NUMBER OF INDIVIDUALS TO SIMULATE
  N	=	10000
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  ### SIMULATE K fixed for individuals, Q for each individual a random draw from
  ###	-	a normal distribution each year
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################

  DAT.1	<- NULL #data.frame(matrix(0,1,2))
  #for(j in 1:2){

  cohort = 1990
  YEARS.SIM = 12

  DAT.fixed.k		<-	data.frame(matrix(0,N,YEARS.SIM))
  DAT.fixed.k[,1]	<-	rnorm(N,mean=x0.true,sd=nu.true)

  k			<- rgamma(N,shape=A.true,rate=B.true)
  gam.all		=	rep(gam.true,YEARS.SIM-1)

  for(i in 2:YEARS.SIM){	### LOOP OVER TIME
    gamma.t	<- gam.all[i-1]
    w		<- rgamma(N,shape=DELTA.true,rate=DELTA.true)
    #	w		<- exp(rnorm(N,0 - sig.true^2/2,sig.true))
    DAT.fixed.k[,i]	<-	exp(-k)*DAT.fixed.k[,i-1] +
      (gamma.t * k^(psi.true-1)*(1-exp(-k)) ) * w
  }

  #DAT.fixed.k	<- DAT.fixed.k[,2:YEARS.SIM]
  hist(rgamma(N,shape=DELTA.true,rate=DELTA.true))

  par(mfrow=c(3,4))
  hist(DAT.fixed.k[,1])
  hist(DAT.fixed.k[,2])
  hist(DAT.fixed.k[,3])
  hist(DAT.fixed.k[,4])
  hist(DAT.fixed.k[,5])
  hist(DAT.fixed.k[,6])
  hist(DAT.fixed.k[,7])
  hist(DAT.fixed.k[,8])
  hist(DAT.fixed.k[,9])
  hist(DAT.fixed.k[,10])
  hist(DAT.fixed.k[,11])
  #	hist(DAT.fixed.k[,12])
  #	hist(DAT.fixed.k[,13])
  #	hist(DAT.fixed.k[,14])
  #	hist(DAT.fixed.k[,15])
  #	hist(DAT.fixed.k[,16])
  #	hist(DAT.fixed.k[,17])
  #	hist(DAT.fixed.k[,18])
  #	hist(DAT.fixed.k[,19])
  #	hist(DAT.fixed.k[,20])
  #
  skew(DAT.fixed.k[,1])
  skew(DAT.fixed.k[,2])
  skew(DAT.fixed.k[,3])
  skew(DAT.fixed.k[,4])
  skew(DAT.fixed.k[,5])
  skew(DAT.fixed.k[,6])
  skew(DAT.fixed.k[,7])
  skew(DAT.fixed.k[,8])
  skew(DAT.fixed.k[,9])
  skew(DAT.fixed.k[,10])
  skew(DAT.fixed.k[,11])
  #

  ### Define a matrix for making

  OBS.high=length(gam.true)

  II	<- matrix(0,OBS.high,OBS.high)
  for(i in 1:OBS.high){
    II[i,i:OBS.high]	<- 1
  }
  Gam		<-	gam.true * II


  library(MASS)
  library(MCMCpack)
  library(cluster)
  require(graphics)

  #### SAMPLES TAKEN FROM THE SIMULATED INDIVIDUALS DATA
  DAT.1	<- NULL #data.frame(matrix(0,1,2))
  #for(j in 1:2){

  YYY			<- input.dat$sample.size[XXX]

  samp.numbs	<- c(0,YYY,YYY,YYY,YYY,YYY,YYY,YYY,0,0)

  for(i in 1:length(samp.numbs)){
    this.year	<- sample(DAT.fixed.k[,i],samp.numbs[i],replace=FALSE)
    TEMP	<- cbind(cohort,i,this.year)
    if(samp.numbs[i]>0){
      DAT.1	<- rbind(DAT.1,TEMP)
    }
  }

  DAT.1				<-	data.frame(DAT.1)
  colnames(DAT.1)[2:3]	<- c("date","X")


  ############################################# CHECK TO MAKE SURE THE EXPECTATION AND VARIANCE FUNCITONS ARE WORKING...
  REAL.mean 	<- mean(DAT.fixed.k)
  REAL.sd		<- sd(DAT.fixed.k)

  c(x0.true,apply(Gam,	2,	EXPECT,	X0=x0.true, ALPHA=A.true,BETA=B.true, PSI=psi.true))
  REAL.mean

  c(0,sqrt(apply(Gam,	2,	VAR,	X0=x0.true, DELTA = DELTA.true,	ALPHA=A.true,	BETA=B.true, PSI=psi.true) -
             apply(Gam,	2,	EXPECT,	X0=x0.true,	ALPHA=A.true,	BETA=B.true, PSI=psi.true)^2))
  REAL.sd



  #which(DAT.fixed.k[,9] > DAT.fixed.k[,10])


  #### MAKE A PLOT OF SOME GROWTH TRAJECTORIES.
  par(mfrow=c(1,1))
  y.lim=c(0,max(DAT.fixed.k[1:20,]))
  for(i in 1:20){
    plot(y=DAT.fixed.k[i,1:10],x=(1:10),ylim=y.lim,type="l")
    par(new=T)
  }


  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ### FITTING PARAMETERS TO THE DATA
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################

  ####################################################################################################
  ## ESTABLISH STARTING CONDITIONS
  ####################################################################################################
  THINNING			<-	10
  burn.in				<-	input.dat$burn.in[XXX]			### Define length of burn-in
  #NUMB				<-	input.dat$second[XXX] 			### Define length of second MCMC
  NUMB					<-	input.dat$final[XXX] +THINNING	### Define length of monitored MCMC
  RESTART 			<-	10000				### Define how often (in iterations) to clean out the monitored MCMC matrix to avoid holding large matrices in memory.
  #Define thining rate of the monitored MCMC
  lags		<-  50000
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################

  ### DATA MANIPULATION

  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ### MAKE A DESIGN MATRIX for comparing mus and sigmas to the dat
  ### FOR EQUALLY SPACED OBSERVATIONS!!!

  ### IDENTIFY WHICH YEARS NEED TO BE MODELED
  OBS.low  <- min(unique(DAT.1$date))
  OBS.high <- max(unique(DAT.1$date))
  OBS.seq	 <- 1:(OBS.high)
  OBS		<-	unique(DAT.1$date)

  #### MAKE A DUMMY MATRIX TO DEAL WITH GAMMA values

  #MAKE A DUMMY VARIABLE THAT CAN MAKE SURE WE DON'T HAVE TIMES WITH NO DATA ACCIDENTALLY CONTRIBUTING TO THE LIKELIHOOD.
  DUMMY						<-	rep(0,OBS.high)
  DUMMY[unique(DAT.1$date)]	<-	1

  X.design	<-	matrix(0,length(DAT.1$date),OBS.high)
  x.n			<-  vector(mode="numeric",length=OBS.high)

  for(i in 1:OBS.high){
    X.design[which(DAT.1$date==OBS.seq[i]),i]	<- 1
  }

  ### SAMPLE SIZES FOR THE AGE CLASSES
  x.n	<-	colSums(X.design)
  #x.n = x.n + 10e-20

  ## VECTOR OF OBSERVED LENGTHS
  x.dat		<-	DAT.1$X

  X.design.NA <-	X.design
  X.design.NA[X.design.NA == 0] = NA

  DAT.means = colMeans(X.design.NA*x.dat,na.rm=T)
  DAT.sd = sd(X.design.NA*x.dat,na.rm=T)
  #########################################################
  ### Define a matrix for making
  II	<- matrix(0,(OBS.high-1),(OBS.high-1))
  for(i in 1:(OBS.high-1)){
    II[i,i:(OBS.high-1)]	<- 1
  }

  #### SET UP MATRICIES TO HOLD ALL OF THE OUTPUT.
  # Parameters
  X0.out		<- matrix(0,burn.in,1)
  DELTA.out	<- matrix(0,burn.in,1)
  ALP.out 	<- matrix(0,burn.in,1)
  BETA.out	<- matrix(0,burn.in,1)
  GAMMA.out	<- matrix(0,burn.in,OBS.high-1)
  PSI.out		<- matrix(0,burn.in,1)
  TAU.out		<- matrix(0,burn.in,1)

  # Latent Variables
  MU.out		<- matrix(0,burn.in,OBS.high)
  SIG.out		<- matrix(0,burn.in,OBS.high)
  #### MAKE counters for the acceptance rate
  accept.x0		<- 0
  accept.alpha	<- 0
  accept.beta		<- 0
  accept.gamma	<- 0
  accept.delta	<- 0
  accept.psi		<- 0
  ####################################################################################################
  # Normal Posterior FOR X0
  x0.bar	<-	x0.true
  x0.sig	<-	0.001

  # FIXED starting standard deviation
  nu.sig	<-	0.0001
  ####################################################################################################
  ## ESTABLISH prior distributions
  # Gamma for delta
  a.delta	<-	0.00001
  b.delta	<-	0.00001
  # GAMMA Prior for ALPHA and BETA
  a.prior	<-	0.001
  b.prior	<-	0.001
  # Normal for the Qs (each year has the same prior)
  gam.bar	<- 100
  gam.sig	<- 100000
  # BETA for PSI
  a.psi	<- 1
  b.psi	<- 1
  # Inverse Gamma for NU
  #	tau.a	<-	0.001
  #	tau.b	<-	0.001
  ####################################################################################################
  #	Proposal Variances
  x0.sd			<- 0.5
  delta.sd		<- 0.10
  alpha.beta.CV	<- 0.15
  alpha.sd		<- 0.03
  beta.sd			<- 0.03
  cor.alp.beta <- 0.5
  cov.A.B		 <-  cor.alp.beta / (alpha.sd*beta.sd)
  #		C
  gamma.sd	<- 0.3
  psi.sd		<- 0.005
  #	tau.sd		<- 0.2
  #	pi.sd		<- 0.5
  ###################################################################################################
  #Starting Conditions
  x0.start		<- 15
  delta.start		<- 100
  alpha.start		<- 100
  beta.start		<- 1000
  gam.start		<- rep(50,OBS.high-1)
  psi.start		<- runif(1)
  #	tau.start		<- 1

  print(NAME)

  ###################################################################################################
  ###################################################################################################
  ## - START THE MCMC SAMPLER
  ###################################################################################################
  ###################################################################################################

  AA <- date()

  for(i in 1:burn.in){
    if(i == 1){
      x0.cur		<- 	x0.start
      delta.cur	<-	delta.start
      alp.cur		<- 	alpha.start
      beta.cur	<- 	beta.start
      gam.cur		<-	gam.start
      psi.cur		<-	psi.start
      #			tau.cur		<-	tau.start
    }
    if(i > 1){
      x0.cur		<- 	X0.out[i-1]
      delta.cur	<-	DELTA.out[i-1]
      alp.cur		<- 	ALP.out[i-1]
      beta.cur	<- 	BETA.out[i-1]
      gam.cur		<-	GAMMA.out[i-1,]
      psi.cur		<-	PSI.out[i-1,]
      #			tau.cur		<-	TAU.out[i-1,]
    }

    ###	Choose a value for X0
    x0.cur		<-	rnorm(1,x0.bar,x0.sig)
    MU.1.cur		<- x0.cur
    SIG.1.cur		<- nu.sig

    ####	Proposals for new parameters
    delta.star		<-	exp(rnorm(1,log(delta.cur),delta.sd))
    delta.star[delta.star < 0 ]	<- delta.cur
    #			w.star = w.true

    #		A.B.star	<- 	mvrnorm(1,c(alp.cur,beta.cur),matrix(c((alp.cur * alpha.beta.CV)^2,0,0,(alp.cur * alpha.beta.CV)^2),2,2))
    #		A.B.star	rnorm(1,c(alp.cur,beta.cur),matrix(c((alpha.sd)^2,0,0,(beta.sd)^2),2,2))
    #		A.B.star	<-	mvrnorm(1,c(log(alp.cur),log(beta.cur)),matrix(c((alpha.sd)^2,cov.A.B,cov.A.B,(beta.sd)^2),2,2))
    #				alp.star	<- A.B.star[1]
    #				beta.star	<- A.B.star[2]
    alp.star	<- exp(rnorm(1,log(alp.cur),alpha.sd))
    beta.star	<- exp(rnorm(1,log(beta.cur),beta.sd))
    #		print(paste("RAW PROPOSAL, ALP.star =",round(alp.star,2),"::: BETA.star =",round(beta.star,2)))
    alp.star[alp.star < 2 | beta.star < 2 ]		<- alp.cur
    beta.star[alp.star == alp.cur | beta.star < 1 ]	<- beta.cur
    gam.star		<-	rep(rnorm(1,gam.cur,gamma.sd),	length(gam.start))		#	mvrnorm(1,gam.cur,diag(rep(gamma.sd^2,9)))
    gam.star[gam.star< 0] <- gam.cur[gam.star< 0]
    psi.star		<-	rnorm(1,psi.cur,psi.sd)
    psi.star[psi.star < 0 | psi.star > 1]	<- psi.cur
    #psi.star = psi.start
    #		tau.star		<-	0.5

    ###################################################################################################
    # METROPOLIS STEP FOR ALPHA and BETA
    Z.cur		<- II*gam.cur

    MU.cur		<- apply(	Z.cur,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, 	PSI = psi.cur)
    MU.star		<- apply(	Z.cur,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.star,	BETA=beta.star, PSI = psi.cur)

    SIG.cur 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, 	PSI = psi.cur) - MU.cur^2)
    SIG.star 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.star,	BETA=beta.star, PSI = psi.cur) - MU.star^2)

    if(is.nan(SIG.star)[1] == FALSE){

      MU.cur 		<-	c(MU.1.cur,MU.cur)
      MU.star		<-	c(MU.1.cur,MU.star)
      SIG.cur		<-	c(SIG.1.cur, SIG.cur)
      SIG.star	<-	c(SIG.1.cur, SIG.star)

      #	SIG.cur[is.nan(SIG.cur)==TRUE]		<- 10e-10
      #	SIG.star[is.nan(SIG.star)==TRUE]	<- 10e-10

      PHI.gam.star	<-	MU.star / SIG.star^2
      THETA.gam.star	<-	PHI.gam.star * MU.star

      PHI.gam.cur		<-	MU.cur / SIG.cur^2
      THETA.gam.cur	<-	PHI.gam.cur * MU.cur

      #	if(is.nan(MU.star)[1] == FALSE){

      ## Compare CURRENT AND PROPOSED VALUES ALPHA and BETA
      TEMP <- log((X.design*x.dat))
      TEMP[is.infinite(TEMP)==TRUE]	<- 0

      CUR 	<-  sum(x.n * (THETA.gam.cur * log(PHI.gam.cur))) - sum(x.n * lgamma(THETA.gam.cur)) +
        sum(TEMP %*% (THETA.gam.cur - 1)) - sum((X.design*x.dat) %*% PHI.gam.cur) +
        (a.prior-1)*log(alp.cur) - b.prior*alp.cur + (a.prior-1)*log(beta.cur) - b.prior*beta.cur +
        log(alp.cur) + log(beta.cur)

      STAR 	<- sum(x.n * (THETA.gam.star * log(PHI.gam.star))) - sum(x.n * lgamma(THETA.gam.star)) +
        sum(TEMP %*% (THETA.gam.star - 1)) - sum((X.design*x.dat) %*% PHI.gam.star) +
        (a.prior-1)*log(alp.star) - b.prior*alp.star + (a.prior-1)*log(beta.star) - b.prior*beta.star +
        log(alp.star) + log(beta.star)

      accept.prop	<- 	min(exp(STAR - CUR),1)
      keep	<-	ceiling(accept.prop - runif(1))
      #		}
    }
    #print("SPOT 1")
    if(keep == 1) {
      alp.cur 	= alp.star; accept.alpha = accept.alpha + 1
      beta.cur 	= beta.star; accept.beta = accept.beta + 1
    }
    #		}
    ###################################################################################################
    # METROPOLIS STEP FOR DELTA
    Z.cur			<- II*gam.cur

    MU.cur		<- apply(	Z.cur,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur)

    SIG.cur 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur,  DELTA=delta.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.cur^2)
    SIG.star 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, 	DELTA=delta.star,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.cur^2)

    if(is.nan(SIG.star)[1] == FALSE){

      MU.cur 		<-	c(MU.1.cur,MU.cur)

      SIG.cur		<-	c(SIG.1.cur, SIG.cur)
      SIG.star	<-	c(SIG.1.cur, SIG.star)

      #	SIG.cur[is.nan(SIG.cur)==TRUE]		<- 10e-10
      #	SIG.star[is.nan(SIG.star)==TRUE]	<- 10e-10

      PHI.gam.star	<-	MU.cur / SIG.star^2
      THETA.gam.star	<-	PHI.gam.star * MU.cur

      PHI.gam.cur		<-	MU.cur / SIG.cur^2
      THETA.gam.cur	<-	PHI.gam.cur * MU.cur

      ## Compare CURRENT AND PROPOSED VALUES for x0 and NU
      TEMP <- log((X.design*x.dat))
      TEMP[is.infinite(TEMP)==TRUE]	<- 0

      CUR 	<-  sum(x.n * (THETA.gam.cur * log(PHI.gam.cur))) - sum(x.n * lgamma(THETA.gam.cur)) +
        sum(TEMP %*% (THETA.gam.cur - 1)) - sum((X.design*x.dat) %*% PHI.gam.cur) +
        (a.delta-1)*log(delta.cur) - b.delta * delta.cur +
        log(delta.cur)

      STAR 	<- sum(x.n * (THETA.gam.star * log(PHI.gam.star))) - sum(x.n * lgamma(THETA.gam.star)) +
        sum(TEMP %*% (THETA.gam.star - 1)) - sum((X.design*x.dat) %*% PHI.gam.star) +
        (a.delta-1)*log(delta.star) - b.delta * delta.star +
        log(delta.star)

      accept.prop	<- 	min(exp(STAR - CUR),1)
      keep	<-	ceiling(accept.prop - runif(1))
    }
    #	print("SPOT 2")
    if(keep == 1) {
      #print("Yes")
      delta.cur		= delta.star; accept.delta = accept.delta +1
    }

    ###################################################################################################
    ###################################################################################################
    # METROPOLIS STEP FOR GAMMAS
    Z.cur			<- II*gam.cur
    Z.star			<- II*gam.star

    MU.cur		<- apply(	Z.cur,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur)
    MU.star		<- apply(	Z.star,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur)

    SIG.cur 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.cur^2)
    SIG.star 	<- sqrt(apply( Z.star,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.star^2)

    if(is.nan(SIG.star)[1] == FALSE){

      MU.cur 	<-	c(MU.1.cur,MU.cur)
      MU.star <-	c(MU.1.cur,MU.star)

      SIG.cur		<-	c(SIG.1.cur, SIG.cur)
      SIG.star	<-	c(SIG.1.cur, SIG.star)

      #	SIG.cur[is.nan(SIG.cur)==TRUE]		<- 10e-10
      #	SIG.star[is.nan(SIG.star)==TRUE]	<- 10e-10

      PHI.gam.star	<-	MU.star / SIG.star^2
      THETA.gam.star	<-	PHI.gam.star * MU.star

      PHI.gam.cur		<-	MU.cur / SIG.cur^2
      THETA.gam.cur	<-	PHI.gam.cur * MU.cur

      ## Compare CURRENT AND PROPOSED VALUES for x0 and NU
      TEMP <- log((X.design*x.dat))
      TEMP[is.infinite(TEMP)==TRUE]	<- 0

      CUR 	<-  sum(x.n * (THETA.gam.cur * log(PHI.gam.cur))) - sum(x.n * lgamma(THETA.gam.cur)) +
        sum(TEMP %*% (THETA.gam.cur - 1)) - sum((X.design*x.dat) %*% PHI.gam.cur)
      #	(1/(2*gam.sig^2))*(gam.cur-gam.bar)^2

      STAR 	<- sum(x.n * (THETA.gam.star * log(PHI.gam.star))) - sum(x.n * lgamma(THETA.gam.star)) +
        sum(TEMP %*% (THETA.gam.star - 1)) - sum((X.design*x.dat) %*% PHI.gam.star)
      #	(1/(2*gam.sig^2))*(gam.star-gam.bar)^2

      accept.prop	<- 	min(exp(STAR - CUR),1)
      keep	<-	ceiling(accept.prop - runif(1))
    }
    #	print("SPOT 2")
    if(keep == 1) {
      #print("Yes")
      #x0.cur		= x0.star; accept.x0 = accept.x0 + 1
      gam.cur		= gam.star; accept.gamma = accept.gamma +1
    }

    ###################################################################################################
    ###################################################################################################
    # METROPOLIS STEP FOR PSI
    Z.cur			<- II*gam.cur

    MU.cur		<- apply(	Z.cur,	2,	EXPECT,		X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur)
    MU.star		<- apply(	Z.cur,	2,	EXPECT,		X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.star)

    SIG.cur 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.cur^2)
    SIG.star 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.star) - MU.star^2)

    if(is.nan(SIG.star)[1] == FALSE){

      MU.cur 		<-	c(MU.1.cur,MU.cur)
      MU.star 	<-	c(MU.1.cur,MU.star)

      SIG.cur		<-	c(SIG.1.cur, SIG.cur)
      SIG.star	<-	c(SIG.1.cur, SIG.star)

      #	SIG.cur[is.nan(SIG.cur)==TRUE]		<- 10e-10
      #	SIG.star[is.nan(SIG.star)==TRUE]	<- 10e-10

      PHI.gam.star	<-	MU.star / SIG.star^2
      THETA.gam.star	<-	PHI.gam.star * MU.star

      PHI.gam.cur		<-	MU.cur / SIG.cur^2
      THETA.gam.cur	<-	PHI.gam.cur * MU.cur

      ## Compare CURRENT AND PROPOSED VALUES for PSI
      TEMP <- log((X.design*x.dat))
      TEMP[is.infinite(TEMP)==TRUE]	<- 0

      CUR 	<-  sum(x.n * (THETA.gam.cur * log(PHI.gam.cur))) - sum(x.n * lgamma(THETA.gam.cur)) +
        sum(TEMP %*% (THETA.gam.cur - 1)) - sum((X.design*x.dat) %*% PHI.gam.cur) +
        (a.psi - 1) * log(psi.cur) +  (b.psi-1) * log(1-psi.cur)
      STAR 	<- sum(x.n * (THETA.gam.star * log(PHI.gam.star))) - sum(x.n * lgamma(THETA.gam.star)) +
        sum(TEMP %*% (THETA.gam.star - 1)) - sum((X.design*x.dat) %*% PHI.gam.star) +
        (a.psi - 1) * log(psi.star) +  (b.psi-1) * log(1-psi.star)

      accept.prop	<- 	min(exp(STAR - CUR),1)
      keep	<-	ceiling(accept.prop - runif(1))
    }
    #print("SPOT 4")
    if(keep == 1) {
      #	print("Yes")
      psi.cur		= psi.star; accept.psi = accept.psi + 1
    }

    ####################################################################################################
    ### SAVE THE current values to the appropriate matrices
    X0.out[i] 		<- x0.cur
    DELTA.out[i]	<- delta.cur
    ALP.out[i] 		<- alp.cur
    BETA.out[i] 	<- beta.cur
    GAMMA.out[i,] 	<- gam.cur
    PSI.out[i]		<- psi.cur

    MU.out[i,]		<- MU.cur
    SIG.out[i,]		<- SIG.cur

    counter	<- c(1,seq(1000,burn.in,by=1000))
    if(max(i==counter)) {
      print(paste(counter[which(counter==i)],date()))

      print(accept.alpha/i)
      print(accept.beta/i)
      print(accept.gamma/i)
      print(accept.delta/i)
      print(accept.psi/i)

    }

  }	 #### END MCMC LOOP.
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  X0.burn.in		<-	X0.out
  DELTA.burn.in	<-	DELTA.out
  ALP.burn.in		<-	ALP.out
  BETA.burn.in	<-	BETA.out
  GAMMA.burn.in	<-	GAMMA.out
  PSI.burn.in		<-	PSI.out

  MU.burn.in		<-	MU.out
  SIG.burn.in		<-	SIG.out

  # Calculate an empirical Var-Cov matrix from the accepted MCMC draws

  TEMP.1 <- cbind(X0.out[(burn.in/4):burn.in],
                  log(ALP.out[(burn.in/4):burn.in]),
                  log(BETA.out[(burn.in/4):burn.in]),
                  log(DELTA.out[(burn.in/4):burn.in]))

  TEMP.2 <- cbind(PSI.out[(burn.in/4):burn.in],
                  GAMMA.out[(burn.in/4):burn.in,1])

  NEW.PROP.1	<- cov(TEMP.1)
  NEW.PROP.2	<- cov(TEMP.2)

  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ## MAKE SOME PLOTS
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################

  BB	<- date()

  print(paste("Start :", AA,":::: Stop", BB))

  par(mfrow=c(3,3))
  x.lim=c(0,burn.in)

  plot(X0.out,xlim=x.lim,pch=".")#ylim=c(0,20))
  abline(h=x0.true,col=2)
  plot(DELTA.out[,1],xlim=x.lim,pch=".")
  abline(h = DELTA.true,col=2)
  plot(ALP.out,xlim=x.lim,pch=".")
  abline(h=A.true,col=2)
  plot(BETA.out,xlim=x.lim,pch=".")
  abline(h=B.true,col=2)
  plot(PSI.out,xlim=x.lim,pch=".")
  abline(h=psi.true,col=2)
  plot(GAMMA.out[,1],xlim=x.lim,pch=".",ylim= c(min(GAMMA.out,gam.true*0.9),max(GAMMA.out,gam.true*1.1)))
  abline(h = gam.true[1],col=2)

  plot(ALP.out/BETA.out,ylim=c(0,0.5),pch=".")
  title(paste("mean K =", round(mean(ALP.out[0.9*i:i]/BETA.out[0.9*i:i]),3),":::: True K =",round(A.true/B.true,3)),line=0.5)
  abline(h= A.true/B.true,col=2)

  plot(1/sqrt(ALP.out),main=paste("CV k",round(mean(1/sqrt(ALP.out)),3),":::: True = ",1/sqrt(A.true)),,pch=".")
  #	main= paste("mean CV =", round(sqrt(ALP.out),3),":::: True CV =",round(1/sqrt(k.A),3)))
  abline(h= 1/sqrt(A.true),col=2)

  plot(sqrt(ALP.out/BETA.out^2),main=paste("SD k",round(sqrt(ALP.out[0.9*i:i]/BETA.out[0.9*i:i]^2),3),":::: True = ",round(sqrt(A.true/B.true^2),3)),pch=".")
  #		main= paste("mean SD k =", round(sqrt(ALP.out[0.9*i:i]/BETA.out[0.9*i:i]^2)),3),":::: True SD k =",round(k.A/k.B^2,3))
  abline(h= sqrt(A.true/B.true^2),col=2)


  ### MEANS....
  #par(mfrow=c(1,1))
  #	y.lim <- c(0,50)
  #	plot(MU.out[,2],ylim=y.lim,col=3,pch=".")
  #		abline(h=REAL.mean[2],lwd=2,col="darkgreen",pch=".")
  #		par(new=T)
  #	plot(MU.out[,3],ylim=y.lim,col=2,pch=".")
  #		abline(h=REAL.mean[3],lwd=2,col=2,pch=".")
  #		par(new=T)
  #	plot(MU.out[,4],ylim=y.lim,col=4,pch=".")
  #		abline(h=REAL.mean[4],lwd=2,col=4)
  #		par(new=T)
  #	plot(MU.out[,5],ylim=y.lim,col=5,pch=".")
  #		abline(h=REAL.mean[5],lwd=2,col=5)
  #		par(new=T)
  #	plot(MU.out[,6],ylim=y.lim,col=6,pch=".")
  #		abline(h=REAL.mean[6],lwd=2,col=6)
  #		par(new=T)
  #	plot(MU.out[,7],ylim=y.lim,col=7,pch=".")
  #		abline(h=REAL.mean[7],lwd=2,col=7)
  #		par(new=T)
  #	plot(MU.out[,8],ylim=y.lim,col=8,pch=".")
  #		abline(h=REAL.mean[8],lwd=2,col=8)
  ##		par(new=T)
  #
  ##### VARIANCES....
  #	par(mfrow=c(1,1))
  #	y.lim <- c(0,REAL.sd[8]*1.1)
  #	plot(SIG.out[,2],ylim=y.lim,col=2,pch=".")
  #		abline(h=REAL.sd[2],lwd=2,col=2)
  #		par(new=T)
  #	plot(SIG.out[,3],ylim=y.lim,col=4,pch=".")
  #		abline(h=REAL.sd[3],lwd=2,col=4)
  #		par(new=T)
  #	plot(SIG.out[,4],ylim=y.lim,col=3,pch=".")
  #		abline(h=REAL.sd[4],lwd=2,col=3)
  #		par(new=T)
  #	plot(SIG.out[,5],ylim=y.lim,col=5,pch=".")
  #		abline(h=REAL.sd[5],lwd=2,col=5)
  #		par(new=T)
  #	plot(SIG.out[,6],ylim=y.lim,col=6,pch=".")
  #		abline(h=REAL.sd[6],lwd=2,col=6)
  #		par(new=T)
  #	plot(SIG.out[,7],ylim=y.lim,col=7,pch=".")
  #		abline(h=REAL.sd[7],lwd=2,col=7)
  #		par(new=T)
  #	plot(SIG.out[,8],ylim=y.lim,col=8,pch=".")
  #		abline(h=REAL.sd[8],lwd=2,col=8)

  ###############################################################################################
  ### HISTOGRAMS OF LENGTHS of Sampled INDIVIDUALS and estimated lengths
  ###############################################################################################
  #
  #all.MU	<- colMeans(MU.out[(burn.in/2):burn.in,])
  #all.SIG	<- colMeans(SIG.out[(burn.in/2):burn.in,])
  #
  #par(mfrow=c(3,4))
  #
  #for(i in 1:OBS.high){
  #
  #	Z <- DAT.1$X[DAT.1$date==i]
  #	if(samp.numbs[i] > 0){
  #		x.lim	<- c(min(Z)*0.99, max(Z)*1.01)
  #		hist(Z,xlim=x.lim,main=paste("Length of age ",i),breaks=20)
  #		par(new=T)
  #		plot(x = seq(x.lim[1], x.lim[2],length=100),
  #			 y = dgamma( seq(x.lim[1], x.lim[2],length=100),
  #								 all.MU[i]^2 / all.SIG[i]^2, all.MU[i] / all.SIG[i]^2),
  #			xlim=x.lim,type="l",col=2,lwd=2,
  #			axes=F,xlab="",ylab="")
  #		}
  #	if(samp.numbs[i]==0){
  #		plot(1:5,1:5,axes=F,col="white",xlab="",ylab="")
  #		}
  #	}
  #
  ###############################################################################################
  ### HISTOGRAMS OF LENGTHS of ALL SIMULATED INDIVIDUALS
  ###############################################################################################
  #
  #par(mfrow=c(3,4))
  #	hist(DAT.fixed.k[,1])
  #	hist(DAT.fixed.k[,2])
  #	hist(DAT.fixed.k[,3])
  #	hist(DAT.fixed.k[,4])
  #	hist(DAT.fixed.k[,5])
  #	hist(DAT.fixed.k[,6])
  #	hist(DAT.fixed.k[,7])
  #	hist(DAT.fixed.k[,8])
  #	hist(DAT.fixed.k[,9])
  #	hist(DAT.fixed.k[,10])
  #	hist(DAT.fixed.k[,11])
  ##	hist(DAT.fixed.k[,12])
  ##	hist(DAT.fixed.k[,13])
  ##	hist(DAT.fixed.k[,14])
  ##	hist(DAT.fixed.k[,15])
  ##	hist(DAT.fixed.k[,16])
  ##	hist(DAT.fixed.k[,17])
  ##	hist(DAT.fixed.k[,18])
  ##	hist(DAT.fixed.k[,19])
  ##	hist(DAT.fixed.k[,20])
  #


  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################

  #  FINAL MCMC

  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  #### SET UP MATRICIES TO HOLD ALL OF THE OUTPUT.
  # Parameters
  X0.out			<- matrix(0,NUMB,1)
  DELTA.out		<- matrix(0,NUMB,1)
  ALP.out 		<- matrix(0,NUMB,1)
  BETA.out		<- matrix(0,NUMB,1)
  GAMMA.out		<- matrix(0,NUMB,OBS.high-1)
  PSI.out			<- matrix(0,NUMB,1)
  #TAU.out	<- matrix(0,burn.in,1)

  # Latent Variables
  MU.out		<- matrix(0,NUMB,OBS.high)
  SIG.out		<- matrix(0,NUMB,OBS.high)

  #Starting Conditions
  x0.start		<- X0.burn.in[burn.in]
  alpha.start		<- ALP.burn.in[burn.in]
  beta.start		<- BETA.burn.in[burn.in]
  gamma.start		<- GAMMA.burn.in[burn.in,]
  delta.start		<- DELTA.burn.in[burn.in]
  psi.start		<- PSI.burn.in[burn.in]

  #### MAKE counters for the acceptance rate
  accept.x0		<- 0
  accept.alpha	<- 0
  accept.beta		<- 0
  accept.gamma	<- 0
  accept.delta	<- 0
  accept.psi		<- 0
  ############################# START MCMC
  for(i in 1:NUMB){
    if(i == 1){
      x0.cur		<- 	x0.start
      delta.cur	<-	delta.start
      alp.cur		<- 	alpha.start
      beta.cur	<- 	beta.start
      gam.cur		<-	gamma.start
      delta.cur	<-	delta.start
      psi.cur		<-	psi.start
      #			tau.cur		<-	tau.start
    }
    if(i > 1){
      x0.cur		<- 	X0.out[i-1]
      delta.cur	<-	DELTA.out[i-1]
      alp.cur		<- 	ALP.out[i-1]
      beta.cur	<- 	BETA.out[i-1]
      gam.cur		<-	GAMMA.out[i-1,]
      psi.cur		<-	PSI.out[i-1,]
      #			tau.cur		<-	TAU.out[i-1,]
    }

    ###	Choose a value for X0
    x0.cur		<-	rnorm(1,x0.bar,x0.sig)
    MU.1.cur		<- x0.cur
    SIG.1.cur		<- nu.sig

    ####	Proposals for new parameters

    ALL.prop.1	<-	mvrnorm(1, c(x0.cur,log(alp.cur),log(beta.cur),log(delta.cur)),NEW.PROP.1)
    ALL.prop.2	<-	mvrnorm(1, c(psi.cur,gam.cur[1]),NEW.PROP.2)

    alp.star	<-	exp(ALL.prop.1[2])
    beta.star	<-	exp(ALL.prop.1[3])
    delta.star	<-	exp(ALL.prop.1[4])
    psi.star	<-	ALL.prop.2[1]
    #psi.star<-  psi.start
    gam.star	<-	ALL.prop.2[2]
    gam.star<- rep(ALL.prop.2[2],length(gam.cur))

    ### ADJUSTMENTS to avoid ridiculous proposed parameters
    delta.star[delta.star < 0 ]	<- delta.cur
    alp.star[alp.star < 2 | beta.star < 2 ]		<- alp.cur
    beta.star[alp.star == alp.cur | beta.star < 1 ]	<- beta.cur
    gam.star[gam.star< 0] <- gam.cur[gam.star< 0]
    psi.star[psi.star < 0 | psi.star > 1]	<- psi.cur
    ###################################################################################################
    # METROPOLIS STEP FOR ALPHA and BETA
    Z.cur		<- II*gam.cur

    MU.cur		<- apply(	Z.cur,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, 	PSI = psi.cur)
    MU.star		<- apply(	Z.cur,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.star,	BETA=beta.star, PSI = psi.cur)

    SIG.cur 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, 	PSI = psi.cur) - MU.cur^2)
    SIG.star 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.star,	BETA=beta.star, PSI = psi.cur) - MU.star^2)

    if(is.nan(SIG.star)[1] == FALSE){

      MU.cur 		<-	c(MU.1.cur,MU.cur)
      MU.star		<-	c(MU.1.cur,MU.star)
      SIG.cur		<-	c(SIG.1.cur, SIG.cur)
      SIG.star	<-	c(SIG.1.cur, SIG.star)

      #	SIG.cur[is.nan(SIG.cur)==TRUE]		<- 10e-10
      #	SIG.star[is.nan(SIG.star)==TRUE]	<- 10e-10

      PHI.gam.star	<-	MU.star / SIG.star^2
      THETA.gam.star	<-	PHI.gam.star * MU.star

      PHI.gam.cur		<-	MU.cur / SIG.cur^2
      THETA.gam.cur	<-	PHI.gam.cur * MU.cur

      #	if(is.nan(MU.star)[1] == FALSE){

      ## Compare CURRENT AND PROPOSED VALUES ALPHA and BETA
      TEMP <- log((X.design*x.dat))
      TEMP[is.infinite(TEMP)==TRUE]	<- 0

      CUR 	<-  sum(x.n * (THETA.gam.cur * log(PHI.gam.cur))) - sum(x.n * lgamma(THETA.gam.cur)) +
        sum(TEMP %*% (THETA.gam.cur - 1)) - sum((X.design*x.dat) %*% PHI.gam.cur) +
        (a.prior-1)*log(alp.cur) - b.prior*alp.cur + (a.prior-1)*log(beta.cur) - b.prior*beta.cur +
        log(alp.cur) + log(beta.cur)

      STAR 	<- sum(x.n * (THETA.gam.star * log(PHI.gam.star))) - sum(x.n * lgamma(THETA.gam.star)) +
        sum(TEMP %*% (THETA.gam.star - 1)) - sum((X.design*x.dat) %*% PHI.gam.star) +
        (a.prior-1)*log(alp.star) - b.prior*alp.star + (a.prior-1)*log(beta.star) - b.prior*beta.star +
        log(alp.star) + log(beta.star)

      accept.prop	<- 	min(exp(STAR - CUR),1)
      keep	<-	ceiling(accept.prop - runif(1))
      #		}
    }
    #print("SPOT 1")
    if(keep == 1) {
      alp.cur 	= alp.star; accept.alpha = accept.alpha + 1
      beta.cur 	= beta.star; accept.beta = accept.beta + 1
    }
    #		}
    ###################################################################################################
    # METROPOLIS STEP FOR DELTA
    Z.cur			<- II*gam.cur

    MU.cur		<- apply(	Z.cur,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur)

    SIG.cur 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur,  DELTA=delta.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.cur^2)
    SIG.star 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, 	DELTA=delta.star,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.cur^2)

    if(is.nan(SIG.star)[1] == FALSE){

      MU.cur 		<-	c(MU.1.cur,MU.cur)

      SIG.cur		<-	c(SIG.1.cur, SIG.cur)
      SIG.star	<-	c(SIG.1.cur, SIG.star)

      #	SIG.cur[is.nan(SIG.cur)==TRUE]		<- 10e-10
      #	SIG.star[is.nan(SIG.star)==TRUE]	<- 10e-10

      PHI.gam.star	<-	MU.cur / SIG.star^2
      THETA.gam.star	<-	PHI.gam.star * MU.cur

      PHI.gam.cur		<-	MU.cur / SIG.cur^2
      THETA.gam.cur	<-	PHI.gam.cur * MU.cur

      ## Compare CURRENT AND PROPOSED VALUES for x0 and NU
      TEMP <- log((X.design*x.dat))
      TEMP[is.infinite(TEMP)==TRUE]	<- 0

      CUR 	<-  sum(x.n * (THETA.gam.cur * log(PHI.gam.cur))) - sum(x.n * lgamma(THETA.gam.cur)) +
        sum(TEMP %*% (THETA.gam.cur - 1)) - sum((X.design*x.dat) %*% PHI.gam.cur) +
        (a.delta-1)*log(delta.cur) - b.delta * delta.cur + log(delta.cur)

      STAR 	<- sum(x.n * (THETA.gam.star * log(PHI.gam.star))) - sum(x.n * lgamma(THETA.gam.star)) +
        sum(TEMP %*% (THETA.gam.star - 1)) - sum((X.design*x.dat) %*% PHI.gam.star) +
        (a.delta-1)*log(delta.star) - b.delta * delta.star + log(delta.star)

      accept.prop	<- 	min(exp(STAR - CUR),1)
      keep	<-	ceiling(accept.prop - runif(1))
    }
    #	print("SPOT 2")
    if(keep == 1) {
      #print("Yes")
      delta.cur		= delta.star; accept.delta = accept.delta +1
    }

    ###################################################################################################
    ###################################################################################################
    # METROPOLIS STEP FOR GAMMAS
    Z.cur			<- II*gam.cur
    Z.star			<- II*gam.star

    MU.cur		<- apply(	Z.cur,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur)
    MU.star		<- apply(	Z.star,	2,	EXPECT,	X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur)

    SIG.cur 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.cur^2)
    SIG.star 	<- sqrt(apply( Z.star,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.star^2)

    if(is.nan(SIG.star)[1] == FALSE){

      MU.cur 	<-	c(MU.1.cur,MU.cur)
      MU.star <-	c(MU.1.cur,MU.star)

      SIG.cur		<-	c(SIG.1.cur, SIG.cur)
      SIG.star	<-	c(SIG.1.cur, SIG.star)

      #	SIG.cur[is.nan(SIG.cur)==TRUE]		<- 10e-10
      #	SIG.star[is.nan(SIG.star)==TRUE]	<- 10e-10

      PHI.gam.star	<-	MU.star / SIG.star^2
      THETA.gam.star	<-	PHI.gam.star * MU.star

      PHI.gam.cur		<-	MU.cur / SIG.cur^2
      THETA.gam.cur	<-	PHI.gam.cur * MU.cur

      ## Compare CURRENT AND PROPOSED VALUES for GAMMA
      TEMP <- log((X.design*x.dat))
      TEMP[is.infinite(TEMP)==TRUE]	<- 0

      CUR 	<-  sum(x.n * (THETA.gam.cur * log(PHI.gam.cur))) - sum(x.n * lgamma(THETA.gam.cur)) +
        sum(TEMP %*% (THETA.gam.cur - 1)) - sum((X.design*x.dat) %*% PHI.gam.cur) -
        (1/(2*gam.sig^2))*(gam.cur-gam.bar)^2

      STAR 	<- sum(x.n * (THETA.gam.star * log(PHI.gam.star))) - sum(x.n * lgamma(THETA.gam.star)) +
        sum(TEMP %*% (THETA.gam.star - 1)) - sum((X.design*x.dat) %*% PHI.gam.star) -
        (1/(2*gam.sig^2))*(gam.star-gam.bar)^2

      accept.prop	<- 	min(exp(STAR - CUR),1)
      keep	<-	ceiling(accept.prop - runif(1))
    }
    #	print("SPOT 2")
    if(keep == 1) {
      #print("Yes")
      #x0.cur		= x0.star; accept.x0 = accept.x0 + 1
      gam.cur		= gam.star; accept.gamma = accept.gamma +1
    }

    ###################################################################################################
    ###################################################################################################
    # METROPOLIS STEP FOR PSI
    Z.cur			<- II*gam.cur

    MU.cur		<- apply(	Z.cur,	2,	EXPECT,		X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur)
    MU.star		<- apply(	Z.cur,	2,	EXPECT,		X0=x0.cur,	ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.star)

    SIG.cur 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.cur) - MU.cur^2)
    SIG.star 	<- sqrt(apply( Z.cur,	2,	VAR,	X0=x0.cur, DELTA = delta.cur, ALPHA=alp.cur,	BETA=beta.cur, PSI = psi.star) - MU.star^2)

    if(is.nan(SIG.star)[1] == FALSE){

      MU.cur 		<-	c(x0.cur,MU.cur)
      MU.star 	<-	c(x0.cur,MU.star)

      SIG.cur		<-	c(SIG.1.cur, SIG.cur)
      SIG.star	<-	c(SIG.1.cur, SIG.star)

      #	SIG.cur[is.nan(SIG.cur)==TRUE]		<- 10e-10
      #	SIG.star[is.nan(SIG.star)==TRUE]	<- 10e-10

      PHI.gam.star	<-	MU.star / SIG.star^2
      THETA.gam.star	<-	PHI.gam.star * MU.star

      PHI.gam.cur		<-	MU.cur / SIG.cur^2
      THETA.gam.cur	<-	PHI.gam.cur * MU.cur

      ## Compare CURRENT AND PROPOSED VALUES for PSI
      TEMP <- log((X.design*x.dat))
      TEMP[is.infinite(TEMP)==TRUE]	<- 0

      CUR 	<-  sum(x.n * (THETA.gam.cur * log(PHI.gam.cur))) - sum(x.n * lgamma(THETA.gam.cur)) +
        sum(TEMP %*% (THETA.gam.cur - 1)) - sum((X.design*x.dat) %*% PHI.gam.cur) +
        (a.psi - 1) * log(psi.cur) +  (b.psi-1) * log(1-psi.cur)
      STAR 	<- sum(x.n * (THETA.gam.star * log(PHI.gam.star))) - sum(x.n * lgamma(THETA.gam.star)) +
        sum(TEMP %*% (THETA.gam.star - 1)) - sum((X.design*x.dat) %*% PHI.gam.star) +
        (a.psi - 1) * log(psi.star) +  (b.psi-1) * log(1-psi.star)

      accept.prop	<- 	min(exp(STAR - CUR),1)
      keep	<-	ceiling(accept.prop - runif(1))
    }
    #print("SPOT 4")
    if(keep == 1) {
      #	print("Yes")
      psi.cur		= psi.star; accept.psi = accept.psi + 1
    }

    ####################################################################################################
    ### SAVE THE current values to the appropriate matrices
    X0.out[i] 		<- x0.cur
    DELTA.out[i]	<- delta.cur
    ALP.out[i] 		<- alp.cur
    BETA.out[i] 	<- beta.cur
    GAMMA.out[i,] 	<- gam.cur
    PSI.out[i]		<- psi.cur

    MU.out[i,]		<- MU.cur
    SIG.out[i,]		<- SIG.cur

    counter	<- c(1,seq(1000,NUMB,by=1000))
    if(max(i==counter)) {
      print(paste(counter[which(counter==i)],date()))
    }

  }	 #### END MCMC LOOP.
  ##################################################
  ####################################################################################################

  ### SAVE THE FINAL MCMC STUFF
  X0.FINAL		<-	X0.out
  DELTA.FINAL		<-	DELTA.out
  ALP.FINAL		<-	ALP.out
  BETA.FINAL		<-	BETA.out
  GAMMA.FINAL		<-	GAMMA.out
  PSI.FINAL		<-	PSI.out

  MU.FINAL		<-	MU.out
  SIG.FINAL		<-	SIG.out

  ###### THIN THE CHAINS
  THIN <- seq(THINNING,NUMB,by=THINNING)

  X0.thin			<- X0.FINAL[THIN]
  ALP.thin		<- ALP.FINAL[THIN]
  BETA.thin		<- BETA.FINAL[THIN]
  GAMMA.thin		<- GAMMA.FINAL[THIN,]
  DELTA.thin		<- DELTA.FINAL[THIN]
  PSI.thin		<- PSI.FINAL[THIN]

  MU.thin			<- MU.FINAL[THIN,]
  SIG.thin		<- SIG.FINAL[THIN,]

  k.bar.thin		<- ALP.thin / BETA.thin
  k.CV.thin		<-	1/ sqrt(ALP.thin)
  ##### MAKE SOME PLOTS of thinned final MCMC output.

  ####################################################################################################
  ####################################################################################################

  ###  MAKE AN OUTPUT FILE OF PLOTS

  setwd("/Users/ole/Documents/Von B Growth/OUTPUT/Observational samples/Plots")

  pdf(file=paste(NAME,".MCMC.output.pdf"),onefile=TRUE,width=7,height=7)

  ### PLOT TRACES of BURN-IN MCMC.

  par(mfrow=c(3,3))
  x.lim=c(0,burn.in)

  plot(X0.burn.in,xlim=x.lim,pch=".",main="BURN IN")#ylim=c(0,20))
  abline(h=x0.true,col=2)
  plot(DELTA.burn.in[,1],xlim=x.lim,pch=".")
  abline(h = DELTA.true,col=2)
  plot(ALP.burn.in,xlim=x.lim,pch=".")
  abline(h=A.true,col=2)
  plot(BETA.burn.in,xlim=x.lim,pch=".")
  abline(h=B.true,col=2)
  plot(PSI.burn.in,xlim=x.lim,pch=".")
  abline(h=psi.true,col=2)
  plot(GAMMA.burn.in[,1],xlim=x.lim,pch=".",ylim= c(min(GAMMA.out,gam.true*0.9),max(GAMMA.out,gam.true*1.1)))
  abline(h = gam.true[1],col=2)

  plot(ALP.burn.in/BETA.burn.in,ylim=c(0,0.5),pch=".")
  title(paste("mean K =", round(mean(ALP.out[0.9*i:i]/BETA.out[0.9*i:i]),3),":::: True K =",round(A.true/B.true,3)),line=0.5)
  abline(h= A.true/B.true,col=2)

  plot(1/sqrt(ALP.out),main=paste("CV k",round(mean(1/sqrt(ALP.out)),3),":::: True = ",1/sqrt(A.true)),,pch=".")
  #	main= paste("mean CV =", round(sqrt(ALP.out),3),":::: True CV =",round(1/sqrt(k.A),3)))
  abline(h= 1/sqrt(A.true),col=2)

  plot(sqrt(ALP.out/BETA.out^2),main=paste("SD k",round(sqrt(ALP.out[0.9*i:i]/BETA.out[0.9*i:i]^2),3),":::: True = ",round(sqrt(A.true/B.true^2),3)),pch=".")
  #		main= paste("mean SD k =", round(sqrt(ALP.out[0.9*i:i]/BETA.out[0.9*i:i]^2)),3),":::: True SD k =",round(k.A/k.B^2,3))
  abline(h= sqrt(A.true/B.true^2),col=2)

  ### PLOT TRACES of THIRD MCMC.
  par(mfrow=c(3,3),mar=c(4,5,2,0.5))
  x.lim=c(0,length(THIN))

  plot(X0.thin,xlim=x.lim,pch=".")#ylim=c(0,20))
  abline(h=x0.true,col=2)
  plot(ALP.thin,xlim=x.lim,pch=".")
  abline(h=A.true,col=2)
  plot(BETA.thin,xlim=x.lim,pch=".")
  abline(h=B.true,col=2)
  plot(PSI.thin,xlim=x.lim,pch=".")
  abline(h=psi.true,col=2)
  plot(GAMMA.thin[,1],xlim=x.lim,pch=".",ylim= c(min(GAMMA.thin,gam.true*0.9),max(GAMMA.thin,gam.true*1.1)))
  abline(h = gam.true[1],col=2)
  plot(DELTA.thin,xlim=x.lim,pch=".")
  abline(h = DELTA.true,col=2)

  plot(ALP.thin/BETA.thin,ylim=c(min(ALP.thin/BETA.thin),max(ALP.thin/BETA.thin)),pch=".")
  title(paste("mean K =", round(mean(ALP.thin/BETA.thin),5),":::: True K =",round(A.true/B.true,3)),line=0.5)
  abline(h= A.true/B.true,col=2)

  plot(1/sqrt(ALP.thin),main=paste("CV k",round(mean(1/sqrt(ALP.thin)),3),":::: True = ",1/sqrt(A.true)),pch=".")
  #	main= paste("mean CV =", round(sqrt(ALP.out),3),":::: True CV =",round(1/sqrt(k.A),3)))
  abline(h= 1/sqrt(A.true),col=2)

  plot(sqrt(ALP.thin/BETA.thin^2),main=paste("SD k",round(mean(sqrt(ALP.thin/BETA.thin^2)),4),":::: True = ",round(sqrt(A.true/B.true^2),3)),pch=".")
  #		main= paste("mean SD k =", round(sqrt(ALP.thin[0.9*i:i]/BETA.thin[0.9*i:i]^2)),3),":::: True SD k =",round(k.A/k.B^2,3))
  abline(h= sqrt(A.true/B.true^2),col=2)

  #### PLOT ACF of MCMC
  par(mfrow=c(3,2),mar=c(4,5,3,0.5))
  x.lim=c(0,length(THIN))

  acf(X0.thin,main="x0")
  acf(ALP.thin,main=expression(alpha))
  acf(BETA.thin,main=expression(beta))
  acf(GAMMA.thin[,1],main=expression(gamma))
  if(max(PSI.thin)!=psi.start){
    acf(PSI.thin,main=expression(psi))
  }
  acf(DELTA.thin,main=expression(delta))

  ##### HISTOGRAM OF MARGINAL DISTRIBUTIONS for non-gamma parameters
  par(mfrow=c(3,3),mar=c(4,5,2,0.5))

  hist(X0.thin, main="x0",xlab="",breaks=20)
  abline(v=x0.true,col=2,lwd=1.5)
  abline(v=mean(X0.thin),lty=2,col=2,lwd=1.5)
  abline(v=median(X0.thin),col=4,lwd=1.5)
  hist(ALP.thin, main=expression(alpha),xlab="",breaks=20)
  abline(v=A.true,col=2,lwd=1.5)
  abline(v=mean(ALP.thin),lty=2,col=2,lwd=1.5)
  abline(v=median(ALP.thin),col=4,lwd=1.5)
  hist(BETA.thin, main=expression(beta),xlab="",breaks=20)
  abline(v=B.true,col=2,lwd=1.5)
  abline(v=mean(BETA.thin),lty=2,col=2,lwd=1.5)
  abline(v=median(BETA.thin),col=4,lwd=1.5)
  hist(DELTA.thin, main=expression(delta),xlab="",breaks=20)
  abline(v=DELTA.true,col=2,lwd=1.5)
  abline(v=mean(DELTA.thin),lty=2,col=2,lwd=1.5)
  abline(v=median(DELTA.thin),col=4,lwd=1.5)
  hist(PSI.thin, main=expression(psi),xlab="",breaks=20)
  abline(v=psi.true,col=2,lwd=1.5)
  abline(v=mean(PSI.thin),lty=2,col=2,lwd=1.5)
  abline(v=median(PSI.thin),col=4,lwd=1.5)
  hist(k.bar.thin, main="mean(k)" ,breaks=20)
  abline(v=mean.k,col=2,lwd=1.5)
  abline(v=mean(k.bar.thin),lty=2,col=2,lwd=1.5)
  abline(v=median(k.bar.thin),col=4,lwd=1.5)
  hist(k.CV.thin, main="CV(k)" ,breaks=20)
  abline(v=1/ sqrt(A.true),col=2,lwd=1.5)
  abline(v=mean(k.CV.thin),lty=2,col=2,lwd=1.5)
  abline(v=median(k.CV.thin),col=4,lwd=1.5)

  ##### HISTOGRAM OF MARGINAL DISTRIBUTIONS for gamma parameters
  par(mfrow=c(3,3),mar=c(4,5,2,0.5))
  for(i in 1:(OBS.high -1)){
    hist(GAMMA.thin[,i], main=paste(expression(gamma),i),xlab="",breaks=20)
    abline(v=gam.all[i],col=2,lwd=1.5)
    abline(v=mean(GAMMA.thin[,i]),lty=2,col=2,lwd=1.5)
    abline(v=median(GAMMA.thin[,i]),col=4,lwd=1.5)
  }

  ### PLOT MEANS....
  par(mfrow=c(1,1))
  y.lim <- c(0,100)
  plot(MU.thin[,2],ylim=y.lim,col="darkgreen",pch=".",ylab="Length at age Mean")
  abline(h=REAL.mean[2],lwd=1.5,col=3,pch=".")
  par(new=T)
  plot(MU.thin[,3],ylim=y.lim,col=2,pch=".",ylab = "")
  abline(h=REAL.mean[3],lwd=1.5,col=2,pch=".")
  par(new=T)
  plot(MU.thin[,4],ylim=y.lim,col=4,pch=".",ylab = "")
  abline(h=REAL.mean[4],lwd=1.5,col=4)
  par(new=T)
  plot(MU.thin[,5],ylim=y.lim,col=5,pch=".",ylab = "")
  abline(h=REAL.mean[5],lwd=1.5,col=5)
  par(new=T)
  plot(MU.thin[,6],ylim=y.lim,col=6,pch=".",ylab = "")
  abline(h=REAL.mean[6],lwd=1.5,col=6)
  par(new=T)

  if(7 <= OBS.high){
    plot(MU.thin[,7],ylim=y.lim,col=7,pch=".",ylab = "")
    abline(h=REAL.mean[7],lwd=1.5,col=7)
    par(new=T)
  }
  if(8 <= OBS.high){
    plot(MU.thin[,8],ylim=y.lim,col=8,pch=".",ylab = "")
    abline(h=REAL.mean[8],lwd=1.5,col=8)
  }
  #		par(new=T)
  legend( x =0, y=min(MU.thin[,2:dim(MU.thin)[2]]),
          lwd =2, col= c("darkgreen",2,4,3,5,6,7,8),cex = 0.7,
          c("Age 1", "Age 2", "Age 3", "Age 4", "Age 5", "Age 6", "Age 7", "Age 8"))

  ####  PLOT VARIANCES....
  y.lim <- c(0,REAL.sd[OBS.high]*1.1)
  plot(SIG.thin[,2],ylim=y.lim,col="darkgreen",pch=".",ylab="Length at age Variance")
  abline(h=REAL.sd[2],lwd=1.5,col="darkgreen")
  par(new=T)
  plot(SIG.thin[,3],ylim=y.lim,col=2,pch=".",ylab = "")
  abline(h=REAL.sd[3],lwd=1.5,col=2)
  par(new=T)
  plot(SIG.thin[,4],ylim=y.lim,col=4,pch=".",ylab = "")
  abline(h=REAL.sd[4],lwd=1.5,col=4)
  par(new=T)
  plot(SIG.thin[,5],ylim=y.lim,col=3,pch=".",ylab = "")
  abline(h=REAL.sd[5],lwd=1.5,col=3)
  par(new=T)
  plot(SIG.thin[,6],ylim=y.lim,col=5,pch=".",ylab = "")
  abline(h=REAL.sd[6],lwd=1.5,col=5)
  par(new=T)
  if(7 <= OBS.high){
    plot(SIG.thin[,7],ylim=y.lim,col=6,pch=".",ylab = "")
    abline(h=REAL.sd[7],lwd=1.5,col=6)
    par(new=T)
  }
  if(8 <= OBS.high){
    plot(SIG.thin[,8],ylim=y.lim,col=7,pch=".",ylab = "")
    abline(h=REAL.sd[8],lwd=1.5,col=7)
    par(new=T)
  }
  if(9 <= OBS.high){
    plot(SIG.thin[,9],ylim=y.lim,col=8,pch=".",ylab = "")
    abline(h=REAL.sd[9],lwd=1.5,col=8)
  }
  legend( x =0, y=min(SIG.thin[,2:dim(SIG.thin)[2]]),
          lwd =2, col= c("darkgreen",2,4,3,5,6,7,8),cex = 0.7,
          c("Age 1", "Age 2", "Age 3", "Age 4", "Age 5", "Age 6", "Age 7", "Age 8"))

  ##############################################################################################
  ## HISTOGRAMS OF LENGTHS of Sampled INDIVIDUALS and estimated lengths
  ##############################################################################################

  all.MU	<- colMeans(MU.thin)
  all.SIG	<- colMeans(SIG.thin)

  par(mfrow=c(3,4))

  for(i in 1:OBS.high){

    Z <- DAT.1$X[DAT.1$date==i]
    if(samp.numbs[i] > 0){
      x.lim	<- c(min(Z)*0.99, max(Z)*1.01)
      hist(Z,xlim=x.lim,main=paste("Length of age ",i-1),breaks=20)
      par(new=T)
      plot(x = seq(x.lim[1], x.lim[2],length=100),
           y = dgamma( seq(x.lim[1], x.lim[2],length=100),
                       all.MU[i]^2 / all.SIG[i]^2, all.MU[i] / all.SIG[i]^2),
           xlim=x.lim,type="l",col=2,lwd=2,
           axes=F,xlab="",ylab="")
    }
    if(samp.numbs[i]==0){
      plot(1:5,1:5,axes=F,col="white",xlab="",ylab="")
    }
  }

  #### MAKE A PLOT OF 10 GROWTH TRAJECTORIES.
  par(mfrow=c(1,1))
  y.lim=c(0,max(DAT.fixed.k[1:OBS.high,]))
  for(i in 1:20){
    plot(y=DAT.fixed.k[i,1:OBS.high],x=(0:(dim(MU.thin)[2]-1)),ylim=y.lim,type="l",
         ylab="",xlab="",axes=F)
    par(new=T)
  }
  ## Add the mean and sd at each age.
  plot(x= (0:(dim(MU.thin)[2]-1)) , y = colMeans(MU.thin), col = 2,lwd=2.5,ylim=y.lim,axes=T,
       pch=21, bg=2, cex= 1.5,xlab="Age",ylab="Length")

  ##############################################################################################
  #### PLOT SKEW
  ##############################################################################################

  par(mfrow=c(1,1))
  skew.data<- vector(mode="numeric",length=length(samp.numbs))
  for(i in 1: (length(samp.numbs))){
    skew.data[i]	<-	skew(DAT.fixed.k[,i])
  }#

  temp 	<- max(c(abs(min(skew.data)),abs(max(skew.data))))
  y.lim = c(-temp*1.1,temp*1.1)
  x.lim = c(0,(length(samp.numbs) -1))

  plot(y=skew.data, x=(0:(length(samp.numbs)-1)),
       xlim = x.lim, ylim=y.lim,xlab="TIME",ylab="SKEW",lwd=2)
  par(new=T)
  plot(y=skew.data[samp.numbs>0], x=(0:(length(samp.numbs)-1))[samp.numbs>0],
       xlim=x.lim, ylim=y.lim, col = 2,xlab="",ylab="",lwd=2)
  abline(h=0,lty=2,lwd=1.5)
  title("Red indicates observations")
  text(samp.numbs,y=y.lim[1],x=(0:(length(samp.numbs)-1)),cex=0.8)

  par(mfrow=c(1,1))
  plot(1:5,1:5,col="white",axes=F,xlab="",ylab="")
  text(x=2, y=5, paste("Accept ALPHA =",round(accept.alpha/ NUMB,3)))
  text(x=2, y=4.5, paste("Accept Beta =",round(accept.beta/ NUMB,3)))
  text(x=2, y=4, paste("Accept Gamma =",round(accept.gamma/ NUMB,3)))
  text(x=2, y=3.5, paste("Accept Delta =",round(accept.delta/ NUMB,3)))
  text(x=2, y=3, paste("Accept Psi =",round(accept.psi/ NUMB,3)))

  dev.off()


  ##################################################################################################
  ##################################################################################################
  ### WRITE OUTPUT TO FILE...
  ##################################################################################################
  ##################################################################################################

  PRINT.TO.FILE	<- cbind( X0.thin, ALP.thin, BETA.thin,
                          DELTA.thin, PSI.thin, k.bar.thin, k.CV.thin)
  colnames(PRINT.TO.FILE)	<- c("X0","ALPHA","BETA","DELTA","PSI","k.bar","k.CV")

  colnames(GAMMA.thin)	<- paste("GAMMA",seq(0,OBS.high-2,1),sep=".")
  colnames(MU.thin)		<- paste("MU",seq(0,OBS.high-1,1),sep=".")
  colnames(SIG.thin)		<- paste("SIG",seq(0,OBS.high-1,1),sep=".")

  PRINT.TO.FILE <- cbind(PRINT.TO.FILE,GAMMA.thin,MU.thin,SIG.thin)


  setwd("/Users/ole/Documents/Von B Growth/OUTPUT/Observational samples/CSV data")

  write.csv(PRINT.TO.FILE,file=paste(NAME,"growth.observational.fits.csv"),row.names=F)



  ##################################### ADD THIS STUFF

}	#### END SCENARIO LOOP.





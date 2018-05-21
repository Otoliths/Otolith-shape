rm = list(ls())
#### INPUTS - k = constant.
# Process
CV <-	0.15
A <-	1 / CV^2
B =	3*A
T =	1
q = 85
x0 = 10

mean.k <- A/B
sd.k <- A/B^2

k.ALPHA <- A
k.BETA <- B

k.ALPHA / k.BETA	### THIS IS THE MEAN k
1/sqrt(k.ALPHA)		### THIS IS THE CV of k

N = 1000000
#####################################################################################################
#####################################################################################################
#####################################################################################################################################
### MAKE Q CONSTANT, k variabile among individuals
#####################################################################################################################################
#####################################################################################################
#####################################################################################################
### SIMULATE K fixed for individuals
QQ = 10

DAT.fixed.k	<-	data.frame(matrix(0,N,QQ))
DAT.rand.k	<-	data.frame(matrix(0,N,QQ))

N.start	<-	10
DAT.fixed.k[,1]	<-	N.start
DAT.rand.k[,1]	<-	N.start

k	<- rgamma(N,shape=k.ALPHA,rate=k.BETA)

for(i in 2:QQ){	### LOOP OVER TIME

  #q	<- rnorm(1,mean = 85, sd = 10)
  DAT.fixed.k[,i]	<-	(q/k)*(1-exp(-k)) + exp(-k)*DAT.fixed.k[,i-1]
}

#for(i in 2:QQ){	### LOOP OVER TIME
#	k	<- rgamma(N,shape=k.ALPHA,rate=k.BETA)
#	DAT.rand.k[,i]	<-	(q/k)*(1-exp(-k)) + exp(-k)*DAT.rand.k[,i-1]
#}

#calcualte some summary statistics:

SD.f	<- sd(DAT.fixed.k)
SD.r	<- sd(DAT.rand.k)

mean.f	<- mean(DAT.fixed.k)
mean.r	<- mean(DAT.rand.k)

CV.f	<- SD.f / mean.f
CV.r	<- SD.r / mean.r
#############################################################################################
# Write down expectation from analytic expression.
#############################################################################################
## EXPECTATION
T=1:5

p.1 = k.BETA * q  / (k.ALPHA -1)
p.2 = (k.BETA^k.ALPHA) * x0 / ( T + k.BETA)^(k.ALPHA)
p.3 = (q * (k.BETA^k.ALPHA)) / ((k.ALPHA-1)*( T + k.BETA)^(k.ALPHA-1))

expect.size <- p.1 + p.2 - p.3

plot(mean.f[2:QQ],expect.size[1:(QQ-1)])

T	<-	1:10
## VARIANCE

p.1 <-	((q^2)*(B^2) / (A-1)) * ((1/(A-2)) - (1/(A-1)))

p.2	<-	(2*(q^2)*(B^A) / ((A-1)*(T+B)^(A-2)) ) * (
  (B/ ((A-1)*(T+B))) -
    (1/(A-2))
)

p.3	<-	((q^2)*(B^A) / (A-1)) *
  ((1 / ((A-2)*((2*T + B)^(A-2))) ) -
     ((B^A) / ((A-1) * (T + B)^(2*A-2) )) )


p.4 <-	(2*q*B^A*x0 / (A-1)) * ( (1 / (T + B)^(A-1)) - ( 1 / (2 * T + B)^(A-1) ) )

p.5	<-	(2*q*B^(A+1)*x0 /((A-1)*(T + B)^(A))) * (( B^(A-1)/(T+B)^(A-1)) - 1)


p.6	<-	x0^2*B^A * ( (1 / (2*T + B)^(A)) - ( B^(A) /(T + B)^(2*A) ) )

p.1; p.2; p.3; p.4; p.5; p.6

sqrt(p.1+p.2+p.3+p.4+p.5+p.6)
SD.f


















##################################
#### PLOT RANDOM K EACH YEAR
##################################
setwd("/Users/ole/Documents/Herring/_Herring UCSC/State Space Growth/Simulation 2-2011")

pdf(file=paste("Q.Fixed; K.CV =",CV.k,"pdf",sep="."),onefile=TRUE,width=7,height=7)

par(mfrow=c(3,4),mar=c(3,4,3,0.5))
hist(k,main="")
title(paste("K =", mean.k,"; CV.k = ", CV.k,"; q = 85"),line=2,cex.main=0.8)
title(paste("Start size = 10mm"),line=0.5)
title( xlab="k",line=2)

for(i in 1:11){
  hist(DAT.fixed.k[,i+1],main="",xlab="")
  if(i==1){
    title(paste("Year =", i+1),line=0.5)
    title("K FIXED FOR IND",line=2)
  }
  if(i>1){
    title(paste("Year =", i+1),line=2)
    title(paste("mean=",round(mean.f[i+1],0),"; sd=",round(SD.f[i+1],0),"; CV=",round(SD.f[i+1]/mean.f[i+1],3)),line=0.5,cex.main=0.7)
    title(xlab="Length(mm)",line=2)
  }
}


par(mfrow=c(3,4),mar=c(3,4,3,0.5))
hist(k,main="")
title(paste("K =", mean.k,"; CV.k = ", CV.k,"; q = 85"),line=2,cex.main=0.8)
title(paste("Start size = 10mm"),line=0.5)
title( xlab="k",line=2)

for(i in 1:seq(5,40,by=5)){
  hist(DAT.fixed.k[,i+1],main="",xlab="")
  if(i==1){
    title(paste("Year =", i+1),line=0.5)
    title("K FIXED FOR IND",line=2)
  }
  if(i>1){
    title(paste("Year =", i+1),line=2)
    title(paste("mean=",round(mean.f[i+1],0),"; sd=",round(SD.f[i+1],0),"; CV=",round(SD.f[i+1]/mean.f[i+1],3)),line=0.5,cex.main=0.7)
    title(xlab="Length(mm)",line=2)
  }
}



# Write down expectation from analytic expression.

## EXPECTATION

p.1 = k.BETA * q  / (k.ALPHA -1)
p.2 = (k.BETA^k.ALPHA) * x0 / ( T + k.BETA)^(k.ALPHA)
p.3 = (q * (k.BETA^k.ALPHA)) / ((k.ALPHA-1)*( T + k.BETA)^(k.ALPHA-1))

expect.size <- p.1 + p.2 - p.3

plot(mean.f[2:11],expect.size[1:10])


B	<- k.BETA
A	<- k.ALPHA
T	<- 1:2


## VARIANCE

p.1 <-	((q^2)*(B^2) / (A-1)) * ((1/(A-2)) - (1/(A-1)))

p.2	<-	(2*(q^2)*(B^A) / ((A-1)*(T+B)^(A-2)) ) * (
  (B/ ((A-1)*(T+B))) -
    (1/(A-2))
)

p.3	<-	((q^2)*(B^A) / (A-1)) *
  ((1 / ((A-2)*((2*T + B)^(A-2))) ) -
     ((B^A) / ((A-1) * (T + B)^(2*A-2) )) )


p.4 <-	(2*q*B^A*x0 / (A-1)) * ( (1 / (T + B)^(A-1)) - ( 1 / (2 * T + B)^(A-1) ) )

p.5	<-	(2*q*B^(A+1)*x0 /((A-1)*(T + B)^(A))) * (( B^(A-1)/(T+B)^(A-1)) - 1)


p.6	<-	x0^2*B^A * ( (1 / (2*T + B)^(A)) - ( B^(A) /(T + B)^(2*A) ) )

p.1
p.2
p.3
p.4
p.5
p.6

sqrt(p.1+p.2+p.3+p.4+p.5+p.6)
SD.f














##################################
#### PLOT RANDOM K EACH YEAR
##################################
par(mfrow=c(3,4),mar=c(3,4,3,0.5))

hist(k,main="")
title(paste("K =", mean.k,"; CV.k = ", CV.k,"; q = 85"),line=2,cex.main=0.8)
title(paste("Start size = 10mm"),line=0.5)
title( xlab="k",line=2)

for(i in 1:11){
  hist(DAT.rand.k[,i+1],main="",xlab="")
  if(i==1){
    title(paste("Year =", i+1),line=0.5)
    title("K RAND FOR IND",line=2)
  }
  if(i>1){
    title(paste("Year =", i+1),line=2)
    title(paste("mean=",round(mean.r[i+1],0),"; sd=",round(SD.r[i+1],0),"; CV=",round(SD.r[i+1]/mean.r[i+1],3)),line=0.5,cex.main=0.7)
    title(xlab="Length(mm)",line=2)
  }
}

dev.off()



#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
### MAKE K CONSTANT, Q variabile among individuals
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
rm(list=ls())

#### INPUTS - k = constant.
# Process

mean.k	<-	0.35

mean.q	= 85
CV.q	<-	0.5
sd.q	<-	mean.q * CV.q

q.ALPHA	<-	mean.q^2 / sd.q^2
q.BETA	<-	mean.q / sd.q^2

N=1000000

#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
### SIMULATE K fixed for individuals

DAT.fixed.q	<-	data.frame(matrix(0,N,12))
DAT.rand.q	<-	data.frame(matrix(0,N,12))

N.start	<-	10
DAT.fixed.q[,1]	<-	N.start
DAT.rand.q[,1]	<-	N.start

q	<- rgamma(N,shape=q.ALPHA,rate=q.BETA)
k	<- mean.k

for(i in 2:12){	### LOOP OVER TIME
  DAT.fixed.q[,i]	<-	(q/k)*(1-exp(-k)) + exp(-k)*DAT.fixed.q[,i-1]
}

for(i in 2:12){	### LOOP OVER TIME
  q	<- rgamma(N,shape=q.ALPHA,rate=q.BETA)
  DAT.rand.q[,i]	<-	(q/k)*(1-exp(-k)) + exp(-k)*DAT.rand.q[,i-1]
}

#calcualte some summary statistics:

SD.f	<- sd(DAT.fixed.q)
SD.r	<- sd(DAT.rand.q)

mean.f	<- mean(DAT.fixed.q)
mean.r	<- mean(DAT.rand.q)

CV.f	<- SD.f / mean.f
CV.r	<- SD.r / mean.r

##################################
#### PLOT RANDOM Q EACH YEAR
##################################

setwd("/Users/ole/Documents/Herring/_Herring UCSC/State Space Growth/Simulation 2-2011")

pdf(file=paste("K.Fixed; q.CV =",CV.q,"pdf",sep="."),onefile=TRUE,width=7,height=7)

par(mfrow=c(3,4),mar=c(3,4,3,0.5))
hist(q,main="")
title(paste("Q =", mean.q,"; CV.q = ", CV.q,"; k=", mean.k),line=2,cex.main=0.8)
title(paste("Start size = 10mm"),line=0.5)
title( xlab="q",line=2)

for(i in 1:11){
  hist(DAT.fixed.q[,i+1],main="",xlab="")
  if(i==1){
    title(paste("Year =", i+1),line=0.5)
    title("Q FIXED FOR IND",line=2)
  }
  if(i>1){
    title(paste("Year =", i+1),line=2)
    title(paste("mean=",round(mean.f[i+1],0),"; sd=",round(SD.f[i+1],0),"; CV=",round(SD.f[i+1]/mean.f[i+1],3)),line=0.5,cex.main=0.7)
    title(xlab="Length(mm)",line=2)
  }
}
##################################
#### PLOT RANDOM Q EACH YEAR
##################################
par(mfrow=c(3,4),mar=c(3,4,3,0.5))

hist(q,main="")
title(paste("Q =", mean.q,"; CV.q = ", CV.q,"; k=", mean.k),line=2,cex.main=0.8)
title(paste("Start size = 10mm"),line=0.5)
title( xlab="k",line=2)

for(i in 1:11){
  hist(DAT.rand.q[,i+1],main="",xlab="")
  if(i==1){
    title(paste("Year =", i+1),line=0.5)
    title("Q RAND FOR IND",line=2)
  }
  if(i>1){
    title(paste("Year =", i+1),line=2)
    title(paste("mean=",round(mean.r[i+1],0),"; sd=",round(SD.r[i+1],0),"; CV=",round(SD.r[i+1]/mean.r[i+1],3)),line=0.5,cex.main=0.7)
    title(xlab="Length(mm)",line=2)
  }
}

dev.off()

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
### MAKE K AND Q variabile among individuals
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
rm(list=ls())
### MAKE k and q independent random variables.

KK	<- c(0.01,0.05,0.10,0.2,0.5)
QQ	<- c(0.01,0.05,0.10,0.2,0.5)

### LOOP OVER POSSIBLE CV COMBINATIONS
for(j in 1:length(KK)){
  for(l in 1: length(QQ)){
    ####

    mean.k	<-	0.35
    CV.k	<-	KK[j]
    sd.k	<-	mean.k * CV.k

    k.ALPHA	<-	mean.k^2 / sd.k^2
    k.BETA	<-	mean.k / sd.k^2

    mean.q	= 85
    CV.q	<-	QQ[l]
    sd.q	<-	mean.q * CV.q

    q.ALPHA	<-	mean.q^2 / sd.q^2
    q.BETA	<-	mean.q / sd.q^2

    N = 1000000
    #####################################################################################################
    #####################################################################################################
    #####################################################################################################
    #####################################################################################################

    DAT.fixed	<-	data.frame(matrix(0,N,12))
    DAT.rand	<-	data.frame(matrix(0,N,12))

    N.start	<-	10
    DAT.fixed[,1]	<-	N.start
    DAT.rand[,1]	<-	N.start

    q	<- rgamma(N,shape=q.ALPHA,rate=q.BETA)
    k	<- rgamma(N,shape=k.ALPHA,rate=k.BETA)

    ### SIMULATE K and q fixed for individuals
    for(i in 2:12){	### LOOP OVER TIME
      DAT.fixed[,i]	<-	(q/k)*(1-exp(-k)) + exp(-k)*DAT.fixed[,i-1]
    }

    ### SIMULATE K and q random each year for individuals
    for(i in 2:12){	### LOOP OVER TIME
      k	<- rgamma(N,shape=k.ALPHA,rate=k.BETA)
      q	<- rgamma(N,shape=q.ALPHA,rate=q.BETA)
      DAT.rand[,i]	<-	(q/k)*(1-exp(-k)) + exp(-k)*DAT.rand[,i-1]
    }

    #calcualte some summary statistics:

    SD.f	<- sd(DAT.fixed)
    SD.r	<- sd(DAT.rand)

    mean.f	<- mean(DAT.fixed)
    mean.r	<- mean(DAT.rand)

    CV.f	<- SD.f / mean.f
    CV.r	<- SD.r / mean.r

    ##################################
    #### PLOT RANDOM Q and RANDOM K EACH INDIVIDUAL
    ##################################

    setwd("/Users/ole/Documents/Herring/_Herring UCSC/State Space Growth/Simulation 2-2011/_K-Q both vary")

    pdf(file=paste("K.CV=",CV.k,"; q.CV =",CV.q,"pdf",sep="."),onefile=TRUE,width=7,height=7)

    par(mfrow=c(3,4),mar=c(3,4,3,0.5))
    hist(q,main="")
    title(paste("CV.q = ", CV.q),line=2,cex.main=0.8)
    title(paste("Start size = 10mm"),line=0.5)
    title( xlab="q",line=2)

    hist(k,main="")
    title(paste("CV.k = ", CV.k),line=2,cex.main=0.8)
    title(paste("Start size = 10mm"),line=0.5)
    title( xlab="k",line=2)

    for(i in 1:10){
      hist(DAT.fixed[,i+1],main="",xlab="")
      if(i==1){
        title(paste("Year =", i+1),line=0.5)
        title("Q and K FIXED FOR IND",line=2,cex.main=0.9)
      }
      if(i>1){
        title(paste("Year =", i+1),line=2)
        title(paste("mean=",round(mean.f[i+1],0),"; sd=",round(SD.f[i+1],0),"; CV=",round(SD.f[i+1]/mean.f[i+1],3)),line=0.5,cex.main=0.7)
        title(xlab="Length(mm)",line=2)
      }
    }
    ##################################
    #### PLOT RANDOM Q and RANDOM K EACH YEAR
    ##################################
    par(mfrow=c(3,4),mar=c(3,4,3,0.5))

    hist(q,main="")
    title(paste("CV.q = ", CV.q),line=2,cex.main=0.8)
    title(paste("Start size = 10mm"),line=0.5)
    title( xlab="q",line=2)

    hist(k,main="")
    title(paste("CV.k = ", CV.k),line=2,cex.main=0.8)
    title(paste("Start size = 10mm"),line=0.5)
    title( xlab="k",line=2)

    for(i in 1:10){
      hist(DAT.rand[,i+1],main="",xlab="")
      if(i==1){
        title(paste("Year =", i+1),line=0.5)
        title("Q and K RAND FOR IND",line=2,cex.main=0.9)
      }
      if(i>1){
        title(paste("Year =", i+1),line=2)
        title(paste("mean=",round(mean.r[i+1],0),"; sd=",round(SD.r[i+1],0),"; CV=",round(SD.r[i+1]/mean.r[i+1],3)),line=0.5,cex.main=0.7)
        title(xlab="Length(mm)",line=2)
      }
    }

    dev.off()

    #### END LOOP ACROSS SCENARIOS.
  }
}








######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################



rm=list(ls())
#### INPUTS - k = constant.
# Process

k	<-	0.35
q	<-  85

N=100000

DAT.fixed	<-	data.frame(matrix(0,N,12))
DAT.rand	<-	data.frame(matrix(0,N,12))

N.start	<-	10
DAT.fixed[,1]	<-	N.start
DAT.rand[,1]	<-	N.start

#q	<- rgamma(N,shape=q.ALPHA,rate=q.BETA)
#k	<- rgamma(N,shape=k.ALPHA,rate=k.BETA)

### SIMULATE K and q fixed for individuals
for(i in 2:12){	### LOOP OVER TIME
  DAT.fixed[,i]	<-	((q/k)*(1-exp(-k)) + exp(-k)*DAT.fixed[,i-1]) + (rnorm(N,0,1*(DAT.fixed[,i-1]))*((q/k)-DAT.fixed[i-1]))
}


SD.f	<- sd(DAT.fixed)
#SD.r	<- sd(DAT.rand)

mean.f	<- mean(DAT.fixed)
#mean.r	<- mean(DAT.rand)

CV.f	<- SD.f / mean.f
#CV.r	<- SD.r / mean.r


##################################

setwd("/Users/ole/Documents/Herring/_Herring UCSC/State Space Growth/Simulation 2-2011")

#pdf(file=paste("K.Fixed; q.CV =",CV.q,"pdf",sep="."),onefile=TRUE,width=7,height=7)

par(mfrow=c(3,4),mar=c(3,4,3,0.5))
hist(q,main="")
title(paste("Q =", q,"; CV.q = ", 0,"; k=", k),line=2,cex.main=0.8)
title(paste("Start size = 10mm"),line=0.5)
title( xlab="q",line=2)

for(i in 1:11){
  hist(DAT.fixed[,i+1],main="",xlab="")
  if(i==1){
    title(paste("Year =", i+1),line=0.5)
    title("Q+K FIXED FOR IND",line=2)
  }
  if(i>1){
    title(paste("Year =", i+1),line=2)
    title(paste("mean=",round(mean.f[i+1],0),"; sd=",round(SD.f[i+1],0),"; CV=",round(SD.f[i+1]/mean.f[i+1],3)),line=0.5,cex.main=0.7)
    title(xlab="Length(mm)",line=2)
  }
}
##################################
quartz()
hist(DAT.fixed[,12]-DAT.fixed[,11])

hist(1000/DAT.fixed[,i-1])


hist(exp(rnorm(N,0,0.01)))

























mean.k	<-	0.35
CV.k	<-	0.01
sd.k	<-	mean.k * CV.k

k.ALPHA	<-	mean.k^2 / sd.k^2
k.BETA	<-	mean.k / sd.k^2

k.ALPHA / k.BETA	### THIS IS THE MEAN k
1/sqrt(k.ALPHA)		### THIS IS THE CV of k

mean.q	= 85
CV.q	<-	0.01
sd.q	<-	mean.q * CV.q

q.ALPHA	<-	mean.q^2 / sd.q^2
q.BETA	<-	mean.q / sd.q^2


















# Observations
alpha	=	-20
bet		=	0.2

Y	=	1/(1+exp(-alpha - bet*(80:120)))

plot(x=80:120,y=Y)

N = 1000000

DAT	<-	data.frame(matrix(0,N,10))



var(DAT)
pred.var	<-	exp(-2*k) * sd(DAT[,1])^2 + tau ^2; print(pred.var)
pred.var.2	<-	exp(-2*k) * sd(DAT[,2])^2 + tau ^2; print(pred.var.2)





OBS	<-	data.frame(matrix(0,N,4))

OBS[,1]	<-
  runif(N)

hist(DAT[,1]	/(1+exp(-alpha - bet*DAT[,1])))

A <- runif(N)
hist(DAT[which(1/(1+exp(-alpha - bet*DAT[,1])) > A),1])












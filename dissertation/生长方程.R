library("fishmethods")
library("FSA")
## Make an example age-length key
data(WR79)
WR.age <- subset(WR79, !is.na(age))      # isolate the age sample
WR.age$LCat <- lencat(WR.age$len,w=5)    # add length intervals (width=5)
raw <- xtabs(~LCat+age,data=WR.age)      # create age-length key
( WR.key <- prop.table(raw, margin=1) )



## Various visualizations of the age-length key
alkPlot(WR.key,"barplot")
alkPlot(WR.key,"barplot",pal="gray")
alkPlot(WR.key,"barplot",showLegend=TRUE)
alkPlot(WR.key,"area")
alkPlot(WR.key,"area",showLegend=TRUE)
alkPlot(WR.key,"area",pal="gray")
alkPlot(WR.key,"lines")
alkPlot(WR.key,"lines",pal="gray")
alkPlot(WR.key,"lines",showLegend=TRUE)
alkPlot(WR.key,"splines")
alkPlot(WR.key,"splines",span=0.2)
alkPlot(WR.key,"splines",pal="gray",showLegend=TRUE)
alkPlot(WR.key,"bubble")
alkPlot(WR.key,"bubble",grid=FALSE)
alkPlot(WR.key,"bubble",grid="blue")
alkPlot(WR.key,"bubble",grid=rgb(0,0,0,0.2),col=rgb(0,0,0,0.5))



data(pinfish)
growth(intype = 1,unit = 1,size = pinfish$sl,age = pinfish$age,
       calctype = 1,wgtby=1,error = 1,Sinf = 200,K = 0.3,t0 = -1)

data(trout)
growhamp(L1=trout$L1,L2=trout$L2,TAL=trout$dt,models=c(1,2,3),
         method=c("Nelder-Mead","Nelder-Mead","L-BFGS-B"),
         varcov=c(TRUE,TRUE,TRUE),
         Linf=list(startLinf=650,lowerLinf=400,upperLinf=800),
         K=list(startK=0.30,lowerK=0.01,upperK=1),
         sigma2_error=list(startsigma2=100,lowersigma2=0.1,uppersigma2=10000),
         sigma2_Linf=list(startsigma2=100,lowersigma2=0.1,uppersigma2=100000),
         sigma2_K=list(startsigma2=0.5,lowersigma2=1e-8,uppersigma2=10))


data(bonito)
temp <- bonito[c(bonito$T2-bonito$T1)>0,]
growthResid(0.19,97.5,lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2 - temp$T1),graph=1)
growthTraject(0.19,97.5,lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2 - temp$T1))

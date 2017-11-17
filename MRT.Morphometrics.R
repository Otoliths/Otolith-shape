library(ade4)
library(vegan)
library(mvpart)
library(MVPARTwrap)
library(FactoMineR)
library(MASS)
source("D:\\Project\\Data\\MRTG.R")
tt=read.csv("D:\\Project\\Data\\Morphometrics.csv",
            header = T)
s.auto <- mvpart(data.matrix(tt[,4:13])~Average.river.gradient+Altitude+Annual.average.temperature+Latitude,tt,margin=0.1,cp=0,xv="pick",xval = nrow(tt),xvmult = 100,which = 4)
s.auto$cptable
rsqrpart(s.auto)
(R2<-abs(diff(s.auto$cptable[,3])))
importance<-VarImp(s.auto, pretty=T)
importance$Relative
rpart.pca(s.auto)
rpart.pca(s.auto,wgt.ave=TRUE,interact=TRUE)
res<-MRT(s.auto,10)
plot(res, Cex=1,widthtree=10, heighttree=10,lwd=2)
importance<-VarImp(s.auto, pretty=F)
importance$Relative

########################
df11<-tt[,4:13]
df1<- decostand(df11, "hellinger")
names(df1)
df2<-tt[,14:17]
names(df2)
tab2 <- data.frame(df1,df2)
dim(tab2)
(num <- c(ncol(df1), ncol(df2)))

# Compute the MFA with multiple plots
# graphics.off()  # Close the previous graphic windows
t2.mfa <- MFA(tab2, group=num, type=c("c","s"), ncp=2,
              name.group=c("Otolith morphological parameters ","Multi-scale environmental factors"))
t2.mfa
#customize image
plot(t2.mfa, choix="var", habillage="group",axes = 1:2,cex=0.8,shadow=T)
(rvp <- t2.mfa$group$RV)
rvp[1,2] <- coeffRV(df11,scale(df2))$p.value
rvp[1,2]
##################
MCC(s.auto,standard = F)
MCC(s.auto,weight = F)
Interaction <- MCC(s.auto ,weight=T)
Interaction$interact
s.auto$frame

#result
rmarkdown::render("MRT.Morphometrics.R",
                 "pdf_document")

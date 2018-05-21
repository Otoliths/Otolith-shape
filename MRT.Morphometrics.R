##FIGURE 4 Regression analysis of 493 otolith shapes from Schizothorax nukiangensis in the Nu-Salween River
library(ade4)
library(vegan)
library(mvpart)
library(MVPARTwrap)
library(FactoMineR)
library(MASS)
source("Data/MRTG.r")
tt <- read.csv("Data/Morphometrics.csv",
              header = T)
s.auto <- mvpart(data.matrix(tt[,4:7])~Average.river.gradient+Altitude+Annual.average.temperature+Latitude,tt,margin = 0.1,cp = 0,xv = "pick",xval = nrow(tt),xvmult = 100,which = 4)
s.auto$cptable
rsqrpart(s.auto)
(R2 <- abs(diff(s.auto$cptable[,3])))
importance <- VarImp(s.auto, pretty = T)
importance$Relative
rpart.pca(s.auto)
rpart.pca(s.auto,wgt.ave = TRUE,interact = TRUE)
res <- MRT(s.auto,10)
##summary.MRT(res)
plot(res, Cex = 1,widthtree = 10, heighttree = 10,lwd = 2)
importance <- VarImp(s.auto, pretty = F)
importance$Relative

########################
##FIGURE 5 Two-dimensional plot of the multiple factor analysis performed for our whole dataset, including multi-scale environmental indices and otolith morphological characters
df11 <- tt[,4:13]
df1 <- decostand(df11, "hellinger")
names(df1)
df2 <- tt[,14:17]
names(df2)
#pop <-tt[,1]
#land.use <-tt[,17]
#df3 <- data.frame(pop,land.use)
library(MASS)
#farms.mca <- mca(df3, abbrev=TRUE)
#names(df3)
tab3 <- data.frame(df1,df2)
dim(tab3)

(num <- c(ncol(df1), ncol(df2)))

# Compute the MFA with multiple plots
# graphics.off()  # Close the previous graphic windows
t3.mfa <- MFA(tab3, group = num, type = c("c","s"), ncp = 2,
              name.group = c("Otolith morphology indices","Environmental factors"))
t3.mfa
#dimdesc(t3.mfa )
palette = palette(c("black","gray50"))
plot(t3.mfa, choix = "var", habillage = "group",axes = 1:2,palette = palette)
plot(t3.mfa, choix = "ind", habillage = "none")
plot(t3.mfa, choix = "ind", habillage = "none", partial = "all")

#plot(t3.mfa, choix="quanti.var")
(rvp <- t3.mfa$group$RV)
rvp[1,2] <- coeffRV(df11,scale(df2))$p.value
rvp[1,2]
##################
MCC(s.auto,standard = F)
MCC(s.auto,weight = F)
Interaction <- MCC(s.auto ,weight = T)
Interaction$interact
s.auto$frame






library(FactoMineR)
data(decathlon)
View(decathlon)
res.pca <- PCA(decathlon,quanti.sup=11:12,quali.sup=13)
res.pca
round(res.pca$eig,2)
barplot(res.pca$eig[,1],main="Eigenvalues",
        names.arg=paste("dim",1:nrow(res.pca$eig)))
round(cbind(res.pca$ind$coord[,1:4],res.pca$ind$cos2[,1:4],
            res.pca$ind$contrib[,1:4]),2)
round(cbind(res.pca$quali.sup$coord[,1:4],res.pca$quali.sup$cos2[,1:4],res.pca$quali.sup$vtest[,1:4]),2)
plot(res.pca,choix="ind",axes=3:4)
plot(res.pca,choix="ind",habillage=13,cex=0.7)

concat.data <- cbind.data.frame(decathlon[,13],res.pca$ind$coord)
ellipse.coord <- coord.ellipse(concat.data,bary=TRUE)
plot.PCA(res.pca,habillage=13,ellipse=ellipse.coord,cex=0.8)
round(cbind(res.pca$var$coord[,1:4],res.pca$var$cos2[,1:4],
           res.pca$var$contrib[,1:4]),2)
round(cbind(res.pca$quanti.sup$coord[,1:4],res.pca$quanti.sup$cos2[,1:4]),2)
plot(res.pca,choix="var",axes=3:4)
dimdesc(res.pca)
dimdesc(res.pca,proba=0.2)
res.pca$call$centre
res.pca$call$ecart.type
round(scale(decathlon[,1:12]),2)
round(cor(decathlon[,1:12]),2)
pairs(decathlon[,c(1,2,6,10)])



#example2
library(FactoMineR)
temperature <- read.table("http://factominer.free.fr/book/temperature.csv", header=TRUE,sep=";",dec=".",row.names=1)
res<-PCA(temperature,ind.sup=24:35,quanti.sup=13:16,quali.sup=17)
plot.PCA(res,choix="ind",habillage=17)
dimdesc(res)
scale(temperature[1:23,1:16])*sqrt(22/23)
cor(temperature[1:23,1:16])
dimdesc(res)



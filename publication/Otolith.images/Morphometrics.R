###FIGURE 3 Reconstruction of the mean otolith shape in Schizothorax nukiangensis, based on Wavelet transformation, differentiated by sampling site

library("shapeR")

#read data and image
shape1 = shapeR("Otolith.images/img_bzlz","BZLZ.39.csv")
shape1 = detect.outline(shape1, threshold = 0.2, write.outline.w.org = TRUE)
shape1 = smoothout(shape1, n = 100)
shape1 = generateShapeCoefficients(shape1)
shape1 = enrich.master.list(shape1)
data1 <- getMeasurements(shape1)
write.csv(data1,file = "Otolith.images/img_bzlz/measure.BZLZ.csv")
plotWaveletShape(shape1, "pop", lwd = 4,lty = 1)
#########################################################
shape2 = shapeR("Otolith.images/img_cwlz","CWLZ.79.csv")
shape2 = detect.outline(shape2, threshold = 0.2, write.outline.w.org = TRUE)
shape2 = smoothout(shape2, n = 100)
shape2 = generateShapeCoefficients(shape2)
shape2 = enrich.master.list(shape2)
data2 <- getMeasurements(shape2)
write.csv(data2,file = "Otolith.images/img_cwlz/measure.CWLZ.csv")
plotWaveletShape(shape2, "pop", lwd = 4,lty = 1,asp = 1)

###########################################
shape3 = shapeR("Otolith.images/img_lkz","LKZ.87.csv")
shape3 = detect.outline(shape3, threshold = 0.2, write.outline.w.org = TRUE)
shape3 = smoothout(shape3, n = 100)
shape3 = generateShapeCoefficients(shape3)
shape3 = enrich.master.list(shape3)
data3 <- getMeasurements(shape3)
write.csv(data3,file = "Otolith.images/img_lkz/measure.LKZ.csv")
plotWaveletShape(shape3, "pop", lwd = 4,lty = 1,asp = 1)

###########################################
shape4 = shapeR("Otolith.images/img_mkz","MKZ.138.csv")
shape4 = detect.outline(shape4, threshold = 0.2, write.outline.w.org = TRUE)
shape4 = smoothout(shape4, n = 100)
shape4 = generateShapeCoefficients(shape4)
shape4 = enrich.master.list(shape4)
data4 <- getMeasurements(shape4)
write.csv(data4,file = "Otolith.images/img_mkz/measure.MKZ.csv")
plotWaveletShape(shape4, "pop", lwd = 4,lty = 1,asp = 1)

###########################################
shape5 = shapeR("Otolith.images/img_spz","SPZ.12.csv")
shape5 = detect.outline(shape5, threshold = 0.2, write.outline.w.org = TRUE)
shape5 = smoothout(shape5, n = 100)
shape5 = generateShapeCoefficients(shape5)
shape5 = enrich.master.list(shape5)
data5 <- getMeasurements(shape5)
write.csv(data5,file = "Otolith.images/img_spz/measure.SPZ.csv")
plotWaveletShape(shape5, "pop", lwd = 4,lty = 1,asp = 1)

###########################################
shape6 = shapeR("Otolith.images/img_ttz","TTZ.72.csv")
shape6 = detect.outline(shape6, threshold = 0.2, write.outline.w.org = TRUE)
shape6 = smoothout(shape6, n = 100)
shape6 = generateShapeCoefficients(shape6)
shape6 = enrich.master.list(shape6)
data6 <- getMeasurements(shape6)
write.csv(data6,file = "Otolith.images/img_ttz/measure.TTZ.csv")
plotWaveletShape(shape6, "pop", lwd = 4,lty = 1,asp = 1)

###########################################
shape7 = shapeR("Otolith.images/img_zyz","ZYZ.66.csv")
shape7 = detect.outline(shape7, threshold = 0.2, write.outline.w.org = TRUE)
shape7 = smoothout(shape7, n = 100)
shape7 = generateShapeCoefficients(shape7)
shape7 = enrich.master.list(shape7)
data7 <- getMeasurements(shape7)
write.csv(data7,file = "Otolith.images/img_zyz/measure.ZYZ.csv")
plotWaveletShape(shape7, "pop", lwd = 4,lty = 1,asp = 1)

########################
#shape8 = shapeR("D:/Documents/R Workplace/Otolith.images/img_all","All.493.csv")
#shape8 = detect.outline(shape8, threshold = 0.2, write.outline.w.org = TRUE)
#shape8 = smoothout(shape8, n = 100)
#shape8 = generateShapeCoefficients(shape8)
#shape8 = enrich.master.list(shape8)
#data8 <-getMeasurements(shape8)
#write.csv(data8,file = "D:/Documents/R Workplace/Otolith.images/img_all/measure.All.csv")
#plotWaveletShape(shape8, "pop", lwd = 4,lty = 1,asp = 1)
########################
a <- read.csv("Otolith.images/img_bzlz/BZLZ.39.csv")
b <- read.csv("Otolith.images/img_cwlz/CWLZ.79.csv")
c <- read.csv("Otolith.images/img_lkz/LKZ.87.csv")
d <- read.csv("Otolith.images/img_mkz/MKZ.138.csv")
e <- read.csv("Otolith.images/img_spz/SPZ.12.csv")
f <- read.csv("Otolith.images/img_ttz/TTZ.72.csv")
g <- read.csv("Otolith.images/img_zyz/ZYZ.66.csv")

newdata <- data.frame(rbind(a,b,c,d,e,f,g))
write.csv(newdata,file = "Otolith.images/img_total/newdata.csv")
shape = shapeR("Otolith.images/img_total","newdata.csv")
shape = detect.outline(shape, threshold = 0.2, write.outline.w.org = TRUE)
shape = smoothout(shape, n = 100)
shape = generateShapeCoefficients(shape)
shape = enrich.master.list(shape)
data <- getMeasurements(shape)
write.csv(data,file = "Otolith.images/img_total/measure.Total.csv")
plotWaveletShape(shape, "pop", show.angle = TRUE,lwd = 2,lty = 1,col=1:7)

#Mean otolith shape based on wavelet reconstruction
#pdf("D:/Documents/R Workplace/length_100_200/otolith.pdf",family="GB1")
##plotWaveletShape(shape1, "pop", show.angle = TRUE, lwd = 2,lty = 1)
#dev.off()
#plotFourierShape(shape, "pop",show.angle = TRUE,lwd=2,lty =2,col=1:7)
#legend(-1, 0.9,c("MKZ", "LKZ", "SPZ","BZLZ", "CWLZ", "LKXZ","ZYZ","TTZ"),col=1:7,lty = 1)
tapply(getMeasurements(shape)$otolith.area, getMasterlist(shape)$pop, mean)
#BZLZ     CWLZ      LKZ      MKZ      SPZ      TTZ      ZYZ
#2.114314 2.501707 1.482099 1.174326 1.808778 1.074255 1.546961
tapply(getMeasurements(shape)$otolith.length, getMasterlist(shape)$pop, mean)
#BZLZ     CWLZ      LKZ      MKZ      SPZ      TTZ      ZYZ
#1.858856 2.063710 1.563033 1.335997 1.711345 1.289928 1.606169
tapply(getMeasurements(shape)$otolith.width, getMasterlist(shape)$pop, mean)
#BZLZ     CWLZ      LKZ      MKZ      SPZ      TTZ      ZYZ
#1.370629 1.563305 1.211755 1.086305 1.314009 1.021240 1.219364
tapply(getMeasurements(shape)$otolith.perimeter, getMasterlist(shape)$pop, mean)
#BZLZ     CWLZ      LKZ      MKZ      SPZ      TTZ      ZYZ
#5.283325 5.948831 4.534124 3.936989 4.927616 3.757546 4.610932
est.list = estimate.outline.reconstruction(shape)
outline.reconstruction.plot(est.list, max.num.harmonics = 15)
shape = stdCoefs(shape, classes = "pop", "Standard_length.mm.", bonferroni = FALSE)
#op <- par(mfrow = c(1, 1))
plotWavelet(shape, level = 5, class.name = "pop", useStdcoef = TRUE)
###############################################
library(gplots)
library(vegan)
pdf("length_100_200/kkk.pdf",family="GB1")
shape <- stdCoefs(shape, classes = "pop", "Standard_length.mm.", bonferroni = FALSE)
cap.res <- capscale(getStdWavelet(shape) ~ getMasterlist(shape)$pop)
anova(cap.res, by = "terms", step = 1000)
eig = eigenvals(cap.res,constrained = T)
eig.ratio = eig/sum(eig)
cluster.plot(scores(cap.res)$sites[,1:2],getMasterlist(shape)$pop,
             xlim = range(scores(cap.res)$sites[,1]),
             ylim = range(scores(cap.res)$sites[,2]),
             xlab = paste("CAP1 (",round(eig.ratio[1]*100,1),"%)",sep = ""),
             ylab = paste("CAP2 (",round(eig.ratio[2]*100,1),"%)",sep = ""),              plotCI = TRUE,
             conf.level = 0.95,las = 1)
abline(v = 0,lty = 3)
abline(h = 0,lty = 3)
dev.off()

####SUPPLEMENTARY

Morphometrics <- read.csv("data/Morphometrics.csv",header = T)
str(Morphometrics)
groups = factor(rep(c("BZLZ","CWLZ","LKZ","MKZ","SPZ","TTZ","ZYZ"),
                    c(39,79,87,138,12,72,66)))



##############relative.otolith.area
Relative.OA <- data.frame(Morphometrics$relative.otolith.area,groups)
colnames(Relative.OA)[1] <- 'relative.otolith.area'
#The mean of each group
aggregate(Relative.OA$relative.otolith.area,by = list(Relative.OA$groups),FUN = mean)
#Standard deviation of each group
aggregate(Relative.OA$relative.otolith.area,by = list(Relative.OA$groups),FUN = sd)
bartlett.test(Relative.OA$relative.otolith.area,Relative.OA$groups)
kruskal.test(Relative.OA$relative.otolith.area~Relative.OA$groups)

m1 <- aov(Relative.OA$relative.otolith.area~Relative.OA$groups)
summary(m1)
TukeyHSD(m1)
re1 <- TukeyHSD(m1)




##############relative.otolith.length
Relative.OL <- data.frame(Morphometrics$relative.otolith.length,groups)
colnames(Relative.OL)[1] <- 'relative.otolith.length'
#The mean of each group
aggregate(Relative.OL$relative.otolith.length,by = list(Relative.OL$groups),FUN = mean)
#Standard deviation of each group
aggregate(Relative.OL$relative.otolith.length,by = list(Relative.OL$groups),FUN = sd)
bartlett.test(Relative.OL$relative.otolith.length,Relative.OL$groups)
kruskal.test(Relative.OL$relative.otolith.length~Relative.OL$groups)
m2 <- aov(Relative.OL$relative.otolith.length~Relative.OL$groups)
summary(m2)
TukeyHSD(m2)
re2 <- TukeyHSD(m2)


##############relative.otolith.width
Relative.OW <- data.frame(Morphometrics$relative.otolith.width,groups)
colnames(Relative.OW)[1] <- 'relative.otolith.width'
#The mean of each group
aggregate(Relative.OW$relative.otolith.width,by = list(Relative.OW$groups),FUN = mean)
#Standard deviation of each group
aggregate(Relative.OW$relative.otolith.width,by = list(Relative.OW$groups),FUN = sd)
bartlett.test(Relative.OW$relative.otolith.width~Relative.OW$groups)
kruskal.test(Relative.OW$relative.otolith.width~Relative.OW$groups)
m3 <- aov(Relative.OW$relative.otolith.width~Relative.OW$groups)
summary(m3)
TukeyHSD(m3)
re3 <- TukeyHSD(m3)

##############relative.otolith.perimeter
Relative.OP <- data.frame(Morphometrics$relative.otolith.perimeter,groups)
colnames(Relative.OP)[1] <- 'relative.otolith.perimeter'
#The mean of each group
aggregate(Relative.OP$relative.otolith.perimeter,by = list(Relative.OP$groups),FUN = mean)
#Standard deviation of each group
aggregate(Relative.OP$relative.otolith.perimeter,by = list(Relative.OP$groups),FUN = sd)
bartlett.test(Relative.OP$relative.otolith.perimeter~Relative.OP$groups)
kruskal.test(Relative.OP$relative.otolith.perimeter~Relative.OP$groups)
m4 <- aov(Relative.OP$relative.otolith.perimeter~Relative.OP$groups)
summary(m4)
TukeyHSD(m4)
re4 <- TukeyHSD(m4)

##############Aspect ratio
Aspect.ratio <- data.frame(Morphometrics$Aspect.ratio,groups)
colnames(Aspect.ratio)[1] <- 'Aspect.ratio'
#The mean of each group
aggregate(Aspect.ratio$Aspect.ratio,by = list(Aspect.ratio$groups),FUN = mean)
#Standard deviation of each group
aggregate(Aspect.ratio$Aspect.ratio,by = list(Aspect.ratio$groups),FUN = sd)
bartlett.test(Aspect.ratio$Aspect.ratio~Aspect.ratio$groups)
kruskal.test(Aspect.ratio$Aspect.ratio~Aspect.ratio$groups)
m5 <- aov(Aspect.ratio$Aspect.ratio~Aspect.ratio$groups)
summary(m5)
TukeyHSD(m5)
re5 <- TukeyHSD(m5)
##############Rectangularity
Rectangularity <- data.frame(Morphometrics$Rectangularity,groups)
colnames(Rectangularity)[1] <- 'Rectangularity'
#The mean of each group
aggregate(Rectangularity$Rectangularity,by = list(Rectangularity$groups),FUN = mean)
#Standard deviation of each group
aggregate(Rectangularity$Rectangularity,by = list(Rectangularity$groups),FUN = sd)
bartlett.test(Rectangularity$Rectangularity~Rectangularity$groups)
kruskal.test(Rectangularity$Rectangularity~Rectangularity$groups)
m6 <- aov(Rectangularity$Rectangularity~Rectangularity$groups)
summary(m6)
TukeyHSD(m6)
re6 <- TukeyHSD(m6)

##############Roundness
Roundness <- data.frame(Morphometrics$Roundness,groups)
colnames(Roundness)[1] <- 'Roundness'
#The mean of each group
aggregate(Roundness$Roundness,by = list(Roundness$groups),FUN = mean)
#Standard deviation of each group
aggregate(Roundness$Roundness,by = list(Roundness$groups),FUN = sd)
bartlett.test(Roundness$Roundness~Roundness$groups)
kruskal.test(Roundness$Roundness~Roundness$groups)
m7 <- aov(Roundness$Roundness~Roundness$groups)
summary(m7)
TukeyHSD(m7)
re7 <- TukeyHSD(m7)

##############Form factor
Form.factor <- data.frame(Morphometrics$Form.factor,groups)
colnames(Form.factor)[1] <- 'Form.factor'
#The mean of each group
aggregate(Form.factor$Form.factor,by = list(Form.factor$groups),FUN = mean)
#Standard deviation of each group
aggregate(Form.factor$Form.factor,by = list(Form.factor$groups),FUN = sd)
bartlett.test(Form.factor$Form.factor~Form.factor$groups)
kruskal.test(Form.factor$Form.factor~Form.factor$groups)
m8 <- aov(Form.factor$Form.factor~Form.factor$groups)
summary(m8)
TukeyHSD(m8)
re8 <- TukeyHSD(m8)

##############Ellipticity
Ellipticity <- data.frame(Morphometrics$Ellipticity,groups)
colnames(Ellipticity)[1] <- 'Ellipticity'
#The mean of each group
aggregate(Ellipticity$Ellipticity,by = list(Ellipticity$groups),FUN = mean)
#Standard deviation of each group
aggregate(Ellipticity$Ellipticity,by = list(Ellipticity$groups),FUN = sd)
bartlett.test(Ellipticity$Ellipticity~Ellipticity$groups)
kruskal.test(Ellipticity$Ellipticity~Ellipticity$groups)
m9 <- aov(Ellipticity$Ellipticity~Ellipticity$groups)
summary(m9)
TukeyHSD(m9)
re9 <- TukeyHSD(m9)

##############Circularity
Circularity <- data.frame(Morphometrics$Circularity,groups)
colnames(Circularity)[1] <- 'Circularity'
#The mean of each group
aggregate(Circularity$Circularity,by = list(Circularity$groups),FUN = mean)
#Standard deviation of each group
aggregate(Circularity$Circularity,by = list(Circularity$groups),FUN = sd)
bartlett.test(Circularity$Circularity~Circularity$groups)
kruskal.test(Circularity$Circularity~Circularity$groups)
m10 <- aov(Circularity$Circularity~Circularity$groups)
summary(m10)
TukeyHSD(m10)
re10 <- TukeyHSD(m10)
##############################
par(mfrow = c(2,4),mar = c(4, 4, 0.5, 0.5))
plot(re2,las = 1)
plot(re3,las = 1)
plot(re4,las = 1)
plot(re5,las = 1)
plot(re6,las = 1)
plot(re7,las = 1)
plot(re9,las = 1)
dev.off()







###FIGURE S1 TukeyHSD-corrected P-values shown as a heatmap. ANOVAs wereused to calculate differences among means for each of the otolith morphological parameters (heatmap columns)

library(readxl)
library(RColorBrewer)
data <- read_excel("data/data.xlsx")
x11 <- data[,2:11]
y11 <- data.matrix(x11)
rownames(y11) = data$Group
heatmap(y11,margins = c(8,1))
library(pheatmap)
windowsFonts(myFont1 = windowsFont("Times New Roman"))
pheatmap(y11,cellwidth = 60, cellheight = 20,fontsize = 12, fontsize_row = 12,display_numbers = TRUE,
         color = colorRampPalette(brewer.pal(n = 8, name = "RdBu"))(100),
         number_format = "%.7f",FontFamily = "myfont1")

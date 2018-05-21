data1 <- read.csv("C:/Users/Administrator/Desktop/data1.csv", header = T)
data2 <- read.csv("C:/Users/Administrator/Desktop/data2.csv", header = T)
fish <- dplyr::inner_join(data1,data2,by = "ID")
write.csv(fish, file = "C:/Users/Administrator/Desktop/fish.csv")
#data <- merge(data1,data2,by="ID",all=T)
#data <- dplyr::bind_rows(data1,data2)
#data <- dplyr::combine(data1,data2)

data(Stickleback)
Stickleback = as.data.frame(Stickleback)
GutContents_SiteA <- import.RInSp(Stickleback, row.names = 1,
                                  info.cols = c(2:13), subset.rows = c("Site", "A"))
PSi <- PSicalc(GutContents_SiteA, exclude = FALSE, replicates = 99)
sumMC.RInSp(PSi)



library("RInSp")
fish <- read.csv("C:/Users/Administrator/Desktop/fish.csv", header = T)
##Variance partition and WIC/TNW
fish <- import.RInSp(fish, col.header = TRUE,
                     row.names = 1, info.cols = 2,
                     data.type = "double")
decomp = Hier2L(fish, factor = 1)
#Calculate index E and Cws


fish_SiteBZLZ <- import.RInSp(fish, row.names = 1,
      info.cols = c(2:12), subset.rows = c("Site", "BZLZ"))
fish_SiteCWLZ <- import.RInSp(fish, row.names = 1,
      info.cols = c(2:12), subset.rows = c("Site", "CWLZ"))
fish_SiteLKZ <- import.RInSp(fish, row.names = 1,
      info.cols = c(2:12), subset.rows = c("Site", "LKZ"))
fish_SiteMKZ <- import.RInSp(fish, row.names = 1,
      info.cols = c(2:12), subset.rows = c("Site", "MKZ"))
fish_SiteSPZ <- import.RInSp(fish, row.names = 1,
      info.cols = c(2:12), subset.rows = c("Site", "SPZ"))
fish_SiteTTZ <- import.RInSp(fish, row.names = 1,
       info.cols = c(2:12), subset.rows = c("Site", "TTZ"))
fish_SiteZYZ <- import.RInSp(fish, row.names = 1,
       info.cols = c(2:12), subset.rows = c("Site", "ZYZ"))

par(mfrow = c(2,4))

PSi_BZLZ <- PSicalc(fish_SiteBZLZ, exclude = FALSE,
                    replicates = 999)
sumMC.RInSp(PSi_BZLZ)
PSi_CWLZ <- PSicalc(fish_SiteCWLZ, exclude = FALSE,
                    replicates = 999)
sumMC.RInSp(PSi_CWLZ)
PSi_LKZ <- PSicalc(fish_SiteLKZ, exclude = FALSE,
                   replicates = 999)
sumMC.RInSp(PSi_LKZ)
PSi_MKZ <- PSicalc(fish_SiteMKZ, exclude = FALSE,
                   replicates = 999)
sumMC.RInSp(PSi_MKZ)
PSi_SPZ <- PSicalc(fish_SiteSPZ, exclude = FALSE,
                   replicates = 999)
sumMC.RInSp(PSi_SPZ)
PSi_TTZ <- PSicalc(fish_SiteTTZ, exclude = FALSE,
                   replicates = 999)
sumMC.RInSp(PSi_TTZ)
PSi_ZYZ <- PSicalc(fish_SiteZYZ, exclude = FALSE,
                   replicates = 999)
sumMC.RInSp(PSi_ZYZ)
dev.off()

#Monte Carlo resampling of WIC/TNW for the discrete case
par(mfrow = c(2,4))


WT_BZLZ = WTdMC(fish_SiteBZLZ, replicates = 10000)
sumMC.RInSp(WT_BZLZ)

WT_CWLZ = WTdMC(fish_SiteCWLZ, replicates = 10000)
sumMC.RInSp(WT_CWLZ)

WT_LKZ = WTdMC(fish_SiteLKZ, replicates = 10000)
sumMC.RInSp(WT_LKZ)

WT_MKZ = WTdMC(fish_SiteMKZ, replicates = 10000)
sumMC.RInSp(WT_MKZ)

WT_SPZ = WTdMC(fish_SiteSPZ, replicates = 10000)
sumMC.RInSp(WT_SPZ)

WT_TTZ = WTdMC(fish_SiteTTZ, replicates = 10000)
sumMC.RInSp(WT_TTZ)

WT_ZYZ = WTdMC(fish_SiteZYZ, replicates = 10000)
sumMC.RInSp(WT_ZYZ)
dev.off()

#######################
BZLZ <- as.data.frame(WT_BZLZ$montecarlo)
CWLZ <- as.data.frame(WT_CWLZ$montecarlo)
LKZ <- as.data.frame(WT_LKZ$montecarlo)
MKZ <- as.data.frame(WT_MKZ$montecarlo)
SPZ <- as.data.frame(WT_SPZ$montecarlo)
TTZ <- as.data.frame(WT_TTZ$montecarlo)
ZYZ <- as.data.frame(WT_ZYZ$montecarlo)

plot(density(BZLZ$WIC/BZLZ$TNW),
     xlab = "WIC/TNW",
     ylab = "density",
     main = "",
     xlim = c(0.99,1.01),
     lwd = 2,
     col = "black")
lines(density(CWLZ$WIC/CWLZ$TNW), col = "blue", lwd = 2)
lines(density(LKZ$WIC/LKZ$TNW), col = "red", lwd = 2)
lines(density(MKZ$WIC/MKZ$TNW), col = "yellow", lwd = 2)
lines(density(SPZ$WIC/SPZ$TNW), col = "brown", lwd = 2)
lines(density(TTZ$WIC/TTZ$TNW), col = "chartreuse4", lwd = 2)
lines(density(ZYZ$WIC/ZYZ$TNW), col = "coral4", lwd = 2)

text.legend = c("BZLZ","CWLZ","LKZ","MKZ",
                "SPZ","TTZ","ZYZ")
col2 <- c("black","blue","red","yellow","brown",
          "chartreuse4","coral4")

legend("topright",legend = text.legend,
       lty = c(1,1,1,1,1,1,1),
       lwd = c(2,2,2,2,2,2,2),
             col = col2,bty = "n")



Eresult <- Eindex(fish_SiteZYZ, index = "saramaki",
                  jackknife = TRUE)


similarity = overlap(fish_SiteBZLZ)




library("WDI")
dat = WDI(indicator='NY.GDP.PCAP.KD', country= "US" , start=1990, end=2017)
write.csv(dat, file = "C:/Users/Administrator/Desktop/dat.csv")
WDIsearch(string ='gdp', field ='name', cache=NULL)
WDIsearch('gdp_country')[350:400,]
#"SE.XPD.TOTL.GD.ZS"
#"Government expenditure on education, total (% of GDP)"

#"SE.TER.CUAT.BA.ZS"
#"Educational attainment, at least Bachelor's or equivalent, population 25+, total (%) (cumulative)"
ind1 = c("gdp_Gov_education" = "SE.XPD.TOTL.GD.ZS",
        "population_Bachelor" = "SE.TER.CUAT.BA.ZS")
#----------------------------------
library("WDI")

Country <- read.csv(file = "C:/Users/Administrator/Desktop/Country.csv")
ind = c("gdp_per_capita" = "NY.GDP.PCAP.KD",
        "population" = "SP.POP.TOTL")
iso <- Country[,2]
gdp_pop = WDI(indicator = ind, country = iso,
          start = 1990, end = 2016)
write.csv(gdp_pop, file = "C:/Users/Administrator/Desktop/gdp_pop.csv")


library(scatterplot3d)
Country <- read.csv(file = "data/Country.csv")
layout(cbind(1:2, 1:2), heights = c(7, 1))
prc <- hsv((prc <- 0.7 * Country$TLCS / diff(range(Country$TLCS))) - min(prc) + 0.3)
s3d <- scatterplot3d(y = Country$Population,
                     z = Country$Records,
                     x = Country$GDP_per_capita,
                     pch = 16,
                     zlim = c(0, 1200,100),
                     angle = -158,
                     highlight.3d = F,
                     type = "h", main = "",
                     xlab = "国家人均GDP(10^4,US dollars)",
                     ylab = "人口数量(10^8)",
                     zlab = "文章数量")

fit <- lm(Country$Records ~ Country$Population+Country$GDP_per_capita)
s3d$plane3d(fit,col = "gray20")
summary(fit)
fit

s3d$points3d(y = Country$Population,
             z = Country$TLCS,
             x = Country$GDP_per_capita,
             col = "blue",
             pch = 21)




s3d$points(y = Country$Population,
           z = Country$TLCS,
           x = Country$GDP_per_capita, pch = 21,bg = prc)

par(mar=c(5, 3, 0, 3))
plot(seq(min(Country$TLCS), max(Country$TLCS), length = 100), rep(0, 100),
     axes = FALSE, xlab = "color code of variable \"Precision\"",
     ylab = "", col = hsv(seq(0.3, 1, length = 100)))
axis(1, at = 4:7, labels = expression(10^4, 10^5, 10^6, 10^7))
dev.off()


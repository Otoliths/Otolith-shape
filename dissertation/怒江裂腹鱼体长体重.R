library("ggplot2")
library("ggExtra")
library("cowplot")
library("ggimage")
library("yyplot")
library("ggridges")
data <- read.csv("./data/Schizothorax nukiangensis.csv")
head(data)
str(data)
#**************************
plot(data[,3:5])
a <- ggplot(data, aes(x = SL, y = Abbreviation)) +
  geom_density_ridges()
b <- ggplot(data, aes(x = BW, y = Abbreviation)) +
  stat_density_ridges()
plot_grid(a, b)
#**********************************
fit1 <- lm(log(data$BW)~log(data$SL),
           data = subset(data,Abbreviation == "MKZ"))
weight <- predict(fit1, se.fit = T,
                  interval = "confidence",level = 0.95)
weight
#length <- log(data$SL[data$Abbreviation == "MKZ"])
plot(log(data$SL),weight$fit[,1],data = subset(data,Abbreviation == "MKZ"),
     xlab = "Stardard Length",ylab = "Body weigth",
     pch = 16,col = 'red')
lines(log(data$SL),weight$fit[,2],
      data = subset(data,Abbreviation == "MKZ"),
      col = 'blue',lty = 2)
lines(log(data$SL),weight$fit[,3],
      data = subset(data,Abbreviation == "MKZ"),
      col = 'blue',lty = 2)
#############################

ggplot(data, aes(SL,BW, colour = Abbreviation)) +
  geom_point() +
  geom_smooth(se = T)
#-----------------------------------
#x <- data$SL
#test <- function(x) {0.00002 * x ^ 2.97889}
#pb <- ggplot(data = data, aes(SL,BW)) +
#      geom_point() +
#      stat_function(fun = test, colour = "red") +
#      geom_text(aes(x = 6.5, y = 2, label = "y=2x+3"), #parse = TRUE) +
#      theme_bw()

#pb <- ggplot(data = data, aes(SL,BW)) +
##      geom_smooth(method = 'loess') +
#      annotate("text", x = 150, y = 1500, parse = F,
#      label = "paste(italic(R) ^ 2, \" = .75\")",size = 5) +
#      theme_bw()
pa <- ggplot(data, aes(SL,BW)) +
      geom_point(size = 2,shape = 21,colour = "black", fill = "#BEBEBE") +
      geom_smooth(span = 0.8, method = 'loess',
                  color = "black",size = 0.5) +
      facet_wrap(~Abbreviation,scales = "free_y") +
      theme_bw() +
      labs(x = "体长(mm)", y = "体重(g)")
#一次性设置字体
g1 <- set_font(pa, family = "Times", size = 2)
#d <- data.frame( x = seq(50,350,100),
#                 y = seq(250,1000,100),
#                 image = "54bf13e7-b15e-4175-862f-f6028c1ee9b3",
#                 size = 10,replace = TRUE
#               )
#d$size = seq(.05, .15, length.out = 8)
pb <- ggplot(data = data, aes(SL,BW)) +
      geom_point(size = 3,shape = 21,
                 colour = "black", fill = "#BEBEBE") +
      geom_smooth(method = 'loess',color = "black",
                  size = 1) +
      geom_phylopic(aes(x = 50,y = 250),
           image = "54bf13e7-b15e-4175-862f-f6028c1ee9b3",              size = 0.1, height = 200,
           alpha = .6,inherit.aes = F) +
      theme_bw() +
      labs(x = "体长(mm)", y = "体重(g)")


pbb <- ggMarginal(pb, x = "SL", y = "BW",
                type = "histogram",
                xparams = list(breaks = seq(10,500,10)),
                yparams = list(breaks = seq(100,600,10))
                )
g2 <- set_font(pb, family = "Times", size = 3)
x11()
plot_grid(g1, g2, align = "h", axis = "b",
          nrow = 1, rel_widths = c(1,1),
          labels = c('(a)', '(b)'))
ggsave("C:/Users/Administrator/Desktop/ppp.tiff",
       width = 16, height = 20, units = "cm")
###------------------------------------------
p1 <- ggplot(data = subset(data,Abbreviation == "MKZ"),
      aes(SL,BW)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      ggtitle("芒宽镇") +
      theme_bw()

p2 <- ggplot(data = subset(data,Abbreviation == "LKZ"), aes(SL,BW)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  theme_bw()

p3 <- ggplot(data = subset(data,Abbreviation == "SPZ"), aes(SL,BW)) +
  geom_point() +
  geom_smooth(method = 'loess',se = F) +
  theme_bw()

p4 <- ggplot(data = subset(data,Abbreviation == "BZLZ"),
             aes(SL,BW)) +
            geom_point() +
            geom_smooth(method = 'loess') +
            theme_bw()

p5 <- ggplot(data = subset(data,Abbreviation == "CWLZ"),
             aes(SL,BW)) +
             geom_point() +
             geom_smooth(method = 'loess') +
             theme_bw()

p6 <- ggplot(data = subset(data,Abbreviation == "ZYZ"),
             aes(SL,BW)) +
             geom_point() +
             geom_smooth(method = 'loess') +
             theme_bw()

p7 <- ggplot(data = subset(data,Abbreviation == "TTZ"),
             aes(SL,BW)) +
             geom_point() +
             geom_smooth(method = 'loess') +
             theme_bw()

p8 <- ggplot(data = data, aes(SL,BW)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      theme_bw()



cowplot::plot_grid(p1, p2, p3, p4, p5, p6,p7,p8,
                   ncol = 3, label_x = "SL(mm)",label_y = "BW(g)")
ggsave("C:/Users/Administrator/Desktop/p1.tiff",
       width = 18, height = 25, units = "cm")
#require(Cairo)
#CairoTIFF("C:/Users/Administrator/Desktop/p1.pdf",
#         3.15, 3.15, encoding = "utf-8") #单位为英寸
#dev.off() #关闭图像设备，同时储存图片
#------------------------------------------


p1 <- ggplot(data, aes(SL,BW)) +
      geom_point() +
      geom_smooth(span = 0.8)
ggMarginal(p1, x = "SL", y = "BW",
           type = "histogram",
           xparams = list(breaks = seq(10,500,10)),
           yparams = list(breaks = seq(100,600,10))
           )



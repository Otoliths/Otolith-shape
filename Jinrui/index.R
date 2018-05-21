### Web of science
### TS= (soil organic carbon* AND (soil erosion* OR Soil denudation*))
### 索引=SCI-EXPANDED, SSCI, A&HCI, CPCI-S, CPCI-SSH, BKCI-S, BKCI-SSH, ESCI, CCR-EXPANDED, IC
### timeline=1900-2017
### records=3162
### load bibliometrix package

library("bibliometrix")

#D2 <- readFiles('./data/Scopus/1997-2003(293).bib')
#M2 <- convert2df(D2, dbsource = "scopus", format = "bibtex")
D <- readFiles('./WOS/savedrecs_1-500.txt',
                './WOS/savedrecs_501-1000.txt',
                './WOS/savedrecs_1001-1500.txt',
                './WOS/savedrecs_1501-2000.txt',
                './WOS/savedrecs_2001-2500.txt',
                './WOS/savedrecs_2501-3000.txt',
                './WOS/savedrecs_3001-3162.txt')

M <- convert2df(D, dbsource = "isi")
class(M)

results <- biblioAnalysis(M, sep = ";")
CountryCollaboration <- results[["CountryCollaboration"]]
write.csv(CountryCollaboration,file = "./data/CountryCollaboration1.csv")

S <- summary(object = results, k = 200, pause = FALSE)
AnnualProduction <- S[["AnnualProduction"]]
write.csv(AnnualProduction,file = "./data/AnnualProduction1.csv")
Countries <- data.frame(results$Countries)
library(plyr)
Countries <- rename(Countries,c(Tab = "Country"))
#write.csv(Countries,file = "./data/Countries.csv")
TCperCountries <- S[["TCperCountries"]]
Country <- cbind(Countries,TCperCountries)
Country <- Country[,-3]
write.csv(Country,file = "./data/Countries.csv")


MostProdCountries <- S[["MostProdCountries"]]
write.csv(MostProdCountries,file = "./data/MostProdCountries1.csv")
TCperCountries <- S[["TCperCountries"]]
write.csv(TCperCountries,file = "./data/TCperCountries1.csv")


data <- M
M1 <- subset(data,PY >= 1990 & PY <= 1996)
M2 <- subset(data,PY >= 1997 & PY <= 2003)
M3 <- subset(data,PY >= 2004 & PY <= 2010)
M4 <- subset(data,PY >= 2011 & PY <= 2017)

#合并data.frame
M <- dplyr::bind_rows(M1,M2,M3,M4)
rownames(M) <- M$SR

NetMatrix1111 <- biblioNetwork(M1, analysis = "co-occurrences",
                            network = "author_keywords", sep = ";")
S1111 <- normalizeSimilarity(NetMatrix1111, type = "association")
net1111 <- networkPlot(NetMatrix1111, normalize = "association",n = dim(NetMatrix1111)[1],
                    Title = "co-occurrence network",type="fruchterman",
                    labelsize = 0.7, halo = FALSE, cluster = "walktrap",remove.isolates=FALSE,
                    remove.multiple=FALSE, noloops=TRUE, weighted=TRUE)
res1111 <- thematicMap(net1111, NetMatrix1111, S1111)
plot(res1111$map)
#=============================
NetMatrix2222 <- biblioNetwork(M2, analysis = "co-occurrences",
                               network = "author_keywords", sep = ";")
S2222 <- normalizeSimilarity(NetMatrix2222, type = "association")
net2222 <- networkPlot(NetMatrix2222, normalize = "association",n = 500,
                       Title = "co-occurrence network",type="fruchterman",
                       labelsize = 0.7, halo = FALSE, cluster = "walktrap",remove.isolates = FALSE,
                       remove.multiple = FALSE, noloops = TRUE, weighted = TRUE)
res2222 <- thematicMap(net2222, NetMatrix2222, S2222)
plot(res2222$map)
#=====================================================
NetMatrix3333 <- biblioNetwork(M3, analysis = "co-occurrences",
                               network = "author_keywords", sep = ";")
S3333 <- normalizeSimilarity(NetMatrix3333, type = "association")
net3333 <- networkPlot(NetMatrix3333, normalize = "association",n = 700,
                       Title = "co-occurrence network",type="fruchterman",
                       labelsize = 0.7, halo = FALSE, cluster = "walktrap",remove.isolates = FALSE,
                       remove.multiple = FALSE, noloops = TRUE, weighted = TRUE)
res3333 <- thematicMap(net3333, NetMatrix3333, S3333)
plot(res3333$map)
#=====================================================
NetMatrix4444 <- biblioNetwork(M4, analysis = "co-occurrences",
                               network = "author_keywords", sep = ";")
S4444 <- normalizeSimilarity(NetMatrix4444, type = "association")
net4444 <- networkPlot(NetMatrix4444, normalize = "association",n = 1000,
                       Title = "co-occurrence network",type="fruchterman",
                       labelsize = 0.7, halo = FALSE, cluster = "walktrap",remove.isolates = FALSE,
                       remove.multiple = FALSE, noloops = TRUE, weighted = TRUE)
res4444 <- thematicMap(net4444, NetMatrix4444, S4444)
plot(res4444$map)









#M4 <- plyr::rbind.fill(M1,M2)
#str(M3)
write.csv(M,file = "./data/M.csv")
M <- read.csv(file = "./data/M.csv",header = T)
# Create a historical citation network
histResults <- histNetwork(M, n = 10, sep = ".")
nethistory <- histPlot(histResults, size = T,
                       label = TRUE, arrowsize = 0.5)



#write.csv(M2,file = "./data/M2.csv")
#write.csv(M3[,1:21],file = "./data/1997.csv")
results <- biblioAnalysis(M, sep = ";")
CountryCollaboration <- results[["CountryCollaboration"]]
write.csv(CountryCollaboration,file = "./data/CountryCollaboration1.csv")
########Lotka’s Law coefficient estimation
L <- lotka(results)
# Observed distribution
Observed = L$AuthorProd[,3]
# Theoretical distribution with Beta = 2
Theoretical = 10^(log10(L$C) - 2*log10(L$AuthorProd[,1]))
plot(L$AuthorProd[,1],Theoretical,type = "l",
     col = "red",ylim = c(0, 1), xlab = "N. of Articles",
     ylab = "Freq. of Authors",main = "")
lines(L$AuthorProd[,1],Observed,col = "blue")
legend(x = "topright",c("Theoretical (Beta=2)","Observed"),
       col = c("red","blue"),lty = c(1,1,1),
       cex = 0.6,bty = "n")
lotkalaw <- data.frame(L$AuthorProd[,1],Theoretical,Observed)
write.csv(lotkalaw,file = "./data/lotkalaw.csv" )
lotkalaw <- read.csv(file = "./data/lotkalaw.csv" ,header = T)
library(ggplot2)
ggplot(lotkalaw) +
  geom_line(aes(Articles,Freq,color = Group),size = 1) +
  theme_bw() +
  labs(x = "N. of Articles", y = "Freq. of Authors") +
  theme(legend.position = c(.92,.9)) +
  theme(legend.key = element_blank()) +
  theme(legend.background = element_blank())




#str(results)

S <- summary(object = results, k = 50, pause = FALSE)
AnnualProduction <- S[["AnnualProduction"]]
write.csv(AnnualProduction,file = "./data/AnnualProduction.csv")
MostProdCountries <- S[["MostProdCountries"]]
write.csv(MostProdCountries,file = "./data/MostProdCountries.csv")
TCperCountries <- S[["TCperCountries"]]
write.csv(TCperCountries,file = "./data/TCperCountries.csv")

plot(x = results, k = 10, pause = FALSE)

### R-list-to-data-frame列表转数据框
library(plyr)
df <- ldply(results, data.frame)
x1 <- as.data.frame(results[["CountryCollaboration"]])
colnames(x1) <- c("Country","Frep")
x2 <- as.data.frame(results[["TotalCitation"]])
y1 <- cbind(x1,x2[,2:3])
#-------------------------------------------
### Country Scientific Collaboration
### Create a country collaboration network
M1 <- metaTagExtraction(M1, Field = "AU_CO", sep = ";")
NetMatrix1 <- biblioNetwork(M1, analysis = "collaboration",
                           network = "countries", sep = ";")
M2 <- metaTagExtraction(M2, Field = "AU_CO", sep = ";")
NetMatrix2 <- biblioNetwork(M2, analysis = "collaboration",
                            network = "countries", sep = ";")
M3 <- metaTagExtraction(M3, Field = "AU_CO", sep = ";")
NetMatrix3 <- biblioNetwork(M3, analysis = "collaboration",
                            network = "countries", sep = ";")
M4 <- metaTagExtraction(M4, Field = "AU_CO", sep = ";")
NetMatrix4 <- biblioNetwork(M4, analysis = "collaboration",
                            network = "countries", sep = ";")
# Plot the network
X11()
par(mfrow = c(2,2),omi = c(0.1,0.1,0.1,0.1))
net1 <- networkPlot(NetMatrix1, normalize = "equivalence",
                    n = dim(NetMatrix1)[1],
                   Title = "1990-1996",
                   type = "sphere", size = TRUE,
                   remove.multiple = FALSE,labelsize = 0.5)
#m <- network(NetMatrix1)
#p <- (NetMatrix1@x/6)
#plot.network(m,label = network.vertex.names(m),
#             usearrows = FALSE,
#             vertex.cex = p)



net2 <- networkPlot(NetMatrix2, n = dim(NetMatrix2)[1],
                    Title = "1997-2003",
                    type = "sphere", size = TRUE,
                    remove.multiple = FALSE,
                    labelsize = 0.5)
net3 <- networkPlot(NetMatrix3, n = dim(NetMatrix3)[1],
                    Title = "2004-2010",
                    type = "sphere", size = TRUE,
                    remove.multiple = FALSE,labelsize = 0.5)
#net4 <- networkPlot(NetMatrix4, n = dim(NetMatrix4)[1],
#                    Title = "Country Collaboration",
#                    type = "vosviewer",
#                    vos.path = "VOSviewer",
#                    size = TRUE,
#                    remove.multiple = FALSE,labelsize = 0.8)
net4 <- networkPlot(NetMatrix4, n = dim(NetMatrix4)[1],
                    Title = "2011-2017",
                    type = "sphere",
                    size = TRUE,
                    #cluster = "none",
                    remove.multiple = T,labelsize = 0.5)
#-------------------------------------------
### University Scientific Collaboration
### Create a University collaboration network
M11 <- metaTagExtraction(M1, Field = "AU_UN", sep = ";")
NetMatrix11 <- biblioNetwork(M11, analysis = "collaboration",
                            network = "universities", sep = ";")
M22 <- metaTagExtraction(M2, Field = "AU_UN", sep = ";")
NetMatrix22 <- biblioNetwork(M22, analysis = "collaboration",
                            network = "universities", sep = ";")
M33 <- metaTagExtraction(M3, Field = "AU_UN", sep = ";")
NetMatrix33 <- biblioNetwork(M33, analysis = "collaboration",
                            network = "universities", sep = ";")
M44 <- metaTagExtraction(M4, Field = "AU_CO", sep = ";")
NetMatrix44 <- biblioNetwork(M44, analysis = "collaboration",
                            network = "universities", sep = ";")
# Plot the network
X11()
par(mfrow = c(2,2),omi = c(0.1,0.1,0.1,0.1))
net11 <- networkPlot(NetMatrix11, normalize = "equivalence",
                    n = dim(NetMatrix11)[1]/10,
                    Title = "1990-1996",
                    type = "sphere", size = TRUE,
                    remove.multiple = FALSE,labelsize = 0.5)
#m <- network(NetMatrix1)
#p <- (NetMatrix1@x/6)
#plot.network(m,label = network.vertex.names(m),
#             usearrows = FALSE,
#             vertex.cex = p)



net22 <- networkPlot(NetMatrix22, n = dim(NetMatrix22)[1]/10,
                    Title = "1997-2003",
                    type = "sphere", size = TRUE,
                    remove.multiple = FALSE,
                    labelsize = 0.5)
net33 <- networkPlot(NetMatrix33, n = dim(NetMatrix33)[1]/10,
                    Title = "2004-2010",
                    type = "sphere", size = TRUE,
                    remove.multiple = FALSE,labelsize = 0.5)
#net4 <- networkPlot(NetMatrix4, n = dim(NetMatrix4)[1],
#                    Title = "Country Collaboration",
#                    type = "vosviewer",
#                    vos.path = "VOSviewer",
#                    size = TRUE,
#                    remove.multiple = FALSE,labelsize = 0.8)
net44 <- networkPlot(NetMatrix44, n = dim(NetMatrix44)[1]/30,
                    Title = "2011-2017",
                    type = "sphere",
                    size = TRUE,
                    #cluster = "none",
                    remove.multiple = T,labelsize = 0.5)
#-------------------------------------------
### Authors' Coupling

NetMatrix111 <- biblioNetwork(M1, analysis = "coupling",
                              network = "authors", sep = ";")
# plot authors' similarity (first 20 authors), using salton similarity index
net111  <- networkPlot(NetMatrix111, normalize = "salton",
                       weighted = T, n = 20,
                       Title = "1990-1997",
                       type = "sphere", size = FALSE,
                       remove.multiple = TRUE)

NetMatrix222 <- biblioNetwork(M2, analysis = "coupling",
                              network = "authors", sep = ";")
# plot authors' similarity (first 20 authors), using salton similarity index
net222  <- networkPlot(NetMatrix222, normalize = "salton",
                       weighted = T, n = 20,
                       Title = "1998-2003",
                       type = "sphere", size = FALSE,
                       remove.multiple = TRUE)
NetMatrix333 <- biblioNetwork(M3, analysis = "coupling",
                              network = "authors", sep = ";")
# plot authors' similarity (first 20 authors), using salton similarity index
net333  <- networkPlot(NetMatrix333, normalize = "equivalence",
                       weighted = T, n = 20,
                       Title = "2003-2010",
                       type = "sphere", size = FALSE,
                       remove.multiple = TRUE)

NetMatrix444 <- biblioNetwork(M4, analysis = "coupling",
                              network = "authors", sep = ";")
# plot authors' similarity (first 20 authors), using salton similarity index
net444  <- networkPlot(NetMatrix444, normalize = "salton",
                       weighted = T, n = 20,
                       Title = "2011-2017",
                       type = "sphere", size = FALSE,
                       remove.multiple = TRUE)


#-------------------------------------------
### Conceptual Structure using keywords
CS1 <- conceptualStructure(M1,field = "TI",
                          minDegree = 4,
                          k.max = 5,
                          stemming = T,
                          labelsize = 10)
CS2 <- conceptualStructure(M2,field = "TI",
                           minDegree = 4,
                           k.max = 5,
                           stemming = T,
                           labelsize = 10)
CS3 <- conceptualStructure(M3,field = "ID",
                           minDegree = 4,
                           k.max = 5,
                           stemming = FALSE,
                           labelsize = 10)
CS4 <- conceptualStructure(M4,field = "ID",
                           minDegree = 4,
                           k.max = 5,
                           stemming = FALSE,
                           labelsize = 10)
####Yearly occurrences of top keywords/terms
library(reshape2)
library(ggplot2)
topKW1 <- KeywordGrowth(M, Tag = "ID", sep = ";",
                       top = 10, cdf = TRUE)

DF1 <- melt(topKW1, id = 'Year')
ggplot(DF1,aes(Year,value, group = variable,
               color = variable)) +
               geom_line(linetype = 2)





results <- biblioAnalysis(M3, sep = ";")
S <- summary(object = results, k = 10, pause = FALSE)
par(mfrow = c(2,2))
plot(x = results, k = 10, pause = F)
dev.off()

NetMatrix <- biblioNetwork(M3, analysis = "coupling", network = "authors", sep = ";")
net <- networkPlot(NetMatrix, normalize = "salton", weighted = T, n = 20, Title = "Authors' Coupling", type = "kamada", size = FALSE,remove.multiple = F)

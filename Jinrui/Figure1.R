#-----------------------------------------------
library(ggplot2)
library(ggimage)
library(cowplot)
library(yyplot)
library(RColorBrewer)
MostProdCountries <- read.csv(file = "./data/MostProdCountries.csv", header = TRUE)
str(MostProdCountries)
print(levels(MostProdCountries$Country))
MostProdCountries$Country <- factor(MostProdCountries$Country,
                                    levels = c("United States", "China", "Spain",                    "Germany",  "France", "Canada", "England", "Australia", "Brazil", "India"))

Country <- ggplot(MostProdCountries, aes(Country,Documents)) +
  geom_col(aes(fill = Collaboration), width = .5) +
  labs(y = "N. of Documents", x = "Countries") +
  theme_bw() +
  coord_flip() +
  geom_flag(y = -50, aes(image = code)) +
  expand_limits(y = -50) +
  theme(legend.position = c(.8,.9)) +
  theme(legend.key = element_blank()) +
  theme(legend.background = element_blank())
#+geom_text(aes(label = Documents,  hjust = 1, color = collaboration))
#+geom_flag(y = -50, aes(image = code)) + expand_limits(y = -50)
#+scale_fill_manual(values = c("SCP" = "black", "MCP" = "grey"))




AnnualProduction <- read.csv(file = "./data/AnnualProduction1.csv", header = TRUE)
AnnualProduction <- ggplot(AnnualProduction,
                           aes(factor(Year),Articles)) +
  labs(x = "Publication Year", y = "N. of Documents ") +
  geom_bar(stat = "identity",width = .7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,hjust = 0.5,
                                   vjust = 0.5)) +
  geom_text(aes(label = Articles, vjust = -0.5),size = 3)

windows()

plot_grid(AnnualProduction, Country, align = "h", axis = "b",
          nrow = 1, rel_widths = c(1,1.2),
          labels = c('(a)', '(b)'))





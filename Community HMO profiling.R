##############################################################################
## Community HMO profiling data analysis for Renwick et al. 2024 
##############################################################################

### clear workspace / variables
graphics.off()
rm(list=ls())

setwd("D:/")

## Load libraries
library(tidyr)
library(ggplot2)
library(eoffice)

## Load HMO profiling data from community samples treated with HMOs
hmo <- read.csv("Renwick et al. 2024 Community HMO profiling.csv")

## Format dataframe
row.names(hmo) <- hmo$Sample
hmo2 <- hmo

hmo2$Community <- hmo$Community[match(row.names(hmo2), hmo$Sample)]
hmo2$Treatment <- hmo$Treatment[match(row.names(hmo2), hmo$Sample)]
hmo2$Time <- hmo$Time[match(row.names(hmo2), hmo$Sample)]

hmo2$Sample <- NULL

hmo3 <- gather(hmo2, HMO, Degradation, -Community, -Treatment, -Time)

## Correct names of HMO structures
hmo3$HMO <- gsub('LNFP.III', "LNFP3", hmo3$HMO)
hmo3$HMO <- gsub('LNFP.II', "LNFP2", hmo3$HMO)
hmo3$HMO <- gsub('LNFP.I', "LNFP1", hmo3$HMO)
hmo3$HMO <- gsub('X2.FL', "2'FL", hmo3$HMO)
hmo3$HMO <- gsub('X3FL', "3-FL", hmo3$HMO)
hmo3$HMO <- gsub('X3.SL', "3'SL", hmo3$HMO)
hmo3$HMO <- gsub('X6.SL', "6'SL", hmo3$HMO)

## Calculate percent degradation
hmo4 <- spread(hmo3, Time, Degradation)
hmo4$'0_2' <- hmo4$'0'
hmo4[,4:10] <- sapply(hmo4[,4:10],as.numeric)
hmo4[,4:9] <- 100 - ((hmo4[,4:9] / hmo4[,10])*100)
hmo4$'0_2' <- NULL

## Reorganize dataframe
hmo5 <- gather(hmo4, Time, Degradation, -Community, -Treatment, -HMO)
hmo5 <- transform(hmo5, Time = as.numeric(Time))
hmo5_pHMOs <- subset(hmo5, hmo5$Treatment == "pHMOs")
hmo5_2FL <- subset(hmo5, hmo5$Treatment == "2FL" & hmo5$HMO == "2'FL")

hmo5_pHMOs$HMO <- factor(hmo5_pHMOs$HMO, levels = c("2'FL", "3'SL", "3-FL", "6'SL", "DFLac","DFLNT",  "DFLNH", "DSLNH", "DSLNT","FDSLNH", "FLNH", "LNFP1", "LNFP2", "LNFP3","LNH", "LNnT",  "LNT", "LSTb", "LSTc", "Fuc", "Sia", "SUM"), ordered = TRUE)

## Plot degradation curves
p <- ggplot(hmo5_pHMOs, aes(x = Time, y = Degradation, group = Community)) + 
  geom_line(aes(color = Community), size = 1) +
  geom_point(aes(color = Community), size = 1.5) +
  facet_wrap(~HMO, scales = "free", ncol = 8) +
  theme_bw() +
  ylab("Percent (%) degraded") +
  xlab("Time (h)") +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, 110)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA), breaks = c(0,10,20,30)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 10, colour = "black"),
        legend.key.size = unit(0.5, 'cm'),
        legend.key.height= unit(0.5, 'cm')) +
  scale_color_manual(values=c("#FF69B4","#800080","#0000FF","#008B8B","#008000","#FF8C00","#FF0000"))

p
topptx(file = "Community.HMOs.pptx", append = TRUE, width = 12, height = 4.5, units = "cm")

q <- ggplot(hmo5_2FL, aes(x = Time, y = Degradation, group = Community)) + 
  geom_line(aes(color = Community), size = 1) +
  geom_point(aes(color = Community), size = 2) +
  facet_wrap(~HMO, scales = "free") +
  theme_bw() +
  ylab("Percent (%) degraded") +
  xlab("Time (h)") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 10),
        strip.text.x = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  scale_color_manual(values=c("#FF69B4","#800080","#0000FF","#008B8B","#008000","#FF8C00","#FF0000"))

q
topptx(file = "Community.HMOs.pptx", append = TRUE, width = 1.5, height = 1.5, units = "cm")
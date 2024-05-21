##################################################################
## Metabonomic data analysis for Renwick et al. 2024
##################################################################

##############################
############################## Section 1: Load libraries and import and tidy data
##############################

## Clear workspace / variables
rm(list=ls())
graphics.off()

setwd("D:/")

## Load libraries
library(ggplot2)
library(eoffice)
library(ggfortify)
library(tidyverse)
library(pheatmap)
library(stats)
library(plyr)

## Load metabonomic datasets
metab <- read.csv("Renwick et al. 2024 Metabonomics.csv")

## Load metadata
metadata <- read.csv("Renwick et al. 2024 Metadata.csv")

## Correct metabolite and HMO names
names(metab)[names(metab) == "X2.Piperidinone"   ] <- "2-Piperidinone" 
names(metab)[names(metab) == "X3.Phenylpropionate"   ] <- "3-Phenylpropionate"   
names(metab)[names(metab) == "X4.Hydroxyphenylacetate"   ] <- "4-Hydroxyphenylacetate"  
names(metab)[names(metab) == "X4.Aminobenzoate"   ] <- "4-Aminobenzoate"   
names(metab)[names(metab) == "X5.Aminopentanoate" ] <- "5-Aminopentanoate"
names(metab)[names(metab) == "Beta.Alanine"  ] <- "Beta Alanine"   
names(metab)[names(metab) == "Nicotinic.acid"   ] <- "Nicotinic acid"
names(metab)[names(metab) == "p.Cresol" ] <- "p-Cresol"

metadata$Treatment <- gsub("2FL", "2'FL", metadata$Treatment, fixed=TRUE)




##############################
############################## Section 2: Normalize data and perform principal component analysis at T60 to evaluate the effect of Treatment
############################## 

## Transfer time information to data file
metab$Time <- metadata$Time[match(metab$NMR_file, metadata$NMR_file)]

## Subset samples to T60
metab60 <- subset(metab, metab$Time == "60")
metab60$Time <- NULL

## Autoscale function
autoscale <- function(data) {
  sample_classes <- data[, 1]
  x <- data[, 2:dim(data)[2]]
  x.sc <- scale(x, center = T, scale = T)
  x.sc <- cbind(sample_classes, x.sc)
  x.sc
}

## Normalize data by auto scaling
metab60_auto <- as.data.frame(autoscale(metab60))
metab60_auto[,2:ncol(metab60_auto)] <- sapply(metab60_auto[,2:ncol(metab60_auto)],as.numeric)

## Transfer sample information to data file
metab60_auto$Sample <- metadata$Sample[match(metab60_auto$sample_classes, metadata$NMR_file)]
metab60_auto$Community <- metadata$Community[match(metab60_auto$Sample, metadata$Sample)]
metab60_auto$Treatment <- metadata$Treatment[match(metab60_auto$Sample, metadata$Sample)]
metab60_auto$Time <- metadata$Time[match(metab60_auto$Sample, metadata$Sample)]
metab60_auto$sample_classes <- NULL

## Perform principle component analysis
pca <- subset(metab60_auto, select = -c(Sample, Time, Treatment, Community))
pca_res <- prcomp(pca, scale. = TRUE)

p <- autoplot(pca_res, data = metab60_auto, colour = "Treatment", size = 5) +
  theme_bw() +
  stat_ellipse(geom = "polygon", lwd = 1, aes(fill = Treatment), alpha = 0.25) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 20))

p

topptx(file="Metabonomics.pptx", append = TRUE, width = 6, height = 4, units = "cm")

## Label points with community
p2 <- p + geom_text(aes(label = Community))
p2
topptx(file="Metabonomics.pptx", append = TRUE, width = 6, height = 4, units = "cm")


##############################
############################## Section 3: Perform multiple ANOVAs followed by Tukey's test
##############################

## Change from wide to long format

metab3 <- gather(metab60_auto, Metabolite, Concentration, -Sample, -Time, -Treatment, -Community)

## Run multiple two-way ANOVAs and paste into dataframe
Metabolites <- unique(metab3$Metabolite)
## Create empty dataframe
TUKEYS <- data.frame(matrix(ncol = 6,nrow = 0, 
                            dimnames = list(NULL, c("Comparison", "Treatment.diff","Treatment.lwr","Treatment.upr","Treatment.p.adj","Metabolite"))))


for (i in 1:length(Metabolites)) {
  print(Metabolites[i])
  metab4 <- subset(metab3, metab3$Metabolite == Metabolites[i])
  res.aov <- aov(Concentration ~ Treatment, data = metab4)
  print(summary(res.aov))
  TukeyHSD(res.aov, conf.level=.95)
  
}

## Here are metabolites significantly affected by treatment
SigMet <- c("Acetate","Alanine","Ethanol","Glycolate","Histidine","Indole","Isoleucine","Lactate","Leucine","Malonate","Methionine","Methylamine","Phenylalanine","Propanol","Succinate","Threonine","Thymine","Trimethylamine","Tyrosine","Uracil","Valine","p-Cresol")
metab5 <- subset(metab3, metab3$Metabolite %in% SigMet)

## Create empty dataframe
tuc3 <- data.frame(matrix(ncol = 6,nrow = 0, 
                            dimnames = list(NULL, c("Comparison", "Treatment.diff","Treatment.lwr","Treatment.upr","Treatment.p.adj","Metabolite"))))

## Run Tukey's
for (i in 1:length(SigMet)) {
  metab6 <- subset(metab5, metab5$Metabolite == SigMet[i])
  res.aov <- aov(Concentration ~ Treatment, data = metab6)
  tuc <- (TukeyHSD(res.aov, which = "Treatment"))
  tuc2 <- as.data.frame(tuc[1:1])
  tuc2$Metabolite <- SigMet[i]
  tuc2 <- tibble::rownames_to_column(tuc2, "Comparison")
  tuc3 <- rbind(tuc3, tuc2)
}

## BH correction of p-values
tuc3$BH <- p.adjust(tuc3$Treatment.p.adj, method = "BH")

write.csv(tuc3, "Significant_metabolites.csv")


##############################
############################## Section 4: Generate heatmap
##############################

metab7 <- metab60
metab7$Sample <- metadata$Sample[match(metab7$NMR_file, metadata$NMR_file)]
row.names(metab7) <- metab7$Sample
metab7 <- subset(metab7, select = -c(NMR_file, Sample))
metab7 <- as.data.frame(t(metab7))
metab7[,1:ncol(metab7)] <- sapply(metab7[,1:ncol(metab7)],as.numeric)

## Prepare annotation data for heatmap
anno <- subset(metab60_auto, select = c(Sample, Community, Treatment))
row.names(anno) <- anno$Sample
anno$Sample <- NULL

## Generate heatmap
pheatmap(metab7, scale = "row", annotation_col = anno,
         annotation_colors = list(Treatment = c(`No Treatment` = "#00BA38", `2'FL` = "#F8766D", pHMOs = "#619CFF"), Community = c(Child_1 = "#FF69B4", Child_2 = "#800080", Child_3 = "#0000FF", Child_4 = "#008B8B", Child_5 = "#008000", Child_6 = "#FF8C00", Child_7 = "#FF0000")),
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu")) (256)),
         fontsize = 8, fontsize_row = 8, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE)

topptx(file = "Metabonomics.pptx", append = TRUE, width = 8, height = 8, units = "cm")





##############################
############################## Section 5: Line graphs for each metabolite and community in supplementary data
##############################

## From wide to long
metab8 <- metab
metab8 <- gather(metab8, Metabolite, Concentration, -Time, -NMR_file)

## Transfer metadata details
metab8$Community <- metadata$Community[match(metab8$NMR_file, metadata$NMR_file)]
metab8$Treatment <- metadata$Treatment[match(metab8$NMR_file, metadata$NMR_file)]

## Order the treatments
metab8$Treatment <- factor(metab8$Treatment, levels = c("No Treatment", "2'FL", "pHMOs"))

### Adding error bars to the graphs
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

metab9 <- data_summary(metab8, varname = "Concentration", groupnames = c("Community", "Metabolite", "Treatment", "Time"))

list <- unique(metab9$Metabolite)

plot1 <- subset(metab9, metab9$Metabolite %in% list[1:24])
plot2 <- subset(metab9, metab9$Metabolite %in% list[24:47])


p <-  ggplot(plot1, aes(x = Time, y = Concentration, group = Treatment, color = Treatment)) + 
  geom_line(size = 0.3) +
  geom_point(size = 0.5, show.legend = FALSE) +
  facet_grid(Metabolite~Community, scales = "free") +
  geom_errorbar(aes(ymin = Concentration-sd, ymax = Concentration+sd), width=10, size=0.7,
                position=position_dodge(0.05)) +
  labs(y = "Concentration (mM)") +
  theme_bw() +
  scale_color_manual(values=c("#00BA38","#F8766D","#619CFF")) +
  scale_x_time(name = "Time (h)", 
               limits = c(0, 60),
               breaks = c(0,20,40,60),
               labels = c("0","20","40","60")) +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 3.5),
        legend.position = "none")

p
ggsave(file = "Metab1.png", plot = p, width=4, height=12)

topptx(file = "Metabonomics.pptx", append = TRUE, width = 4, height = 12, units = "cm")


q <-  ggplot(plot2, aes(x = Time, y = Concentration, group = Treatment, color = Treatment)) + 
  geom_line(size = 0.3) +
  geom_point(size = 0.5, show.legend = FALSE) +
  facet_grid(Metabolite~Community, scales = "free") +
  geom_errorbar(aes(ymin = Concentration-sd, ymax = Concentration+sd), width=10, size=0.7,
                position=position_dodge(0.05)) +
  labs(y = "Concentration (mM)") +
  theme_bw() +
  scale_color_manual(values=c("#00BA38","#F8766D","#619CFF")) +
  scale_x_time(name = "Time (h)", 
               limits = c(0, 60),
               breaks = c(0,20,40,60),
               labels = c("0","20","40","60")) +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 6),
        strip.text.y = element_text(size = 3.5))
        #legend.position = "none")

q
ggsave(file = "Metab3.png", plot = q, width=4, height=11.5)
graphics.off()


topptx(file = "Metabonomics.pptx", append = TRUE, width = 4, height = 11.5, units = "cm")




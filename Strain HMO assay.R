#########################################################################################
## Growth curve and HMO degradation profile analysis for Renwick et al. 2024 
#########################################################################################

##############################
############################## Section 1: Load libraries and import and tidy data
##############################

## Clear workspace / variables
rm(list=ls())
graphics.off()

setwd("D:/")

## Load libraries
library(stringr)
library(rstatix)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(zoo)
library(factoextra)
library(ggfortify)
library(eoffice)
library(dplyr)

## Load optical density readings 
growth <- read.csv("Renwick et al. 2024 Strain OD readings.csv")

## Load HMO profiling data
hmo1 <- read.csv("Renwick et al. 2024 Strain HMO profiling 1.csv")
hmo2 <- read.csv("Renwick et al. 2024 Strain HMO profiling 2.csv")

## Simplify strain names
info <- read.csv("Strain_info.csv")

## Transfer simplified strain names to file with OD readings
growth$Simplified <- info$Simplified[match(growth$SpeciesI, info$SpeciesI)]
growth$SpeciesI <- growth$Simplified
growth$Simplified <- NULL

## Update taxonomic names
growth$SpeciesI <- gsub("Catabacter hongkongensis", "Christensenella hongkongensis", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("[Clostridium] aldenense", "Enterocloster aldenense", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("[Clostridium] bolteae", "Enterocloster bolteae", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("[Clostridium] celerecrescens", "Lacrimispora celerecrescens", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("[Clostridium] citroniae", "Enterocloster citroniae", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("[Clostridium] clostridioforme", "Enterocloster clostridioforme", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("[Clostridium] lavalense", "Enterocloster lavalense", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("Bacteroides vulgatus", "Phocaeicola vulgatus", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("Lactobacillus rhamnosus", "Lacticaseibacillus rhamnosus", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("[Clostridium] saccharogumia", "Erysipelatoclostridium saccharogumia", growth$SpeciesI, fixed=TRUE)
growth$SpeciesI <- gsub("Bacteroides dorei", "Phocaeicola dorei", growth$SpeciesI, fixed=TRUE)

info$Simplified <- gsub("Catabacter hongkongensis", "Christensenella hongkongensis", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("[Clostridium] aldenense", "Enterocloster aldenense", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("[Clostridium] bolteae", "Enterocloster bolteae", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("[Clostridium] celerecrescens", "Lacrimispora celerecrescens", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("[Clostridium] citroniae", "Enterocloster citroniae", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("[Clostridium] clostridioforme", "Enterocloster clostridioforme", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("[Clostridium] lavalense", "Enterocloster lavalense", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("Bacteroides vulgatus", "Phocaeicola vulgatus", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("Lactobacillus rhamnosus", "Lacticaseibacillus rhamnosus", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("[Clostridium] saccharogumia", "Erysipelatoclostridium saccharogumia", info$Simplified, fixed=TRUE)
info$Simplified <- gsub("Bacteroides dorei", "Phocaeicola dorei", info$Simplified, fixed=TRUE)




##############################
############################## Section 2: Generate all growth curves
##############################

### Continuing from Section 1

list <- c(1,2,3,4)

for (i in 1:length(list)) {
  
  N <- list(i)
  
  SpeciesIsolates <- sort(unique(growth$SpeciesI))
  
  if (N == 1) {
    SpeciesIsolates <- SpeciesIsolates[1:91]
  } else if (N == 2) {
    SpeciesIsolates <- SpeciesIsolates[92:182]
  } else if (N == 3) {
    SpeciesIsolates <- SpeciesIsolates[183:273]
  } else if (N == 4) {
    SpeciesIsolates <- SpeciesIsolates[274:330]
  }
  
  if (N == 4) {H <- 9.75} else {H <- 13}
  
  # Subset growth into each SpeciesIsolate
  growth_sub <- subset(growth, growth$SpeciesI %in% SpeciesIsolates)
  
  # Create new column for each Treatment-Replicate
  growth_sub$TR <- paste(growth_sub$SpeciesI, growth_sub$Treatment, growth_sub$Replicate, sep = ".")
  
  ## New loop to generate rolling mean for each TR
  TreatRep <- unique(growth_sub$TR)
  
  # New empty dataframe
  growth_sub_rm <- data.frame(matrix(ncol = 10, nrow = 0, dimnames = list(NULL, c("Plate","Well","SpeciesIsolate","Strain","Treatment","Replicate","time","OD","TR","Roll"))))
  
  for (j in 1:length(TreatRep)) {
    
    df <- subset(growth_sub, growth_sub$TR == TreatRep[j])
    df$Roll = rollmean(df$OD, 5, align = "left", fill = NA)
    growth_sub_rm <- rbind(growth_sub_rm, df)
    
  }
  
  ## Now to plot the growth curves
  p2 <- ggplot() + 
    geom_line(growth_sub_rm, mapping = aes(time, Roll, colour = Treatment, group = TR), size = 1.1) + 
    facet_wrap(~SpeciesI, scales = "free", ncol = 7) +
    scale_x_time(name = "Time (h)", 
                 limits = c(0, 50),
                 breaks = c(0, 12, 24, 36, 48),
                 labels = c("0","12","24","36","48")) +
    scale_y_continuous(name = "OD600") +
    theme_bw() + 
    theme(axis.text = element_text(size=10),
          plot.title = element_text(size = 25),
          axis.title = element_text(size=30),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 17),
          strip.text = element_text(size = 8)) +
    scale_color_manual(values = c('#00BA38','#619CFF')) +
    scale_shape_manual(values = c(19, 19))
  
  ggsave(file = paste0("Growth.curves", N, ".png"), plot = p2, width = 9.75, height = H, units = "in", path = "Curves/", scale = 2)
  
  graphics.off()
  
}




##############################
############################## Section 3: Calculate areas under the curves (AUCs), fold changes, and significance using multiple t-tests
##############################

## AUC calculation function  
EmpiricalAreaUnderCurve <- function(data_t, data_n) {
  x <- data_t
  y <- data_n
  n <- length(x)
  auc_e <- sum((x[2:n] - x[1:n-1]) * (y[2:n] + y[1:n-1]) /  2)
  return(auc_e)
}

growth$STR <- paste(growth$SpeciesI, growth$Treatment, growth$Replicate, sep = "_")
str <- sort(unique(growth$STR))

auc <-  data.frame(matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("STR", "AUC"))))

for (i in 1:length(str)) {
  growth_sub <- subset(growth, growth$STR == str[i])
  auc2 <-  data.frame(matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("STR", "AUC"))))
  auc2[1,1] <- str[i]
  auc2[1,2] <- EmpiricalAreaUnderCurve(data_t = growth_sub$time, data_n = growth_sub$OD)
  auc <- rbind(auc, auc2)
}

## Transfer details over to auc dataframe
auc$SpeciesI <- growth$SpeciesI[match(auc$STR, growth$STR)]
auc$Treatment <- growth$Treatment[match(auc$STR, growth$STR)]
auc$STR <- NULL

SpeciesIsolate <- unique(auc$SpeciesI)

auc_mean <-  data.frame(matrix(ncol = 4, nrow = 0, dimnames = list(NULL, c("Mean", "SD", "SpeciesI", "Treatment"))))

for (i in 1:length(SpeciesIsolate)) {
  auc2 <- subset(auc, auc$SpeciesI == SpeciesIsolate[i])
  auc2 <- auc2 %>% 
    dplyr::group_by(Treatment) %>%
    dplyr::summarise(Mean = mean(AUC))
  auc2$SpeciesI <- SpeciesIsolate[i]
  auc_mean <- rbind(auc_mean, auc2)
}

## Change from long to wide to calculate fold change
auc_fold <- spread(auc_mean, Treatment, Mean)

auc_fold$AUCf<- auc_fold$pHMOs / auc_fold$`No Treatment`
auc_fold$Log2AUCf <- log2(auc_fold$AUCf) 

## Multiple t-tests to determine significant changes
auc_t <- auc %>%
  group_by(SpeciesI) %>%
  t_test(AUC ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

## transfer fold change calculation over to dataframe with t-test results
auc_t$Log2AUCf <- auc_fold$Log2AUCf[match(auc_t$SpeciesI, auc_fold$SpeciesI)]
auc_final <- subset(auc_t, select = c(SpeciesI, Log2AUCf, p.adj, p.adj.signif))
auc_final$Log2AUCf <- round(auc_final$Log2AUCf, 2)
auc_final$p.adj <- round(auc_final$p.adj , 6)

write.csv(auc_final, "Strain.AUC.foldchange.csv")




##############################
############################## Section 4: Calculate mean highest OD, fold changes, and significance using multiple t-tests with FDR-adjusted p-values
##############################

##Calculate OD mean
od <- growth %>%
  dplyr::group_by(STR) %>%
  dplyr::slice_max(OD, n = 3) %>%
  dplyr::summarise(MeanODh = mean(OD))

## Transfer details over to od dataframe
od$SpeciesI <- growth$SpeciesI[match(od$STR, growth$STR)]
od$Treatment <- growth$Treatment[match(od$STR, growth$STR)]
od$STR <- NULL

SpeciesIsolate <- unique(od$SpeciesI)

od_mean <-  data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("Species_Isolate", "Treatment", "Mean_OD"))))

for (i in 1:length(SpeciesIsolate)) {
  od2 <- subset(od, od$SpeciesI == SpeciesIsolate[i])
  od2 <- od2 %>% 
    dplyr::group_by(Treatment) %>%
    dplyr::summarise(Mean = mean(MeanODh))
  od2$SpeciesI <- SpeciesIsolate[i]
  od_mean <- rbind(od_mean, od2)
}

## Change from long to wide to calculate fold change
od_fold <- spread(od_mean, Treatment, Mean)

od_fold$ODf<- od_fold$pHMOs / od_fold$`No Treatment`
od_fold$Log2ODf <- log2(od_fold$ODf) 

## Multiple t-tests to determine significant changes
od_t <- od %>%
  group_by(SpeciesI) %>%
  t_test(MeanODh ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

## transfer fold change calculation over to dataframe with t-test results
od_t$Log2ODf <- od_fold$Log2ODf[match(od_t$SpeciesI, od_fold$SpeciesI)]
od_final <- subset(od_t, select = c(SpeciesI, Log2ODf, p.adj, p.adj.signif))
od_final$Log2ODf <- round(od_final$Log2ODf, 2)
od_final$p.adj <- round(od_final$p.adj , 6)

write.csv(od_final, "Strain.OD.foldchange.csv")


##############################
############################## Section 5: Calculate HMO degradation profiles
##############################

## HMO profiles were collected in three sets. In the first set ("Strain_HMO_assay_profiling_1.csv") T0 samples were collected for every strain. In the second and third set (which have been combined into ("Strain_HMO_assay_profiling_2.csv")) T0 was collected for each plate. As such, calculations will be performed using different scripts.


## first set ("Strain_HMO_assay_profiling_1.csv")
## Change from wide to long and then back to wide
hmo1_2 <- gather(hmo1, HMO, Concentration, -Plate, -Strain, -time)
hmo1_2 <- spread(hmo1_2, time, Concentration)

## Calculate percent degradation
hmo1_2$Difference <- hmo1_2$`0` - hmo1_2$`48`
hmo1_2$Degradation <- ((hmo1_2$Difference/hmo1_2$`0`)*100)


## second and third set (combined into ("Strain_HMO_assay_profiling_2.csv")
## Create empty dataframe
hmo2_final <- hmo2[0,]

## Loop through plates and strains to create a dataframe similar to hmo1
Plates <- unique(hmo2$Plate)

for (a in 1:length(Plates)) {
  hmo2_sub <- subset(hmo2, hmo2$Plate == Plates[a])
  hmo2_0 <- subset(hmo2_sub, hmo2_sub$time == 0)
  hmo2_48 <- subset(hmo2_sub, hmo2_sub$time == 48)
  Strains <- unique(hmo2_48$Strain)
  
  for (b in 1:length(Strains)) {
    hmo2_48_strain <- subset(hmo2_48, hmo2_48$Strain == Strains[b])
    hmo2_02 <- hmo2_0
    hmo2_02$Strain <- hmo2_48_strain$Strain
    hmo2_combined <- rbind(hmo2_48_strain, hmo2_02)
    hmo2_final <- rbind(hmo2_final, hmo2_combined)
  }
}

## Change from wide to long and then back to wide
hmo2_final2 <- gather(hmo2_final, HMO, Concentration, -Plate, -Strain, -time)
hmo2_final2 <- spread(hmo2_final2, time, Concentration)
## Calculate percent degradation
hmo2_final2$Difference <- hmo2_final2$`0` - hmo2_final2$`48`
hmo2_final2$Degradation <- ((hmo2_final2$Difference/hmo2_final2$`0`)*100)

## Combine hmo1 and hmo2 data
hmo <- rbind(hmo1_2, hmo2_final2)

## Polish data 
hmo$HMO <- gsub('LNFP.III', "LNFP3", hmo$HMO)
hmo$HMO <- gsub('LNFP.II', "LNFP2", hmo$HMO)
hmo$HMO <- gsub('LNFP.I', "LNFP1", hmo$HMO)
hmo$HMO <- gsub('X2.FL', "2'FL", hmo$HMO)
hmo$HMO <- gsub('X3FL', "3-FL", hmo$HMO)
hmo$HMO <- gsub('X3.SL', "3'SL", hmo$HMO)
hmo$HMO <- gsub('X6.SL', "6'SL", hmo$HMO)

hmo$Degradation[hmo$Degradation > 100] <- 100
hmo$Degradation[hmo$Degradation < 0] <- 0

hmo <- subset(hmo, select = c(Strain, HMO, Degradation))

hmo$SpeciesI <- growth$SpeciesI[match(hmo$Strain, growth$Strain)]
hmo$Phylum <- info$Phylum[match(hmo$SpeciesI, info$Simplified)]

## Generate heatmaps for all taxa within each phyla for supplementary data
HMO_phylum <- subset(hmo, hmo$Phylum == "Firmicutes")

list <- unique(HMO_phylum$SpeciesI)
list <- sort(list)
list2 <- list[1:82]
list3 <- list[83:164]
list4 <- list[165:204]

HMO_phylum4 <- subset(HMO_phylum, HMO_phylum$SpeciesI %in% list4)

p <- ggplot(HMO_phylum4, aes(x=HMO, y=SpeciesI, fill=Degradation)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = '#0c65fa') +
  theme_bw() +
  #scale_fill_continuous(limits=c(0, 100)) +
  scale_y_discrete(limits = rev(unique(sort(HMO_phylum4$SpeciesI)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p

topptx(file = "HMO Degradation.pptx", append = TRUE, width = 8, height = 30, units = "cm")

## Determine top HMO degrading strains for final heatmap
hmo2 <- spread(hmo, HMO, Degradation)
hmo2$total <- rowSums(hmo2[,c(4:25)])
hmo3 <- hmo2 %>% 
  arrange(desc(total)) %>% slice(1:100)

hmo4 <- subset(hmo, hmo$SpeciesI %in% unique(hmo3$SpeciesI))

p <- ggplot(hmo4, aes(x=HMO, y=SpeciesI, fill=Degradation)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = '#0c65fa') +
  theme_bw() +
  scale_y_discrete(limits = rev(unique(sort(hmo$SpeciesI)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p

topptx(file = "HMO Degradation.pptx", append = TRUE, width = 10, height = 30, units = "cm")



metrics <- cbind(od_final, auc_final)

metrics2 <- subset(metrics, metrics$SpeciesI %in% unique(hmo3$SpeciesI))
metrics3 <- subset(metrics2, select = c(SpeciesI, Log2AUCf, Log2ODf))
metrics3 <- gather(metrics3, Metric, Foldchange, -SpeciesI)

p <- ggplot(metrics3, aes(x=Metric, y=SpeciesI, fill=Foldchange)) +
  geom_tile(color = "white",
            lwd = 1.3,
            linetype = 1) +
  coord_fixed() +
 # geom_text(aes(label = Foldchange), color = "black", size = 1) +
  scale_fill_gradient2(low = "#147d41", mid = "white", high = '#c61a11', midpoint = 0) +
 # scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  theme_bw() +
  scale_y_discrete(limits = rev(unique(sort(hmo$SpeciesI)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p

topptx(file = "HMO Degradation.pptx", append = TRUE, width = 8, height = 30, units = "cm")



##############################
############################## Section 6:Determine most degraded HMO structure
##############################

## Determine most degraded HMO structure
for (c in hmo2[,4:23]) {
  print(sum(c))
}



##############################
############################## Section 7: Spearman correlations among AUC and ODh fold changes and HMO degradation profiles
##############################

hmo4 <- hmo

hmo4$Log2AUCf <- auc_final$Log2AUCf[match(hmo4$SpeciesI, auc_final$SpeciesI)]
hmo4$AUC.sig <- auc_final$p.adj.signif[match(hmo4$SpeciesI, auc_final$SpeciesI)]

hmo4$Log2ODf <- od_final$Log2ODf[match(hmo4$SpeciesI, od_final$SpeciesI)]
hmo4$OD.sig <- od_final$p.adj.signif[match(hmo4$SpeciesI, od_final$SpeciesI)]

spear <- spread(hmo4, HMO, Degradation)
row.names(spear) <- spear$SpeciesI
spear <- subset(spear, select = -c(Strain, SpeciesI, Phylum, AUC.sig, OD.sig))

## Perform Spearman correlation
res <- rcorr(as.matrix(spear), type = "spearman")

## Plot to visualize significance
corrplot(res$r, type = "upper", order = "hclust", hclust.method = "ward.D2",
         p.mat = res$p, sig.level = 0.05, insig = "blank", method = 'color')

## Function for formatting the correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

## Format correlation matrix
res2 <- flattenCorrMatrix(res$r, res$P)
res2$p <- format(res2$p, scientific = TRUE)

## Subseting to the significant and strong correlations
res3 <- subset(res2, res2$cor >= 0.3 | res2$cor <= -0.3 | res2$p <= 0.05)

write.csv(res3, "Spearman_correlations.csv", row.names = FALSE)





##############################
############################## Section 8: Investigate strain-level heterogeneity with kmeans clustering
##############################

## Transfer genus information
k <- spread(hmo4, HMO, Degradation)
k$Species <- info$Species[match(k$SpeciesI, info$Simplified)]

## Subset to include only genera with three or more strains
count <- as.data.frame(table(k$Species))
count <- subset(count, count$Freq >= 3)
k2 <- subset(k, k$Species %in% unique(count$Var1))

## Subset to perform kmeans using only growth curve data
k_gc <- subset(k2, select = c(SpeciesI, Log2AUCf, Log2ODf))

# ## Check normality of growth curve data
# shapiro.test(k_gc$Log2AUCf)
# shapiro.test(k_gc$Log2ODf)
# 
# ## Normalize growth curve data
# for (i in 2:ncol(k_gc)) {
#   k_gc[,i] <- ((k_gc[,i] - mean(k_gc[,i])) / sd(k_gc[,i]) )
# }

## Determine optimal number of clusters for growth curve data
kmean_withinss <- function(k) {
  cluster <- kmeans(k_gc[,c(2:ncol(k_gc))], k, nstart = 100)
  return (cluster$tot.withinss)
}

# Set maximum cluster 
max_k <- 25

## Run algorithm over a range of k 
wss <- sapply(2:max_k, kmean_withinss)
elbow <- data.frame(2:max_k, wss)

## Plot an elbow diagram
ggplot(elbow, aes(x = X2.max_k, y = wss)) +
      geom_point(size = 5) +
      geom_line(size = 2.5) 

#### Optimal number of clusters looks like 6 clusters

## Perform kmeans clustering
k_gc_clusters <- kmeans(k_gc[,c(2:3)], 6, nstart = 100)
k_gc_clusters$size

## Visualize clustering
p <- fviz_cluster(k_gc_clusters, data = k_gc[,c(2:(ncol(k_gc)))], 
             geom = c("point"), 
             ellipse.type = "euclid",
             ggtheme = theme_bw(),
             main = "")

p

topptx(file="Kmeans.pptx", append = TRUE, width = 6, height = 4, units = "cm")

## Summarize into a dataframe
k_final <- as.data.frame(k_gc_clusters$cluster)
k_final <- cbind(k_final, k_gc)
k_final$Species <- k$Species[match(k_final$SpeciesI, k$SpeciesI)]

## Perform kmeans clustering using only HMO degradation data
k_hmo <- subset(k2, select = -c(Strain, Fuc, Sia, SUM, Phylum, Log2AUCf, Log2ODf, AUC.sig, OD.sig, Species))

## Check normality of growth curve data
shapiro.test(k_hmo$`2'FL`)
ggdensity(k_hmo$`2'FL`)

## Normalize data
for (i in 2:(ncol(k_hmo))) {
  k_hmo[,i] <- ((k_hmo[,i] - mean(k_hmo[,i])) / sd(k_hmo[,i]) )
}

## Determine optimal number of clusters for growth curve data
kmean_withinss <- function(k) {
  cluster <- kmeans(k_hmo[,c(1:(ncol(k_hmo)-1))], k, nstart = 100)
  return (cluster$tot.withinss)
}

# Set maximum cluster 
max_k <- 25

## Run algorithm over a range of k 
wss <- sapply(2:max_k, kmean_withinss)
elbow <- data.frame(2:max_k, wss)

## Plot an elbow diagram
ggplot(elbow, aes(x = X2.max_k, y = wss)) +
  geom_point(size = 5) +
  geom_line(size = 2.5) 

#### Optimal number of clusters looks like 6 clusters

## Perform kmeans clustering
k_hmo_clusters <- kmeans(k_hmo[,c(2:(ncol(k_hmo)))], 6, nstart = 100)
k_hmo_clusters$size

## Visualize clustering
p <- fviz_cluster(k_hmo_clusters, data = k_hmo[,c(2:(ncol(k_hmo)))], 
             geom = c("point"), 
             ellipse.type = "euclid",
             ggtheme = theme_bw(),
             main = "")

p

topptx(file="Kmeans.pptx", append = TRUE, width = 6, height = 4, units = "cm")

## Summarize into a dataframe
k_hmo2 <- as.data.frame(k_hmo_clusters$cluster)
k_final <- cbind(k_final, k_hmo2)
k_final2 <- subset(k_final, select = c(Species, SpeciesI, `k_gc_clusters$cluster`, `k_hmo_clusters$cluster`))

write.csv(k_final2, "Strain.kmeans.clustering.csv", row.names = FALSE)

## After currating kmean clustering groups, generate figures
k3 <- read.csv("Strain.kmeans.clustering_currated.csv")

k4 <- gather(k3, Metric, Cluster, -Species, -SpeciesI)

k.gc <- subset(k4, k4$Metric == "GC.clusters")
k.deg <- subset(k4, k4$Metric == "Deg.clusters")

gc <- ggplot(k.gc, aes(x=Metric, y=SpeciesI, fill=Cluster)) +
  geom_tile(color = "white",
            lwd = 1.3,
            linetype = 1) +
  coord_fixed() +
  scale_fill_gradientn(colours = c("#660000","#CC0000","#FF9900","#FFFF66")) +
  theme_bw() +
  scale_y_discrete(limits = rev(unique(sort(k.gc$SpeciesI)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gc

topptx(file = "HMO Degradation.pptx", append = TRUE, width = 8, height = 30, units = "cm")


deg <- ggplot(k.deg, aes(x=Metric, y=SpeciesI, fill=Cluster)) +
  geom_tile(color = "white",
            lwd = 1.3,
            linetype = 1) +
  coord_fixed() +
  scale_fill_gradientn(colours = c("#000066","#9933FF","#00CCFF")) +
  theme_bw() +
  scale_y_discrete(limits = rev(unique(sort(k.gc$SpeciesI)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

deg

topptx(file = "HMO Degradation.pptx", append = TRUE, width = 8, height = 30, units = "cm")
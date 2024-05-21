#############################################################################
### 16S rRNA gene sequencing data analysis for Renwick et al. 2024 
#############################################################################


##############################
############################## Section 1: Load libraries and import data, remove low abundance and prevalence ASVs, and replace zeros
##############################

### Clear workspace / variables
rm(list=ls())
graphics.off()

setwd("D:/")

## Load libraries
library(Biostrings)
library(msa)
library(ape)
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(eoffice)
library(Maaslin2)
library(dplyr)
library(vegan)

## Load ASV table
reads <- read.csv("Renwick et al. 2024 16S rRNA reads.csv")

## Load cell count data
cellcounts <- read.csv("Renwick et al. 2024 Community Cell Counts.csv")

## Load metadata
metadata <- read.csv("Renwick et al. 2024 Metadata.csv")


##############################
############################## Section 2: Assess baseline alpha-diversity using Chao1 and Shannon diversity measures
##############################

## Create phyloseq object with unfiltered data

## Organize data
reads_a <- subset(reads, select = -c(Kingdom, Phylum, Class, Order, Family, Genus, sequences))
row.names(reads_a) <- reads_a$ASV
reads_a <- subset(reads_a, select = unique(metadata$Sample))
asv_mat <- as.matrix(reads_a)

## Create matrix of taxonomy
tax <- subset(reads, select = c(Kingdom, Phylum, Class, Order, Family, Genus, ASV))
row.names(tax) <- tax$ASV
tax_mat <- as.matrix(tax)

row.names(metadata) <- metadata$Sample

## Create phyloseq object
OTU <- otu_table(asv_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, samples)

## Subset phyloseq object to T0 samples
ps0 <- subset_samples(ps, Time == 0)

## Calculate alpha-diversity metrics
alpha <- estimate_richness(ps, measures = c("Chao1", "Shannon"))
alpha <- subset(alpha, select = c(Chao1, Shannon))
alpha$Community <- metadata$Community[match(rownames(alpha), rownames(metadata))]
alpha2 <- gather(alpha, Metric, Value, -Community)

p <- ggplot(alpha2, aes(x = Community, y = Value, fill = Community)) +
  geom_boxplot(outlier.size = 0.7, lwd = 0.7) +
  geom_jitter(color = "black", size = 0.7, alpha = 0.6) +
  theme_bw() +
  facet_wrap(~Metric, scales = "free") +
  labs(y = "Alpha Diversity Metric") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        legend.position = "none",
        strip.text.x = element_text(size = 20)) +
  scale_fill_manual(values = c("#FF69B4","#800080","#0000FF","#008B8B","#008000","#FF8C00","#FF0000")) 

p

topptx(file = "Metataxonomics.pptx", append = TRUE, width = 1.25, height = 1, units = "cm")




##############################
############################## Section 3: Filer data and normalize genus level relative abundances to cell counts to produce absolute abundance profiles
##############################

## Transform reads to relative abundance
reads2 <- reads %>% tibble::column_to_rownames("ASV")
reads3 <- subset(reads2, select = -c(Kingdom, Phylum, Class, Order, Family, Genus, sequences))
reads_ra <- apply(reads3, 2, function(x){x/sum(x)})
reads_ra <- as.data.frame(reads_ra)

## Filter low abundance taxa
abund = 0.00001 # 0.001% relative abundance
prev = 0.01    # 1% prevalence 
reads4 <- as.data.frame(reads_ra[apply(reads_ra, 1, function(x) mean(x >= abund) >= prev),], check.names = FALSE)

## Subset to HMO-treated samples up to T60
reads5 <- subset(reads4, select = unique(metadata$Sample))

## Aggregate at the genus-level
reads5$Genus <- reads$Genus[match(rownames(reads5), reads$ASV)]
reads6 <- aggregate(. ~ Genus, reads5, sum)
reads6 <- reads6 %>% tibble::column_to_rownames("Genus")

## Normalize relative abundances to cell counts to produce absolute abundances
absol_g <- dplyr::mutate(reads6, across(1:ncol(reads6), ~.x*cellcounts[cellcounts$Sample == cur_column(),2]))




##############################
############################## Section 4: Maaslin2 anaylsis
##############################

metadata$Treatment <- gsub("No Treatment", "NT", metadata$Treatment, fixed=TRUE)
row.names(metadata) <- metadata$Sample

## Run Maaslin2
fit_data_ixn = Maaslin2(input_data     = absol_g, 
                        input_metadata = metadata, 
                        normalization  = "none",
                        transform = "none",
                        output         = "maaslin2", 
                        fixed_effects  = c("Treatment"),
                        reference = c("Treatment,NT"),
                        random_effects = c("Community"),
                        plot_scatter = TRUE)

## Extract Maaslin2 results and subset to genera significantly impacted by pHMOs
maaslin_results <- fit_data_ixn$results
sig <- subset(maaslin_results, maaslin_results$qval < 0.05)

## Correct taxonomic names
sig$feature <- gsub("Escherichia.Shigella", "Escherichia-Shigella", sig$feature, fixed=TRUE)
sig$feature <- gsub("UCG.", "UCG-", sig$feature, fixed=TRUE)
sig$feature <- gsub("X.", "[", sig$feature, fixed=TRUE)
sig$feature <- gsub("..", "] ", sig$feature, fixed=TRUE)
sig$feature <- gsub(".", " ", sig$feature, fixed=TRUE)

write.csv(sig, "Significant_maaslin_results.csv")

## Determine genera positively and negatively impacted by pHMOs
sig_pHMOs <- subset(sig, sig$value == "pHMOs")

sig_pHMOs <- sig_pHMOs |> mutate(effect = case_when(
  coef >= 1 ~ "positive",
  coef <= 1 ~ "negative",
  TRUE ~ "undetermined"
))

q <- ggplot(sig_pHMOs, aes(x = reorder(feature, +coef), y = coef, fill = effect)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(y = "MaAsLin2 Coefficient", x = "Genus") +
  theme(axis.title = element_text(size = 10),
        legend.position = "none") +
  scale_fill_manual(values = c("#808080","#619CFF"))

q

topptx(file = "Metataxonomics.pptx", append = TRUE, width = 4, height = 7, units = "cm")


## Generate box plots 
## Prepare absolute abundance data for plotting box plots
absol3 <- absol_g %>% 
  rownames_to_column() %>%
  gather(Sample, Abundance, -rowname)

absol3$Treatment <- metadata$Treatment[match(absol3$Sample, metadata$Sample)]
absol3$Treatment <- gsub("2FL", "2'FL", absol3$Treatment, fixed=TRUE)
absol3$Treatment <- gsub("NT", "No Treatment", absol3$Treatment, fixed=TRUE)
absol3$Treatment <- factor(absol3$Treatment, levels = c("No Treatment","2'FL","pHMOs"), ordered = TRUE)

## Fix genera names
absol3$rowname <- gsub("Escherichia.Shigella", "Escherichia-Shigella", absol3$rowname, fixed=TRUE)
absol3$rowname <- gsub("UCG.", "UCG-", absol3$rowname, fixed=TRUE)
absol3$rowname <- gsub("X.", "[", absol3$rowname, fixed=TRUE)
absol3$rowname <- gsub("..", "] ", absol3$rowname, fixed=TRUE)
absol3$rowname <- gsub(".", " ", absol3$rowname, fixed=TRUE)

## Select genera for plotting
sig_pHMOs_pos <- subset(sig_pHMOs, sig_pHMOs$coef > 0)
list_pos <- unique(sig_pHMOs_pos$feature)
sig_pHMOs_pos_plot <- subset(absol3, absol3$rowname %in% list_pos)

sig_pHMOs_neg <- subset(sig_pHMOs, sig_pHMOs$coef < 0)
list_neg <- unique(sig_pHMOs_neg$feature)
sig_pHMOs_neg_plot <- subset(absol3, absol3$rowname %in% list_neg)

## Create box plots

## All in supplementary data
r <- ggplot(sig_pHMOs_pos_plot, aes(x = Treatment, y = Abundance, fill = Treatment)) +
  facet_wrap(~rowname, scales = "free", ncol = 11) +
  geom_boxplot(outlier.size = 0.2, lwd = 0.5) +
  theme_bw() +
  labs(x = "Treatment", y = "Absolute Abundance (CFU/mL)", color = "Treatment") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 6, angle = 25, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 7),
        legend.position = "none") +
  scale_fill_manual(values = c("#00BA38","#F8766D","#619CFF")) 

r 

topptx(file = "Metataxonomics.pptx", append = TRUE, width = 13, height = 5, units = "cm")

s <- ggplot(sig_pHMOs_neg_plot, aes(x = Treatment, y = Abundance, fill = Treatment)) +
  facet_wrap(~rowname, scales = "free", ncol = 11) +
  geom_boxplot(outlier.size = 0.2, lwd = 0.5) +
  theme_bw() +
  labs(x = "Treatment", y = "Absolute Abundance (CFU/mL)", color = "Treatment") +
  theme(axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 6, angle = 25, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 7),
        legend.position = "none") +
  scale_fill_manual(values = c("#00BA38","#F8766D","#619CFF")) 

s

topptx(file = "Metataxonomics.pptx", append = TRUE, width = 8, height = 1.5, units = "cm")

## Top 16 positively affected genera
sig_pHMOs_pos <- sig_pHMOs %>% arrange(desc(coef)) %>% slice(1:16)
list_pos <- unique(sig_pHMOs_pos$feature)
sig_pHMOs_pos_plot <- subset(absol3, absol3$rowname %in% list_pos)

t <- ggplot(sig_pHMOs_pos_plot, aes(x = Treatment, y = Abundance, fill = Treatment)) +
  facet_wrap(~rowname, scales = "free", ncol = 4) +
  geom_boxplot(outlier.size = 0.2, lwd = 0.5) +
  theme_bw() +
  labs(x = "Treatment", y = "Absolute Abundance (CFU/mL)", color = "Treatment") +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 7)) +
  scale_fill_manual(values = c("#00BA38","#F8766D","#619CFF")) 

t 

topptx(file = "Metataxonomics.pptx", append = TRUE, width = 8, height = 7, units = "cm")




##############################
############################## Section 4: Generate PCoA plots of UniFrac distances to compare effects of Treatment and Community
##############################

## Create phyloseq object of genus-level absolute abundance data with phylogenetic tree

## Organize dataframes
asv_mat <- as.matrix(absol_g)

## Create matrix of taxonomy
tax <- subset(reads, select = c(Kingdom, Phylum, Class, Order, Family, Genus))
tax <- tax[!duplicated(tax$Genus),]
row.names(tax) <- tax$Genus
tax_mat <- as.matrix(tax)

## Prepare metadata
row.names(metadata) <- metadata$Sample

## Create phyloseq object
OTU <- otu_table(asv_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(metadata)
ps_HMO <- phyloseq(OTU, TAX, samples)

## Create phylogenetic tree
seq <- reads[!duplicated(reads$Genus),]
seqmat <- as.matrix(subset(seq, select = c(Genus, sequences)))
seqset <- DNAStringSet(seq$sequences)
names(seqset) <- seqmat[,1]
alignedSet <- msa(seqset,"ClustalW")
alignedBin <- as.DNAbin(alignedSet)
ddist <- dist.dna(x = alignedBin)
phy_tree <- njs(ddist)

## Add phylogenetic tree
ps_HMO2 <- merge_phyloseq(ps_HMO, phy_tree)

## Remove T0 samples
ps_HMO3 <- subset_samples(ps_HMO2, Time %in% c(12,24,36,48,60))

## Perform ordination
uni_w <- ordinate(ps_HMO3, method = "PCoA", "unifrac", weighted = TRUE)

## Colours
Com_colours <- c("#FF69B4","#800080","#0000FF","#008B8B","#008000","#FF8C00","#FF0000")
Tre_colours <- c("#F8766D","#00BA38","#619CFF")

## Plot weighted UniFrac distances on PCoA

## Coloured by community
plot_ordination(ps_HMO3, uni_w, type = "samples", color = "Community") +
  geom_point(size = 3) +
  scale_colour_manual(values = Com_colours) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15)) +
  labs(title = "Weighted UniFrac") 

topptx(file = "Metataxonomics.pptx", append = TRUE, width = 5, height = 4, units = "cm")

## Coloured by treatment
plot_ordination(ps_HMO3 , uni_w, type = "samples", color = "Treatment") +
  geom_point(size = 3) +
  scale_colour_manual(values = Tre_colours) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15)) +
  labs(title = "Weighted UniFrac") 

topptx(file = "Metataxonomics.pptx", append = TRUE, width = 5, height = 4, units = "cm")

## Coloured by time
plot_ordination(ps_HMO3, uni_w, type = "samples", color = "Time") +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15)) +
  labs(title = "Weighted UniFrac") 

topptx(file = "Metataxonomics.pptx", append = TRUE, width = 5, height = 4, units = "cm")


## Calculate permanova

# all non-0 time points
uni_w_dist <- distance(ps_HMO3, method = "unifrac", weighted = TRUE)
adonis2(uni_w_dist ~ sample_data(ps_HMO3)$Community)
adonis2(uni_w_dist ~ sample_data(ps_HMO3)$Treatment)
adonis2(uni_w_dist ~ sample_data(ps_HMO3)$Time)

# check other time points
ps_HMO_sub <- subset_samples(ps_HMO3, Time == 60)
uni_w_dist <- distance(ps_HMO_sub, method = "unifrac", weighted = TRUE)
adonis2(uni_w_dist ~ sample_data(ps_HMO_sub)$Community)
adonis2(uni_w_dist ~ sample_data(ps_HMO_sub)$Treatment)



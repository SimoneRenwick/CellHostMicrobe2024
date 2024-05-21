##################################################################
### DADA2 processing of fastq files for Renwick et al. 2024
##################################################################

##############################
############################## Section 1: Process fastq files with DADA2
##############################

### For each set of fastq runs, perform the following functions and save a separate .RDS file

## Load libraries
library(dada2)
library(tidyverse)
library(stringr)

## Set path to fastq files
path <- "E:/BaseSpace/"
list.files(path)

## Read in fastq files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

## Filter and trim 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)

## Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

## Sample inferences
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

## Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

## Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, "E:/BaseSpace/seqtab1.rds")


##############################
############################## Section 1: Merge .rds files and annotate taxonomy using SILVA database file
##############################

### Continue from Section 1

## Merge multiple runs
st1 <- readRDS("C:/Users/srenwick/Sequencing files/seqtab1.rds")
st2 <- readRDS("C:/Users/srenwick/Sequencing files/seqtab2.rds")
st3 <- readRDS("C:/Users/srenwick/Sequencing files/seqtab3.rds")
st4 <- readRDS("C:/Users/srenwick/Sequencing files/seqtab4.rds")
st5 <- readRDS("C:/Users/srenwick/Sequencing files/seqtab5.rds")
st6 <- readRDS("C:/Users/srenwick/Sequencing files/seqtab6.rds")
st7 <- readRDS("C:/Users/srenwick/Sequencing files/seqtab7.rds")

st.all <- mergeSequenceTables(st1, st2, st3, st4, st5, st6, st7)

## Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/srenwick/Sequencing files/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa2 <- as.data.frame(taxa)

## Update "NA"s in taxa file to unclassified
for (a in 2:6) {
  for (i in 1:nrow(taxa2)) {
    if (is.na(taxa2[i,a])) { 
      if (grepl("_unclassified", taxa2[i,a-1])) {
        taxa2[i,a] <- taxa2[i,a-1]
      }
      else {
        taxa2[i,a] <- paste(taxa2[i,a-1], "_unclassified", sep = "")
      }
    }
  }
}

## Create final ASV table
seqtab.nochim2 <- as.data.frame(t(seqtab.nochim))
reads <- merge(seqtab.nochim2, taxa2, by = "row.names", all = TRUE)
names(reads)[names(reads) == "Row.names"] <- "sequences"
reads <- reads %>% relocate(sequences, .after = last_col())

reads$sequence_length <- str_count(reads$sequences)
df <- as.data.frame(table(reads$sequence_length))
reads <- subset(reads, sequence_length > 289 & sequence_length < 304)
reads$sequence_length <- NULL

reads$ASV1 <- "ASV_"
reads$ASV2 <- 1:nrow(reads)
reads$ASV <- paste(reads$ASV1, reads$ASV2, sep = "")
reads <- subset(reads, select = -c(ASV1, ASV2))
reads <- reads %>% relocate(ASV)

reads_new_names <- as.data.frame(t(reads))
row.names(reads_new_names) <- gsub("-", ".", row.names(reads_new_names))
row.names(reads_new_names) <- gsub("S.3.2FL", "S3.2FL", row.names(reads_new_names))
row.names(reads_new_names) <- gsub("S.3.pHMO", "S3.pHMO", row.names(reads_new_names))
row.names(reads_new_names) <- gsub("S.3.Con", "S3.Con", row.names(reads_new_names))

reads2 <- as.data.frame(t(reads_new_names))
write.csv(reads2, "Renwick et al. 2024 16S rRNA reads.csv", row.names = FALSE)

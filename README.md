# CellHostMicrobe2024
R scripts for performing bioinformatics analyses in Renwick et al. 2024 Cell Host Microbe - Modulating the developing gut microbiota with 2’-fucosyllactose and pooled human milk oligosaccharides.

Description of scripts:

"Community HMO profiling.R" was used to process microbial community HMO degradation data ("Renwick et al. 2024 Community HMO profiling.csv") to produce Figure 1.

"Fastq processing using DADA2.R" was used to process 16S rRNA gene sequencing fastq files available on NCBI database with BioProject ID PRJNA1094576, and SRA accession SAMN40701351 to SAMN40701728 using the DADA2 package (Callahan et al. 2016 Nat Methods) to produce the file "Renwick et al. 2024 16S rRNA reads.csv".

"Metataxonomics.R" was used to filter and normalize the 16S rRNA reads ("Renwick et al. 2024 16S rRNA reads.csv"), calculate alpha and beta diversities, and conduct multivariate analysis using the MaAsLin2 package to produce Figure 2 and Table S1. 

"Metabonomics.R" was used to process, normalize, and analyze the NMR metabonomic data ("Renwick et al. 2024 Metabonomics.csv") to produce Figures 3 and S1 and Table S2.

"Strain HMO assay.R" was used to generate and analyze growth curves for each strain treated with pHMOs in monoculture assays from OD readings ("Renwick et al. 2024 Strain OD readings") to produce Figures 4 and S2 and Table S4. Furthermore, this script was used to combine and analyze the strain HMO degradation data ("Renwick et al. 2024 Strain HMO profiling 1.csv" and "Renwick et al. 2024 Strain HMO profiling 2.csv") to produce Figures 5 and S3 and determine any significant correlations, producing Figure S4 and "Spearman_correlations.csv". Lastly, this script also clustered strain growth and HMO degradation data using the Kmeans method to produce Figures 6 and S5 and Table S5.

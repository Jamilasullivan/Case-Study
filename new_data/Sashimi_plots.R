## PACKCAGES ###################################################################

<<<<<<< HEAD
#install.packages("BiocManager")
#BiocManager::install("ggbio")
#BiocManager::install("GenomicRanges")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("GenomicAlignments")
#install.packages("ggplot2")

library(ggbio)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
=======
install.packages("BiocManager")
BiocManager::install("ggbio")
BiocManager::install("GenomicRanges")
install.packages("ggplot2")

library(ggbio)
library(GenomicRanges)
>>>>>>> d6aebcb34edb1fea96fb5cbb9875fab72afcea6f
library(ggplot2)
library(Rsamtools)
library(rtracklayer)

## SET WORKING DIRECTORY #######################################################

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data")

## READING IN BAM FILES FOR 4 DIFFERENT SAMPLES ################################

bam_files <- c("SRR12650992_sorted.bam","SRR12650993_sorted.bam","SRR12650994_sorted.bam","SRR12650995_sorted.bam")

# Read each BAM file
alignments_list <- lapply(bam_files, function(bam_file) {
  readGAlignments(bam_file)
})

## READ IN GFF FILE FROM MISO ##################################################

gff_file <- "genes.gff"
gff_data <- import.gff(gff_file)

## CREATING SASHIMI PLOT #######################################################

# Combine alignments into a single object
combined_alignments <- do.call(c, alignments_list)

# Create the sashimi plot
sashimi_plot <- autoplot(combined_alignments, 
                         layout = "sashimi",
                         aes(color = factor(rep(1:4, each = length(combined_alignments) / 4)))) +
  scale_color_manual(values = c("blue", "red", "green", "purple")) + 
  theme_minimal()

print(sashimi_plot)






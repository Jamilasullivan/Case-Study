## LOAD PACKAGES ###############################################################

#BiocManager::install("Gviz")

library(Gviz)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)

## SET WORKING DIRECTORY #######################################################

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data")

## EXON SKIPPING (SE) RESULTS ##################################################

txdb <- makeTxDbFromGFF("Mus_musculus.GRCm39.113.gtf", format = "gtf")
seqlevels(txdb)
txdb_chr19 <- keepSeqlevels(txdb, "19", pruning.mode = "coarse")
seqlevels(txdb_chr19)

# Load BAM files
bam_file_92 <- readGAlignments("SRR12650992_sorted.bam")
bam_file_93 <- readGAlignments("SRR12650993_sorted.bam")
bam_file_94 <- readGAlignments("SRR12650994_sorted.bam")
bam_file_95 <- readGAlignments("SRR12650995_sorted.bam")

# Create DataTrack objects from the BAM files
covTrack92 <- DataTrack(range = bam_file_92, genome = "hg38", type = "coverage", chromosome = "chr19", name = "Control")
covTrack93 <- DataTrack(range = bam_file_93, genome = "hg38", type = "coverage", chromosome = "chr19", name = "Control")
covTrack94 <- DataTrack(range = bam_file_94, genome = "hg38", type = "coverage", chromosome = "chr19", name = "Control")
covTrack95 <- DataTrack(range = bam_file_95, genome = "hg38", type = "coverage", chromosome = "chr19", name = "Control")


region_chr <- "chr19"
region_start <- 5816185
region_end <- 5816344

## CREAT GENE AND AXIS TRACKS ##################################################

?GeneRegionTrack
geneTrack <- GeneRegionTrack(txdb_chr19,
                             chromosome = region_chr,
                             start = region_start,
                             end = region_end,
                             transcriptAnnotation = "symbol",
                             name = "Gene")

axisTrack <- GenomeAxisTrack()

covTrack1 <- DataTrack(range = "CTR_sample1.bam",
                       genome = "hg38",
                       type = "coverage",
                       chromosome = region_chr,
                       name = "Control")

covTrack2 <- DataTrack(range = "AIRP_sample1.bam",
                       genome = "hg38",
                       type = "coverage",
                       chromosome = region_chr,
                       name = "AIRP")

plotTracks(list(axisTrack, geneTrack, covTrack1, covTrack2),
           from = region_start,
           to = region_end,
           transcriptAnnotation = "symbol",
           main = "Sashimi-like Plot (R-only)")

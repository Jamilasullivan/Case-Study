## installing and loading necessary packages ###################################

# only run the following installation lines on first application of the script
#install.packages("BiocManager")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("apeglm", force = TRUE)
#BiocManager::install("EnhancedVolcano")

# Loading packages

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)    # Mouse annotation package
library(AnnotationDbi)
library(EnhancedVolcano)

## DATA ORGANISATION ###########################################################

## set working directory

setwd(
  "C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data"
)

## load count data

counts <- read.csv("featurecounts.csv", row.names = 1) # read in counts file and set gene IDs as row names
duplicated(rownames(counts)) # check for duplicates
rownames(counts)[duplicated(rownames(counts))] # check duplicate names
#View(counts) # view counts table
head(counts) # get first columns of the table
dim(counts) # get count of rows

## load metadata info

metadata <- read.csv("Metadata2.csv", row.names = 1) # read in metadata with sample IDs as row names
metadata$Condition <- as.factor(metadata$Condition) # make condition a factor variable
levels(metadata$Condition) # get factor levels
metadata$Condition <- factor(metadata$Condition, levels = c("Control", "Case")) # order factor levels
metadata$Condition # get factors
head(metadata) # get first rows of the data frame
str(metadata) # get info on all

# view both

#View(counts)
#View(metadata)

# ordering row names in metadata to match counts data

metadata <- metadata[order(rownames(metadata)), ] # ordering the metadata rows by sample name to match the counts data
#View(metadata) # check the metadata

## DESeq2 ANALYSIS #############################################################

## create a DESeq object and import the sample information and count data

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Condition) # creating the DESeq object
dim(dds) # checking row numbers of DESeq object matches that of the original table

## set the reference for the treatment factor

dds$Condition <- factor(dds$Condition, levels = c("Control", "Case")) # setting factor levels (again. not needed)

## filter the genes

keep <- rowSums(counts(dds) >= 10) >= min(table(metadata$Condition)) # when set to 5 there was a lot of NAs for adjusted p value
dds <- dds[keep, ] # only keeping the rows that correlate to counts above 10
dim(dds) # checking row numbers
duplicated(rownames(dds)) #check for duplicates
rownames(dds)[duplicated(rownames(dds))] # get duplicates names

#View(dds)

## perform statistical tests to identify DEGs

dds <- DESeq(dds) # statistical analysis. Normalisation, dispersion estimation, fitting negative binomial model, Wald test for DEGs, multiple testing correction

# changing ensembl IDs to gene symbols in dds object

rownames(dds) # checking row names
ensembl_ids <- rownames(dds) # extracting row names

gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
) # changing ENSEMBL IDs to gene symbols

head(gene_symbols) # checking gene symbols
rownames(dds) <- ifelse(!is.na(gene_symbols), gene_symbols, rownames(dds)) # changing the row names to the gene symbols
head(rownames(dds)) # checking the row names

# viewing results

deseq_results <- results(dds) # extracting DEGs
deseq_results # visualising
dim(deseq_results) # checking row number

## Change DESeq object to an R object (dataframe)

deseq_results <- as.data.frame(deseq_results) # changes results into a data frame
class(deseq_results) # check it's a data frame
head(deseq_results) # get first 6 rows
#View(deseq_results)
sum(is.na(deseq_results$padj)) # check there is no NAs for adjusted p value

## order the results table by increasing p-value

deseq_results_ordered <- deseq_results[order(deseq_results$padj), ] # order results by adjusted p value
head(deseq_results_ordered) # get first 6 rows to check the above worked

duplicated(rownames(deseq_results_ordered)) # checks for duplicated row names
rownames(deseq_results_ordered)[duplicated(rownames(deseq_results_ordered))] # shows what duplicates are

head(deseq_results_ordered) # checking the data

anyDuplicated(rownames(deseq_results_ordered)) # this tells me there are no duplicates in the row names

## make some queries about specific genes

deseq_results_ordered["Cd44", ] # there
deseq_results_ordered["Siglecf", ] # there. The gene of interest.

## extract the most differentially expressed genes due to the treatment
## select genes with a significant change in gene expression (adjusted p value <0.05)
## and log2foldchange <1 and >1

## Step 1: filter based on adjusted p value

deseq_filtered <- deseq_results_ordered %>% filter(deseq_results_ordered$padj < 0.05) # gives genes below 0.05 adjusted p value

## Step 2: filter based on fold changes

deseq_filtered <- deseq_filtered %>% filter(abs(deseq_filtered$log2FoldChange) > 1) # gives genes above a log2fold change of 1

dim(deseq_results_ordered) # checks data
dim(deseq_filtered) # checks data

## Step 3: make queries

## Save the deseq result. We will save both the original data(res) and the filtered one(hits)

write.csv(deseq_results, "DESeq_all_deseq_results_sifglecf.csv") # saves DESeq results as a csv
write.csv(deseq_filtered, "DESeq_filtered_deseq_results_siglecf.csv") # saves filtered DESeq results as a csv

## save normalised read counts

deseq_normalised_counts <- counts(dds, normalized = TRUE) # creates a new object where the counts data is normalised
head(deseq_normalised_counts) # checks data
write.csv(deseq_normalised_counts, "DESeq_normalised_counts_siglecf.csv") # saves normalised counts as a csv

## VISUALISATION ###############################################################

## Dispersion estimates

plotDispEsts(dds) # plots the dispersion estimates

## PCA

# variance stabilisation transformation

deseq_vsd <- vst(dds, blind = F) # creates an object of variance stabilised data

# use transformed values to create a PCA plot

plotPCA(deseq_vsd, intgroup = c("Condition"), pcsToUse = 1:2) # PC1 explains 98% of variance
plotPCA(deseq_vsd, intgroup = c("Condition"), pcsToUse = 2:3) # PC3 explains 0% of variance

## Heatmaps
# R packgage: pheatmap

# Heatmap of sample-to-sample distance matrix (with clustering) based on the normalised counts.

# generate the distance matrix

deseq_sampleDists <- dist(t(assay(deseq_vsd)))
deseq_sampleDistMatrix <- as.matrix(deseq_sampleDists)
colnames(deseq_sampleDistMatrix)
deseq_sampleDistMatrix

# set colour scheme

colours <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)

# generate the heatmap
samples_heatmap <- pheatmap(
  deseq_sampleDistMatrix,
  clustering_distance_rows = deseq_sampleDists,
  clustering_distance_cols = deseq_sampleDists,
  col = colours,
  border_color = "transparent",
  cluster_cols = T,
  cluster_rows = T,
  fontsize = 8
)

# Heatmap of log transformed normalised counts. We will use the top 10 genes.

colours2 <- colorRampPalette(brewer.pal(9, "Reds"))(255)

top_hits <- deseq_results_ordered[order(deseq_results_ordered$padj), ][1:10, ]
top_hits <- row.names(top_hits)
top_hits

deseq_rld <- rlog(dds, blind = FALSE)

pheatmap(
  assay(deseq_rld)[top_hits, ],
  cluster_rows = T,
  show_rownames = T,
  cluster_cols = T,
  fontsize = 8,
  color = colours2
)

# testing the top 30 genes

top_30 <- deseq_results[order(deseq_results$padj), ][1:30, ]
top_30 <- row.names(top_30)
top_30

pheatmap(
  assay(deseq_rld)[top_30, ],
  show_rownames = T,
  fontsize = 8,
  color = colours2,
  show_colnames = F
)

# MA plot

plotMA(dds, ylim = c(-2, 2), alpha = 0.05)

design(dds) # checks the design
results(dds) # tells me what is being compared
resLFC <- lfcShrink(dds, coef = "Condition_Case_vs_Control", type = "apeglm") # remove the noise
plotMA(resLFC, ylim = c(-2, 2), alpha = 0.05)
print(plotMA(resLFC, ylim = c(-2, 2), alpha = 0.05))
str(resLFC)

## other volcano

deseq_results$symbol <- mapIds(
  org.Mm.eg.db,
  keys = rownames(deseq_results),
  keytype = "SYMBOL",
  column = "SYMBOL"
)
deseq_filtered
deseq_results

EnhancedVolcano(
  resLFC,
  x = "log2FoldChange",
  y = "padj",
  lab = deseq_results$symbol,
  selectLab = top_hits,
  title = NULL,
  subtitle = NULL,
  labSize = 5,
  FCcutoff = 2,
  pCutoff = 0.05
) # looks wrong because the negative binomial distribution model in DESeq cannot handle only 4 samples. Limma does a much better job for fewer samples.

resLFC

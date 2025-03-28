## installing and loading necessary packages ###################################

# only run the next two lines on first application of the script
#install.packages("BiocManager")
#BiocManager::install("org.Mm.eg.db")

# Loading packages

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)  # Mouse annotation package
library(AnnotationDbi)

## DATA ORGANISATION ###########################################################

## set working directory

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data")

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

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~Condition) # creating the DESeq object
dim(dds) # checking row numbers of DESeq object matches that of the original table

## set the reference for the treatment factor

dds$Condition <- factor(dds$Condition, levels = c("Control","Case")) # setting factor levels (again. not needed)

## filter the genes 

keep <- rowSums(counts(dds) >=10) >= min(table(metadata$Condition)) # when set to 5 there was a lot of NAs for adjusted p value
dds <- dds[keep,] # only keeping the rows that correlate to counts above 10
dim(dds) # checking row numbers
duplicated(rownames(dds)) #check for duplicates
rownames(dds)[duplicated(rownames(dds))] # get duplicates names

#View(dds)

## perform statistical tests to identify DEGs

dds <- DESeq(dds) # statistical analysis. Normalisation, dispersion estimation, fitting negative binomial model, Wald test for DEGs, multiple testing correction. 
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

deseq_results_ordered <- deseq_results[order(deseq_results$padj),] # order results by adjusted p value
head(deseq_results_ordered) # get first 6 rows to check the above worked

# gene names from ensemble to mouse gene symbols

gene_symbols <- mapIds(
  org.Mm.eg.db, 
  keys = rownames(deseq_results_ordered), 
  column = "SYMBOL", 
  keytype = "ENSEMBL",
  multiVals = "first"
) # mapping ensemble IDs to gene symbols

print(gene_symbols) # looking at mapped genes

duplicated(rownames(deseq_results_ordered)) # checks for duplicated row names
rownames(deseq_results_ordered)[duplicated(rownames(deseq_results_ordered))] # shows what duplicates are

deseq_results_ordered$gene_symbol <- ifelse(
  !is.na(gene_symbols), 
  gene_symbols, 
  rownames(deseq_results_ordered)
) # adding the gene symbols as an extra column

head(deseq_results_ordered) # checking the column has been added

deseq_results_ordered$gene_symbol <- make.names(deseq_results_ordered$gene_symbol, unique = TRUE) # making the duplicated gene names unique

head(deseq_results_ordered) # checking the data 

make_unique_with_underscore <- function(x) {
  # Initialize a vector to store unique names
  unique_names <- character(length(x))
  
  # Create a named counter to track occurrences
  counter <- list()
  
  for (i in seq_along(x)) {
    name <- x[i]
    
    # Check if the name has been seen before
    if (name %in% names(counter)) {
      counter[[name]] <- counter[[name]] + 1
      unique_names[i] <- paste0(name, "_", counter[[name]])
    } else {
      counter[[name]] <- 1  # Start counting from 1
      unique_names[i] <- name
    }
  }
  
  return(unique_names)
}# function to name duplicates with a _ between the gene name and the number

test <- c("test_run", "test.run", "testrun", "testrun", "testrun") # testing the function
test <- make_unique_with_underscore(test) 
test[grep("\\_", test)] 

deseq_results_ordered$gene_symbol <- make_unique_with_underscore(deseq_results_ordered$gene_symbol) # using function to change duplicates to being separated with '_' rather than '.'. This is because some mouse gene symbol contain .s so it makes it difficult to tell the difference
head(deseq_results_ordered)

deseq_results_ordered$gene_symbol[grep("\\_", deseq_results_ordered$gene_symbol)] # check for the new values

anyDuplicated(deseq_results_ordered$gene_symbol) # this tells me there are no duplicates in the gene symbols
anyDuplicated(rownames(deseq_results_ordered)) # this tells me there are no duplicates in the row names (ENSEMBL IDs)

rownames(deseq_results_ordered) <- deseq_results_ordered$gene_symbol # changes the row names of the dataframe to the gene_symbol column
head(deseq_results_ordered) # checks the results

identical(deseq_results_ordered$gene_symbol, rownames(deseq_results_ordered)) # tells me that all of these row names match perfectly to the gene_symbol column

deseq_results_ordered <- deseq_results_ordered[, !colnames(deseq_results_ordered) %in% "gene_symbol"] # removes the gene_symbol column
head(deseq_results_ordered) # checks the results

## make some queries about specific genes

deseq_results_ordered["Cd44",] # there
deseq_results_ordered["Nrxn1",] # not there 
deseq_results_ordered["Dscam",] # not there
deseq_results_ordered["Mblm1",] # not there
deseq_results_ordered["Sox2",] # not there

## extract the most differentially expressed genes due to the treatment
## select genes with a significant change in gene expression (adjusted p value <0.05)
## and log2foldchange <1 and >1

## Step 1: filter based on adjusted p value

filtered <- deseq_results_ordered %>% filter(deseq_results_ordered$padj < 0.05) # gives genes below 0.05 adjusted p value

## Step 2: filter based on fold changes

filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1) # gives genes above a log2fold change of 1

dim(deseq_results_ordered) # checks data
dim(filtered) # checks data

## Step 3: make queries

## Save the deseq result. We will save both the original data(res) and the filtered one(hits) 

write.csv(deseq_results, "all_deseq_results_sifglecf.csv") # saves DESeq results as a csv
write.csv(filtered, "filtered_deseq_results_siglecf.csv") # saves filtered DESeq results as a csv

## save normalised read counts 

normalised_counts <- counts(dds,normalized = TRUE) # creates a new object where the counts data is normalised
head(normalised_counts) # checks data 
write.csv(normalised_counts, "normalised_counts_siglecf.csv") # saves normalised counts as a csv

## VISUALISATION ###############################################################

## Dispersion estimates

plotDispEsts(dds)

## PCA

# variance stabilisation transformation

vsd <- vst(dds, blind = F)

# use transformed values to create a PCA plot

plotPCA(vsd, intgroup=c("Condition"), pcsToUse = 1:2) # PC1 explains 98% of variance
plotPCA(vsd, intgroup=c("Condition"), pcsToUse = 2:3) # PC3 explains 0% of variance

## Heatmaps
# R packgage: pheatmap

# Heatmap of sample-to-sample distance matrix (with clustering) based on the normalised counts.

# generate the distance matrix

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix)
sampleDistMatrix

# set colour scheme

colours <- colorRampPalette(rev(brewer.pal(9,"Reds")))(255)

# generate the heatmap
samples_heatmap <- pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
                            clustering_distance_cols = sampleDists, col = colours, border_color = "transparent", cluster_cols = T, cluster_rows = T, fontsize = 8)

# Heatmap of log transformed normalised counts. We will use the top 10 genes. 

colours2 <- colorRampPalette(brewer.pal(9,"Reds"))(255)

top_hits <- deseq_results[order(deseq_results$padj),][1:10,]
top_hits <- row.names(top_hits)
top_hits

rld <- rlog(dds,blind=FALSE)

pheatmap(assay(rld)[top_hits,], cluster_rows = T, show_rownames = T, cluster_cols = T, fontsize = 8, color = colours2)

annot_info <- as.data.frame(colData(dds)[,c('group','time_point')])
pheatmap(assay(rld)[top_hits,], cluster_rows = T, show_rownames = T, cluster_cols = T,
         annotation_col = annot_info, fontsize = 8, color = colours2,)

pheatmap(assay(rld)[top_hits,], cluster_rows = T, show_rownames = T, cluster_cols = T,
         annotation_col = annot_info, fontsize = 8, color = colours2, kmeans_k = 4) # k-means cluster analysis of the same data as the heatmap above

pheatmap(assay(rld)[top_hits,], cluster_rows = T, show_rownames = T, cluster_cols = T,
         annotation_col = annot_info, fontsize = 8, color = colours2, cutree_rows = 4, cutree_cols = 5, show_colnames = FALSE)

# testing the top 30 genes

top_30 <- deseq_results[order(deseq_results$padj),][1:30,]
top_30 <- row.names(top_30)
top_30

pheatmap(assay(rld)[top_30,], cluster_rows = T, show_rownames = T, cluster_cols = T,
         annotation_col = annot_info, fontsize = 8, color = colours2, cutree_rows = 4, cutree_cols = 5, show_colnames = F)
# Heatmap of Z scores. We will use the top 10 genes. 

cal_z_score <- function(x) {(x-mean(x))/sd(x)}
zscore_all <- t(apply(normalised_counts,1,cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset, cluster_rows = T, cluster_cols = T, annotation_col = annot_info, fontsize = 8, cutree_rows = 4, cutree_cols = 5, show_colnames = FALSE, color = colours2)

# MA plot 

plotMA(dds,ylim = c(-2,2), alpha = 0.05)

resLFC <- lfcShrink(dds,coef="group_infected_vs_mock", type = "apeglm") # remove the noise
plotMA(resLFC,ylim=c(-2,2), alpha = 0.05)

## other volcano 

deseq_results$symbol <- mapIds(org.Hs.eg.db, keys = rownames(deseq_results), keytype = "SYMBOL", column = "SYMBOL")
filtered
deseq_results

EnhancedVolcano(resLFC, x = "log2FoldChange", y = "padj", lab = deseq_results$symbol, selectLab = top_hits, title = NULL, subtitle = NULL, labSize = 5, FCcutoff = 2, pCutoff = 0.05)






  
  
  
  
  
  
  
  
  
  


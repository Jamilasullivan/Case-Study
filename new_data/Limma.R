## PACKAGE INSTALLATION AND LOADING ############################################

## only run the following two commands on first run of the script

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")

## loading packages

library(edgeR)
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(org.Mm.eg.db)    # Mouse annotation package
library(AnnotationDbi)
library(pheatmap)
library(ComplexHeatmap)

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

## SETUP #######################################################################

## create DGEList object

limma_results <- DGEList(counts) # identifies differentially expressed genes#
head(limma_results)
limma_results
#View(results)

## PREPROCESSING ###############################################################

limma_results <- calcNormFactors(limma_results) # calculates normalisation factors 
limma_results$samples
limma_results

## filtering results

cutoff <- 10 
drop <- which(apply(cpm(limma_results),1,max) < cutoff)
limma_filtered <- limma_results[-drop,]
dim(limma_filtered)
limma_filtered

#set sample names

sample_names <- colnames(counts)
sample_names

## Voom transformation to make data suitable for Limma

limma_design <- model.matrix(~ metadata$Condition) # creates the design of the analysis
limma_filtered <- voom(limma_filtered, limma_design, plot = TRUE) # log transformation and estimates mean-variance trends

## fitting the linear model

limma_fit <- lmFit(limma_filtered, limma_design) # fits linear mode for each gene

## emperical Baysian smoothing

limma_fit <- eBayes(limma_fit)

## EXTRACTING DIFFERENTIALLY EXPRESSED GENES #####################

limma_results <- topTable(limma_fit, coef = 2, adjust.method = "BH", number = Inf) # extracts a table of the top ranked genes from a linear model fit
#View(limma_results)

limma_significant_genes <- limma_results[limma_results$adj.P.Val < 0.05, ] # genes with an adjusted p value of below 0.05

write.csv(limma_results, "limma_differential_expression_results.csv", row.names = TRUE) # saves the results in a csv file

## RENAMING ENSEMBL IDS AS GENE SYMBOLS ###########################

limma_ensembl_ids <- rownames(limma_significant_genes)

limma_gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = limma_ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
) # changing ENSEMBL IDs to gene symbols

limma_gene_symbols <- as.data.frame(limma_gene_symbols)
#View(limma_gene_symbols)

make_unique_with_underscore <- function(x) {
  # Initialize a vector to store unique names
  unique_names <- character(length(x))
  
  # Create a named counter to track occurrences
  counter <- list()
  
  for (i in seq_along(x)) {
    name <- x[i]
    
    # Skip NA values by leaving them unchanged
    if (is.na(name)) {
      unique_names[i] <- NA
    } else {
      # Check if the name has been seen before
      if (name %in% names(counter)) {
        counter[[name]] <- counter[[name]] + 1
        unique_names[i] <- paste0(name, "_", counter[[name]])
      } else {
        counter[[name]] <- 1  # Start counting from 1
        unique_names[i] <- name
      }
    }
  }
  
  return(unique_names)
} # this function reassigns any duplicate values with a numbr following an underscore. e.g. *_1, *_2, *_3

limma_gene_symbols <- make_unique_with_underscore(limma_gene_symbols$limma_gene_symbols) # carries out the function on the gene symbols of the data frame
limma_gene_symbols[grep("\\_", limma_gene_symbols)] # this checks what duplicates were found


head(limma_gene_symbols) # checking gene symbols
#View(limma_gene_symbols)
rownames(limma_significant_genes) <- ifelse(!is.na(limma_gene_symbols), limma_gene_symbols, rownames(limma_significant_genes)) # changing the row names to the gene symbols
#View(limma_significant_genes)

duplicated(rownames(limma_gene_symbols)) # check for duplicates
rownames(limma_significant_genes)[duplicated(rownames(limma_gene_symbols))] # check duplicate names

head(rownames(limma_significant_genes)) # checking the row names
#View(limma_results)
head(limma_significant_genes) 

## ordering significant genes

limma_significant_genes_ordered <- limma_significant_genes[order(limma_significant_genes$adj.P.Val),]

head(limma_significant_genes_ordered)
#View(limma_significant_genes_ordered)

## changing the name in the normal results to make a volcano plot from

limma_ensembl_ids <- rownames(limma_results)

limma_gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = limma_ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
) # changing ENSEMBL IDs to gene symbols

limma_gene_symbols <- as.data.frame(limma_gene_symbols)
#View(limma_gene_symbols)

make_unique_with_underscore <- function(x) {
  # Initialize a vector to store unique names
  unique_names <- character(length(x))
  
  # Create a named counter to track occurrences
  counter <- list()
  
  for (i in seq_along(x)) {
    name <- x[i]
    
    # Skip NA values by leaving them unchanged
    if (is.na(name)) {
      unique_names[i] <- NA
    } else {
      # Check if the name has been seen before
      if (name %in% names(counter)) {
        counter[[name]] <- counter[[name]] + 1
        unique_names[i] <- paste0(name, "_", counter[[name]])
      } else {
        counter[[name]] <- 1  # Start counting from 1
        unique_names[i] <- name
      }
    }
  }
  
  return(unique_names)
} # this function reassigns any duplicate values with a numbr following an underscore. e.g. *_1, *_2, *_3

limma_gene_symbols <- make_unique_with_underscore(limma_gene_symbols$limma_gene_symbols) # carries out the function on the gene symbols of the data frame
limma_gene_symbols[grep("\\_", limma_gene_symbols)] # this checks what duplicates were found

head(limma_gene_symbols) # checking gene symbols
#View(limma_gene_symbols)
rownames(limma_results) <- ifelse(!is.na(limma_gene_symbols), limma_gene_symbols, rownames(limma_results)) # changing the row names to the gene symbols
#View(limma_results)

duplicated(rownames(limma_gene_symbols)) # check for duplicates
rownames(limma_significant_genes)[duplicated(rownames(limma_gene_symbols))] # check duplicate names

head(rownames(limma_significant_genes)) # checking the row names
#View(limma_results)
head(limma_results)

## VISUALISING THE RESULTS #####################################################

## volcano plot

limma_results$significant <- limma_results$adj.P.Val < 0.05 & abs(limma_results$logFC) > 1 # filtering by adjusted p value and log2fold change

dev.off() # resetting the graphics system. Ensures that the following plot works.

write.csv(limma_results, "limma_differential_expression_results.csv", row.names = TRUE) # saves the results in a csv file

## VISUALISING THE RESULTS #####################################################

## VOLCANO PLOTS ###############################################################

limma_results$significant <- limma_results$adj.P.Val < 0.05 & abs(limma_results$logFC) > 1 # filtering by adjusted p value and log2fold change

limma_results # true and false column for significance of genes

dev.off() # resetting the graphics system. Ensures that the following plot works.

ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme_minimal()

## order genes and extract top 10

limma_results_ordered <- limma_results[order(limma_results$adj.P.Val),]

limma_top_hits <- limma_results_ordered[order(limma_results_ordered$adj.P.Val), ][1:200, ] # top 200 significant genes by adjusted p value
limma_top_hits <- row.names(limma_top_hits)
limma_top_hits

EnhancedVolcano(
  limma_results,
  x = 'logFC',
  y = 'adj.P.Val',
  lab = rownames(limma_results),
  selectLab = limma_top_hits,
  FCcutoff = 1,
  pCutoff = 0.05
)

limma_results

## HEATMAPS ####################################################################

class(limma_filtered)
limma_filtered <- as.data.frame(limma_filtered)

limma_ensembl_ids <- rownames(limma_filtered)

limma_gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = limma_ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
) # changing ENSEMBL IDs to gene symbols

limma_gene_symbols <- as.data.frame(limma_gene_symbols)
#View(limma_gene_symbols)

limma_gene_symbols <- make_unique_with_underscore(limma_gene_symbols$limma_gene_symbols) # carries out the function on the gene symbols of the data frame
limma_gene_symbols[grep("\\_", limma_gene_symbols)] # this checks what duplicates were found

head(limma_gene_symbols) # checking gene symbols
#View(limma_gene_symbols)
rownames(limma_filtered) <- ifelse(!is.na(limma_gene_symbols), limma_gene_symbols, rownames(limma_filtered)) # changing the row names to the gene symbols
#View(limma_results)

duplicated(rownames(limma_gene_symbols)) # check for duplicates
rownames(limma_filtered)[duplicated(rownames(limma_gene_symbols))] # check duplicate names

head(rownames(limma_filtered)) # checking the row names
#View(limma_results)
head(limma_filtered)

# Load required packages
library(limma)
library(pheatmap)
library(RColorBrewer)

# Get topTable from limma results
limma_results <- topTable(limma_fit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")

# Filter for significant genes (adjusted p-value < 0.05)
sig_genes <- limma_results[limma_results$adj.P.Val < 0.05, ]

# Select top 50 most significant genes
top_genes <- head(rownames(sig_genes), 30)

# Subset expression matrix by top genes
# (Make sure 'limma_filtered' contains log2 expression values and has matching rownames)
heatmap_matrix <- limma_filtered[top_genes, , drop = FALSE]

# Optional: Scale expression values across genes for better visualization
heatmap_matrix <- t(scale(t(heatmap_matrix)))

# Set up color palette
colours <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)

# Clear any open graphical devices (optional but prevents issues)
if (dev.cur() != 1) dev.off()

# Draw the heatmap
pheatmap(
  heatmap_matrix,
  col = colours,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10
)

## beginning of pathway analysis

library(clusterProfiler)
library(org.Mm.eg.db)

ego <- enrichGO(gene = rownames(limma_significant_genes), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP")
dotplot(ego)

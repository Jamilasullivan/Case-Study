## PACKAGE INSTALLATION AND LOADING ############################################

## only run the following two commands on first run of the script

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("limma")
#BiocManager::install("ComplexHeatmap")

## loading packages

library(edgeR)
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(org.Mm.eg.db)    # Mouse annotation package
library(AnnotationDbi)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)

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

expression_matrix <- limma_filtered$E
#View(expression_matrix)

## fitting the linear model

limma_fit <- lmFit(limma_filtered, limma_design) # fits linear mode for each gene

## emperical Baysian smoothing

limma_fit <- eBayes(limma_fit)

## EXTRACTING DIFFERENTIALLY EXPRESSED GENES ###################################

limma_results <- topTable(limma_fit, coef = 2, adjust.method = "BH", number = Inf) # extracts a table of the top ranked genes from a linear model fit. These are the results you need for volcano plots!
#View(limma_results)

limma_significant_genes <- limma_results[limma_results$adj.P.Val < 0.05 & abs(limma_results$logFC > 1), ] # genes with an adjusted p value of below 0.05

write.csv(limma_results, "limma_differential_expression_results.csv", row.names = TRUE) # saves the results in a csv file

write.csv(limma_significant_genes, "limma_significant_adjp0.05_logFC1.csv", row.names = TRUE) # saves the results in a csv file

################################################################################
################################################################################
################### CORRECT TO HERE ############################################
################################################################################
################################################################################

## ordering significant genes

limma_significant_genes_ordered <- limma_significant_genes[order(limma_significant_genes$adj.P.Val),]

head(limma_significant_genes_ordered)
#View(limma_significant_genes_ordered)

## RENAMING ENSEMBL IDS AS GENE SYMBOLS ###########################

limma_ensembl_ids <- rownames(limma_significant_genes_ordered)

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
rownames(limma_significant_genes_ordered) <- ifelse(!is.na(limma_gene_symbols), limma_gene_symbols, rownames(limma_significant_genes_ordered)) # changing the row names to the gene symbols
#View(limma_significant_genes)

duplicated(rownames(limma_gene_symbols)) # check for duplicates
rownames(limma_significant_genes)[duplicated(rownames(limma_gene_symbols))] # check duplicate names

head(rownames(limma_significant_genes_ordered)) # checking the row names
#View(limma_results)
head(limma_significant_genes_ordered) 

## VOLCANO PLOTS ###############################################################

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

################################################################################
######################## VISUALISING THE RESULTS ###############################
################################################################################

######################### VOLCANO PLOTS ########################################

limma_results$significant <- limma_results$adj.P.Val < 0.05 & abs(limma_results$logFC) > 1 # filtering by adjusted p value and log2fold change

limma_results # true and false column for significance of genes

dev.off() # resetting the graphics system. Ensures that the following plot works.

ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme_minimal()

############################ enhanced volcano plot #############################

## order genes and extract top 10 for labels

limma_results_ordered <- limma_results[order(limma_results$adj.P.Val),]

limma_top_hits <- limma_significant_genes_ordered[1:10, ] # top 30 significant genes by adjusted p value
limma_top_hits <- row.names(limma_top_hits)
limma_top_hits

dev.off()

EnhancedVolcano(
  limma_results,
  x = 'logFC',
  y = 'adj.P.Val',
  lab = rownames(limma_results),
  selectLab = limma_top_hits,
  labSize = 5,
  FCcutoff = 1,
  pCutoff = 0.05
)

############################## HEATMAPS ########################################

annot_info <- as.data.frame(metadata$Condition)

annotation_colors <- list(Group = c("Case" = "green", "Control" = "orange"))

dim(limma_significant_genes_ordered)

top_genes <- head(rownames(limma_significant_genes_ordered), 10)  # Change 30 to however many genes you want
top_genes

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

heatmap_matrix <- limma_filtered[top_genes, ]

heatmap_matrix <- as.matrix(heatmap_matrix)
#View(heatmap_matrix)

head(rownames(heatmap_matrix)) # checking the row names
#View(limma_results)
head(heatmap_matrix)

colours <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)

dev.off()

pheatmap(heatmap_matrix,
         col = colours,
         show_colnames = TRUE,
         show_rownames = TRUE,
         cutree_rows = 4,
         cutree_cols = 2,
         clustering_method = "complete",
         fontsize_row = 6,
         fontsize_col = 10,
         annotation_col = annot_info,
         annotation_colors = annotation_colors
         ) # creates a heatmap of the top 10 genes

########## heatmap of all genes. Visulaising the overall sprea #################

top_genes <- head(rownames(limma_significant_genes_ordered), 500)  # Change 1047 (all filtered significant genes) to however many genes you want.
top_genes

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

heatmap_matrix <- limma_filtered[top_genes, ]

heatmap_matrix <- as.matrix(heatmap_matrix)
#View(heatmap_matrix)

head(rownames(heatmap_matrix)) # checking the row names
#View(limma_results)
head(heatmap_matrix)

colours2 <- colorRampPalette(c("red", "black", "green"))(255)

pheatmap(heatmap_matrix,
         col = colours2,
         show_colnames = TRUE,
         show_rownames = FALSE,
         fontsize_row = 6,
         fontsize_col = 10,
         annotation_col = annot_info,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE, 
         cluster_rows = TRUE
) # this will take a while

## heatmap of z scores

z_scores <- t(scale(t(expression_matrix)))
z_scores

limma_ensembl_ids <- rownames(z_scores)

limma_gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = limma_ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
) # changing ENSEMBL IDs to gene symbols

limma_gene_symbols <- as.data.frame(limma_gene_symbols)
View(limma_gene_symbols)

limma_gene_symbols <- make_unique_with_underscore(limma_gene_symbols$limma_gene_symbols) # carries out the function on the gene symbols of the data frame
limma_gene_symbols[grep("\\_", limma_gene_symbols)] # this checks what duplicates were found

head(limma_gene_symbols) # checking gene symbols
#View(limma_gene_symbols)
rownames(z_scores) <- ifelse(!is.na(limma_gene_symbols), limma_gene_symbols, rownames(z_scores)) # changing the row names to the gene symbols
#View(limma_results)

duplicated(rownames(limma_gene_symbols)) # check for duplicates

zscore_subset <- z_scores[top_genes,]

pheatmap(heatmap_matrix,
         col = colours,
         show_colnames = TRUE,
         show_rownames = TRUE,
         cutree_rows = 4,
         cutree_cols = 2,
         clustering_method = "complete",
         fontsize_row = 8,
         fontsize_col = 10,
         annotation_col = annot_info,
         annotation_colors = annotation_colors
) # creates a heatmap of the top 10 genes

###################### BOX PLOT FOR SELECTED GENES #############################

gene <- "ENSMUSG00000039013" # this should be the gene you want to look at. It must match the rowname in the expression matrix

#view(expression_matrix) # unhash to check the structure of the expression matrix

"ENSMUSG00000039013" %in% rownames(expression_matrix)  # Checking the gene you want is in the expression matrix. TRUE means the gene exists.


gene_boxplot <- data.frame(Expression = expression_matrix[gene, ],
                 Condition = metadata$Condition) # creating the gene boxplot object

ggplot(gene_boxplot, aes(x = Condition, y = Expression)) +
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 3) +  # Customize boxplot
  labs(title = paste(gene, "Expression"),
       x = "Condition", 
       y = "Expression Level") + # in this case, this is by log2-transformed data
  theme_minimal() +  # Cleaner background
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and bold the title
    axis.text = element_text(size = 12),  # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    panel.grid.major = element_line(color = "gray90", size = 0.5)  # Subtle grid lines
  )


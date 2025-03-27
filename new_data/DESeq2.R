## installing and loading necessary packages ###################################

install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)  # Mouse annotation package
library(AnnotationDbi)

## set working directory

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA")

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
V#iew(metadata)

# ordering row names in metadata to match counts data

metadata <- metadata[order(rownames(metadata)), ] # ordering the metadata rows by sample name to match the counts data
#View(metadata) # check the metadata

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
    
    # Check if the name is duplicated
    if (name %in% names(counter)) {
      counter[[name]] <- counter[[name]] + 1
      unique_names[i] <- paste0(name, "_", counter[[name]])
    } else {
      counter[[name]] <- 0
      unique_names[i] <- name
    }
  }
  
  return(unique_names)
} # function to name duplicates with a _ between the gene name and the number

test <- c("help_me", "help.me", "helpme", "helpme", "helpme")
test <- make_unique_with_underscore(test)
test[grep("\\_", test)] 


deseq_results_ordered$gene_symbol <- make_unique_with_underscore(deseq_results_ordered$gene_symbol) # using function to change duplicates to being separated with '_' rather than '.'. This is because some mouse gene symbol contain .s so it makes it difficult to tell the difference
head(deseq_results_ordered)

deseq_results_ordered$gene_symbol[grep("\\_", deseq_results_ordered$gene_symbol)] # check for the new values

## make some queries 

## is MBD2 differentially expressed?

deseq_results["MBD2",]
deseq_results["ACE2",]
deseq_results["IL6",]

## extract the most differentially expressed genes due to the treatment
## select genes with a significant change in gene expression (adjusted p value <0.05)
## and log2foldchange <1 and >1

## Step 1: filter based on adjusted p value

filtered <- deseq_results %>% filter(deseq_results$padj < 0.05)

## Step 2: filter based on fold changes

filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1)

dim(deseq_results)
dim(filtered)

## Step 3: make queries




###############################################################################

## testing stuff

## variability in the data

summary(rowSums(counts(dds)))

summary(deseq_results$pvalue)
summary(deseq_results$padj)








  
  
  
  
  
  
  
  
  
  


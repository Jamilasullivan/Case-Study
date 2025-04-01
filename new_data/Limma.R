## PACKAGE INSTALLATION AND LOADING ############################################

## only run the following two commands on first run of the script

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("limma")

## loading packages

library(edgeR)
library(limma)

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

## converting counts data for use with Limma 

#counts <- voom(counts, design, plot = TRUE)

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

results <- DGEList(counts) # identifies differentially expressed genes
head(results)
results
#View(results)

## PREPROCESSING ###############################################################

results <- calcNormFactors(results) # calculates normalisation factors 
results$samples

## filtering results

cutoff <- 10 
drop <- which(apply(cpm(results),1,max) < cutoff)
filtered <- results[-drop,]
dim(filtered)

#set sample names

sample_names <- colnames(counts)
sample_names

## Voom transformation to make data suitable for Limma

design <- model.matrix(~ metadata$Condition) # creates the design of the analysis
filtered <- voom(filtered, design, plot = TRUE) # log transformation and estimates mean-variance trends

## fitting the linear model

fit <- lmFit(filtered, design) # fits linear mode for each gene

## emperical Baysian smoothing

fit <- eBayes(fit)

## EXTRACTING DIFFERENTIALLY EXPRESSED GENES #####################

results <- topTable(fit, coef = 2, adjust.method = "BH", number = Inf) # extracts a table of the top ranked genes from a linear model fit

significant_genes <- results[results$adj.P.Val < 0.05, ]

write.csv(results, "limma_differential_expression_results.csv", row.names = TRUE)

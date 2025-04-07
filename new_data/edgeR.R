## PACKAGE INSTALLATION AND LOADING ############################################

## only run the following two commands on first run of the script

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

## loading packages

library(edgeR)
library(ggplot2)

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

## ANALYSIS ####################################################################

## create DGEList object

results <- DGEList(counts) # identifies differentially expressed genes
head(results)
results
#View(results)

## PREPROCESSING ###############################################################

results <- calcNormFactors(results) # calculates normalisation factors 
results$samples

results <- estimateDisp(results)

## ANALYSIS ####################################################################

et <- exactTest(results)
results <- topTags(et, n=Inf)
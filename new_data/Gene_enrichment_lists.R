## LOADING PACKAGES ############################################################

library(dplyr)

## READING IN DATA #############################################################

## set working directory

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data")

## making object for testing

limma_results <- read.csv("limma_differential_expression_results.csv", row.names = 1) # read in the limma results 
summary(limma_results)
#view(limma_results)

# organise the data by adjusted p-value

limma_results_ordered <- limma_results[order(limma_results$adj.P.Val), ] 
#view(limma_results_ordered)

# filter by log2fold change

limma_results_ordered <- limma_results_ordered %>% filter(abs(limma_results_ordered$logFC) > 1) 

# filter the results by an adjusted p value of 0.05

filtered_0.05 <- limma_results_ordered %>% filter(limma_results_ordered$adj.P.Val < 0.05) 

filtered_0.01 <- limma_results_ordered %>% filter(limma_results_ordered$adj.P.Val < 0.01) 

filtered_0.005 <- limma_results_ordered %>% filter(limma_results_ordered$adj.P.Val < 0.005) 

filtered_0.001 <- limma_results_ordered %>% filter(limma_results_ordered$adj.P.Val < 0.001) 

## extract gene lists for Gene Profiler

genes_0.05 <- rownames(filtered_0.05)
genes_0.01 <- rownames(filtered_0.01)
genes_0.005 <- rownames(filtered_0.005)
genes_0.001 <- rownames(filtered_0.001)

## save gene lists as text files

writeLines(genes_0.05, "genes_0.05.txt")
writeLines(genes_0.01, "genes_0.01.txt")
writeLines(genes_0.005, "genes_0.005.txt")
writeLines(genes_0.001, "genes_0.001.txt")

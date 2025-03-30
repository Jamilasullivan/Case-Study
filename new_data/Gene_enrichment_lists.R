## READING IN DATA #############################################################

## set working directory

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data")

## making object for testing

deseq_results <- read.csv("all_deseq_results_siglecf.csv", row.names = 1)
summary(deseq_results)
#view(deseq_results)

# organise the data by adjusted p-value

deseq_results_ordered <- deseq_results[order(deseq_results$padj), ] 
#view(deseq_results_ordered)

# filter by log2fold change

deseq_results_ordered <- deseq_results_ordered %>% filter(abs(deseq_results_ordered$log2FoldChange) > 1) 

# filter the results by an adjusted p value of 0.05

filtered_0.05 <- deseq_results_ordered %>% filter(deseq_results_ordered$padj < 0.05) 

filtered_0.01 <- deseq_results_ordered %>% filter(deseq_results_ordered$padj < 0.01) 

filtered_0.005 <- deseq_results_ordered %>% filter(deseq_results_ordered$padj < 0.005) 

filtered_0.001 <- deseq_results_ordered %>% filter(deseq_results_ordered$padj < 0.001) 

## extract gene lists for Gene Profiler

genes_0.05 <- rownames(filtered_0.05)
genes_0.01 <- rownames(filtered_0.01)
genes_0.005 <- rownames(filtered_0.005)
genes_0.001 <- rownames(filtered_0.001)

class(genes_0.05)

## save gene lists as text files

writeLines(genes_0.05, "genes_0.05.txt")
writeLines(genes_0.01, "genes_0.01.txt")
writeLines(genes_0.005, "genes_0.005.txt")
writeLines(genes_0.001, "genes_0.001.txt")

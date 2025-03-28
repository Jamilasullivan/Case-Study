## installing and loading necessary packages ###################################

# only run the following two 'if' commands the first time you run the script

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}

# intalling necessary programmes

BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("pathview")
BiocManager::install("ReactomePA")
BiocManager::install("ggplot2")
install.packages("tibble")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(cowplot)
library(ggplot2)
library(pathview)
library(ReactomePA)
library(dplyr)
library(tibble)

## set working directory

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT101. Biocomputing and Big Data Handling/Mini Project Bioinformatics/Data sorting/Data_sorting/data_analysis")

## making object for testing

deseq_results <- read.csv("all_deseq_results.csv", row.names = 1)
summary(deseq_results)
deseq_filtered <- deseq_results %>% filter(deseq_results$padj < 0.05)
deseq_filtered <- deseq_filtered %>% filter(abs(deseq_filtered$log2FoldChange) > 0.05)
deseq_filtered <- deseq_filtered[order(deseq_filtered$padj),]
deseq_filtered
dim(deseq_filtered)
summary(deseq_filtered)

genes_to_test <- rownames(deseq_filtered)
print(genes_to_test)  

#### GO ENRICHMENT ANALYSIS ###########################################

# biological processes - 6
ego_BP <- enrichGO(gene = genes_to_test,
                   universe = names(genes_to_test),
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP")
head(ego_BP)

# cellular component - 6
ego_CC <- enrichGO(gene = genes_to_test,
                   universe = names(genes_to_test),
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "CC")
head(ego_CC)

# molecular function - 6
ego_MF <- enrichGO(gene = genes_to_test,
                   universe = names(genes_to_test),
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "MF")
head(ego_MF)

# plot results
plot_egoBP <- plot(barplot(ego_BP, showCategory = 20, font.size = 5))
plot_egoCC <- plot(barplot(ego_CC, showCategory = 20, font.size = 5))
plot_egoMF <- plot(barplot(ego_MF, showCategory = 20, font.size = 5))

combined_ego_plot <- plot_grid(plot_egoCC, plot_egoBP, plot_egoMF, ncol = 1)
print(combined_ego_plot)

## ggplot of less data 

# Extract and filter the top 5 terms based on p.adjust for each ontology
data_BP <- as.data.frame(ego_BP)
data_BP$Ontology <- "BP"
top_BP <- data_BP[order(data_BP$p.adjust), ][1:7, ]  # Top 5 by adjusted p-value

data_CC <- as.data.frame(ego_CC)
data_CC$Ontology <- "CC"
top_CC <- data_CC[order(data_CC$p.adjust), ][1:7, ]

data_MF <- as.data.frame(ego_MF)
data_MF$Ontology <- "MF"
top_MF <- data_MF[order(data_MF$p.adjust), ][1:7, ]

# Combine top terms from all categories
combined_data <- rbind(top_BP, top_CC, top_MF)

# Select relevant columns (Description, Count, Ontology)
plot_data <- combined_data[, c("Description", "Count", "Ontology")]

plot_data$Ontology <- factor(plot_data$Ontology, levels = c("MF", "BP", "CC"))

ggplot(plot_data, aes(x = reorder(Description, Count), y = Count, fill = Ontology)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "GO Term Enrichment",x = NULL, y = "Gene Count") +
  theme_minimal_grid() +
  scale_fill_manual(values = c("BP" = "#202495", "MF" = "#af2338", "CC" = "#149413")) +
  coord_flip()+
  theme()

#### GENE SET ENRICHMENT ANALYSIS ############################

#data organisation
deseq_filtered2 <- deseq_filtered[order(-deseq_filtered$log2FoldChange),]
deseq_filtered2
summary(deseq_filtered2)

## extract stat column

gene_list <- deseq_filtered2$log2FoldChange
names(gene_list) <- rownames(deseq_filtered2)
as.data.frame(gene_list)
gene_list

# gene set enrichment analysis

# BP - 6
gse_BP <- gseGO(gene_list,
                ont = "BP",
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db")
head(gse_BP)

# CC - 6 
gse_CC <- gseGO(gene_list,
                ont = "CC",
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db")
head(gse_CC)

## MF - 6
gse_MF <- gseGO(gene_list,
                ont = "MF",
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db")
head(gse_MF)

## plotting ##

#GSAE plots

# These may vary slightly for each different run of the code dues to adjusted p-values 

gseaplot(gse_BP, geneSetID = 1) # mostly downregulated # BP - hormone metabolic process   
gseaplot(gse_BP, geneSetID = 2) # mostly downregulated # BP - organic acid biosynthetic process
gseaplot(gse_BP, geneSetID = 3) # mostly downregulated # BP - carboxylic acid biosynthetic process 
gseaplot(gse_BP, geneSetID = 4) # mostly downregulated # BP - cellular lipid catabolic process  
gseaplot(gse_BP, geneSetID = 5) # mostly downregulated # BP - digestion
gseaplot(gse_BP, geneSetID = 6) # mostly downregulated # BP - organic acid catabolic process

gseaplot(gse_CC, geneSetID = 1) # mostly downregulated # CC - cluster of actin-based cell projections
gseaplot(gse_CC, geneSetID = 2) # mostly downregulated # CC - brush border
gseaplot(gse_CC, geneSetID = 3) # mostly downregulated # CC - endoplasmic reticulum lumen
gseaplot(gse_CC, geneSetID = 4) # mostly downregulated # CC - peroxisome
gseaplot(gse_CC, geneSetID = 5) # mostly downregulated # CC - microbody
gseaplot(gse_CC, geneSetID = 6) # mostly downregulated # CC - primary lysosome 

gseaplot(gse_MF, geneSetID = 1) # mostly downregulated # MF - endopeptidase activity
gseaplot(gse_MF, geneSetID = 2) # mostly downregulated # MF - carboxylic acid binding
gseaplot(gse_MF, geneSetID = 3) # mostly downregulated # MF - monooxygenase activity
gseaplot(gse_MF, geneSetID = 4) # mostly downregulated # MF - monocarboxylic acid binding
gseaplot(gse_MF, geneSetID = 5) # mostly downregulated # MF - glycosaminoglycan binding 
gseaplot(gse_MF, geneSetID = 6) # mostly downregulated # MF - metallopeptidase activity

gene_list[1]
gene_list[3825]

# Plotting with ggplot2


# bar plot# bar plotrank()

data_BP <- as.data.frame(gse_BP)
data_BP$Ontology <- "BP"
top_BP <- data_BP[order(data_BP$p.adjust), ][1:6, ]
top_BP

data_CC <- as.data.frame(gse_CC)
data_CC$Ontology <- "CC"
top_CC <- data_CC[order(data_CC$p.adjust), ][1:6, ]
top_CC

data_MF <- as.data.frame(gse_MF)
data_MF$Ontology <- "MF"
top_MF <- data_MF[order(data_MF$p.adjust), ][1:6, ]
top_MF

# Combine top terms from all categories
combined_data2 <- rbind(top_CC, top_MF, top_BP)
combined_data2

# Select relevant columns (Description, setSize, Ontology)
plot_data2 <- combined_data2[, c("Description", "setSize", "Ontology")]

plot_data2

ggplot(plot_data2, aes(x = reorder(Description, setSize), y = setSize, fill = Ontology)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "GO Term Enrichment",x = NULL, y = "Gene Count") +
  theme_minimal_grid() +
  scale_fill_manual(values = c("BP" = "#202495", "MF" = "#af2338", "CC" = "#149413")) +
  coord_flip()+
  theme()

#### GO ENRICHED PATHWAYS ####################################################

gene_id <-bitr(rownames(deseq_filtered), 
               fromType = "SYMBOL", 
               toType = "ENTREZID", 
               OrgDb= "org.Hs.eg.db")

gene_id <- gene_id$ENTREZID

KEGG <- enrichKEGG(gene = gene_id,
                   organism = "hsa")

KEGG_plot <- KEGG[, c("Description", "Count", "p.adjust")]

head(KEGG)
nrow(KEGG)

ggplot(KEGG_plot, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "GO Term Enrichment",x = NULL, y = "Gene Count") +
  theme_minimal_grid() +
  scale_fill_gradient(low = "#003efc", high = "#C1CFFB") +
  coord_flip()+
  theme()

## GENE SET PATHWAY ENRICHMENT #################################################

# Create the gene list: log2 fold changes as values, Entrez IDs as names
summary(deseq_filtered)
gse_genes <- deseq_filtered$log2FoldChange

# Ensure the gene list has Entrez gene IDs as names
names(gse_genes) <- rownames(deseq_filtered)
gse_genes

# Sort the gene list in decreasing order (important for GSEA)
gse_genes <- sort(gse_genes, decreasing = TRUE)
summary(gse_genes)
gse_genes

gene_id <-bitr(rownames(deseq_filtered), 
               fromType = "SYMBOL", 
               toType = "ENTREZID",  
               OrgDb= "org.Hs.eg.db") # 5.9% of genes lost leaving 6343 genes

deseq_filtered

deseq_filtered <- rownames_to_column(deseq_filtered, "SYMBOL")  # Convert rownames to a column
deseq_filtered <- left_join(deseq_filtered, gene_id, by = "SYMBOL")  # Join to match ENTREZID

# Remove rows with missing ENTREZID
deseq_filtered <- deseq_filtered %>% filter(!is.na(ENTREZID))

# Create the named vector for GSEA
gse_genes <- deseq_filtered$log2FoldChange
names(gse_genes) <- deseq_filtered$ENTREZID
gse_genes <- sort(gse_genes, decreasing = TRUE)
gse_genes
summary(gse_genes)

gse_KEGG <- gseKEGG(geneList = gse_genes,
                    organism = "hsa")

head(gse_KEGG)
nrow(gse_KEGG)
gse_KEGGplot <- as.data.frame(gse_KEGG)

ggplot(gse_KEGGplot, aes(x = reorder(Description, setSize), y = setSize, fill = p.adjust)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "GO Term Enrichment",x = NULL, y = "Gene Count") +
  theme_minimal_grid() +
  scale_fill_gradient(low = "#003efc", high = "#C1CFFB") +
  coord_flip()+
  theme()

head(gse_KEGG, n=41)

gseaplot(gse_KEGG, geneSetID = 1)
gseaplot(gse_KEGG, geneSetID = 2)
gseaplot(gse_KEGG, geneSetID = 3)
gseaplot(gse_KEGG, geneSetID = 4)
gseaplot(gse_KEGG, geneSetID = 5)
gseaplot(gse_KEGG, geneSetID = 6)
gseaplot(gse_KEGG, geneSetID = 7)
gseaplot(gse_KEGG, geneSetID = 8)
gseaplot(gse_KEGG, geneSetID = 9)
gseaplot(gse_KEGG, geneSetID = 10)
gseaplot(gse_KEGG, geneSetID = 11)
gseaplot(gse_KEGG, geneSetID = 12)
gseaplot(gse_KEGG, geneSetID = 13)
gseaplot(gse_KEGG, geneSetID = 14)
gseaplot(gse_KEGG, geneSetID = 15)
gseaplot(gse_KEGG, geneSetID = 16)
gseaplot(gse_KEGG, geneSetID = 17)
gseaplot(gse_KEGG, geneSetID = 18)
gseaplot(gse_KEGG, geneSetID = 19)
gseaplot(gse_KEGG, geneSetID = 20)
gseaplot(gse_KEGG, geneSetID = 21)
gseaplot(gse_KEGG, geneSetID = 22)
gseaplot(gse_KEGG, geneSetID = 23)
gseaplot(gse_KEGG, geneSetID = 24)
gseaplot(gse_KEGG, geneSetID = 26)
gseaplot(gse_KEGG, geneSetID = 27)
gseaplot(gse_KEGG, geneSetID = 28)
gseaplot(gse_KEGG, geneSetID = 29)
gseaplot(gse_KEGG, geneSetID = 30)
gseaplot(gse_KEGG, geneSetID = 31)
gseaplot(gse_KEGG, geneSetID = 32)
gseaplot(gse_KEGG, geneSetID = 33)
gseaplot(gse_KEGG, geneSetID = 34)
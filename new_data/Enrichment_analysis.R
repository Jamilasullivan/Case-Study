## installing and loading necessary packages ###################################

# only run the following two 'if' commands the first time you run the script

#if (!requireNamespace("BiocManager", quietly = TRUE)){
#  install.packages("BiocManager")
#}

#if (!requireNamespace("cowplot", quietly = TRUE)) {
#  install.packages("cowplot")
#}

# intalling necessary programmes

#BiocManager::install("clusterProfiler")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("enrichplot")
#BiocManager::install("pathview")
#BiocManager::install("ReactomePA")
#BiocManager::install("ggplot2")
#BiocManager::install("ggarchery")
#install.packages("tibble")
#BiocManager::install("enrichplot")
# loading packages

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(cowplot)
library(ggplot2)
library(pathview)
library(ReactomePA)
library(dplyr)
library(tibble)
library(enrichplot)
library(GOplot)
library(enrichplot)

## set working directory

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data")

## making object for testing

limma_results <- read.csv("limma_differential_expression_results.csv", row.names = 1) # read in results file
limma_filtered <- limma_results %>% filter(limma_results$adj.P.Val < 0.05) # filter to results with an adjusted p value below 0.05
limma_filtered <- limma_filtered %>% filter(abs(limma_filtered$logFC) > 1) # filter to results with a log 2 fold change of +/- 1.
limma_filtered <- limma_filtered[order(limma_filtered$adj.P.Val),] # order the filter results by adjusted p value
limma_filtered
dim(limma_filtered)
summary(limma_filtered)

genes_to_test <- rownames(limma_filtered)
print(genes_to_test)  

################################################################################
########################## GO ENRICHMENT ANALYSIS ##############################
################################################################################

## enrichGO looks for overrepresented terms in significant genes

###### biological processes ######
ego_BP <- enrichGO(gene = genes_to_test,
                   universe = names(genes_to_test),
                   keyType = "SYMBOL",
                   OrgDb = org.Mm.eg.db,
                   ont = "BP")
head(ego_BP)
#view(ego_BP)

##### cellular component #####
ego_CC <- enrichGO(gene = genes_to_test,
                   universe = names(genes_to_test),
                   keyType = "SYMBOL",
                   OrgDb = org.Mm.eg.db,
                   ont = "CC")
head(ego_CC)
#view(ego_CC)

##### molecular function #####
ego_MF <- enrichGO(gene = genes_to_test,
                   universe = names(genes_to_test),
                   keyType = "SYMBOL",
                   OrgDb = org.Mm.eg.db,
                   ont = "MF")
head(ego_MF)
#view(ego_MF)

ego_all <- enrichGO(gene = genes_to_test,
                    universe = names(genes_to_test),
                    keyType = "SYMBOL",
                    OrgDb = org.Mm.eg.db)

head(ego_all)
#view(ego_all)

############################# plot results #####################################
plot_egoBP <- plot(barplot(ego_BP, showCategory = 20, font.size = 5))
plot_egoCC <- plot(barplot(ego_CC, showCategory = 20, font.size = 5))
plot_egoMF <- plot(barplot(ego_MF, showCategory = 20, font.size = 5))

combined_ego_plot <- plot_grid(plot_egoCC, plot_egoBP, plot_egoMF, ncol = 1)
print(combined_ego_plot)

## barplot

barplot(ego_all, showCategory=20)

## dotplot

dotplot(ego_all, showCategory=20)

## GOplot

goplot(ego_all)

## enrichment map 

cnetplot(ego_all, showCategory=20)

## GO netwrok visualisation

cnetplot(ego_all, categorySize="pvalue", foldChange=geneList)

## circular plot

circosplot(ego_all)

# Install enrichplot if not already installed
install.packages("enrichplot")

# Load the package
library(enrichplot)

install.packages("circlize")
library(circlize)

packageVersion("enrichplot")
# To update the package from GitHub if necessary (requires devtools)
devtools::install_github("YuLab-SMU/enrichplot")
library(enrichplot)

# Update enrichplot and clusterProfiler packages
install.packages("enrichplot")
install.packages("clusterProfiler")


# Load the libraries
library(enrichplot)
library(clusterProfiler)

remove.packages("enrichplot")
# Install devtools if it's not already installed
install.packages("devtools")

# Load devtools
library(devtools)

# Install the latest version of enrichplot from GitHub
devtools::install_github("YuLab-SMU/enrichplot")

# Check the version of enrichplot
packageVersion("enrichplot")

# Install remotes if you don't have it already
install.packages("remotes")

# Install enrichplot from GitHub using remotes
remotes::install_github("YuLab-SMU/enrichplot", force = TRUE)

 
##### ggplot of less data #####

# Extract and filter the top 5 terms based on p.adjust for each ontology
data_BP <- as.data.frame(ego_BP)
data_BP$Ontology <- "BP"
top_BP <- data_BP[order(data_BP$p.adjust), ][1:7, ]  # Top 7 by adjusted p-value

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

################################################################################
####################### GENE SET ENRICHMENT ANALYSIS ###########################
################################################################################

## rank-based enrichment analysis on all genes involved

limma_results2 <- read.csv("limma_differential_expression_results.csv", row.names = 1) # read in results file

#data organisation
limma_filtered2 <- limma_results2[order(-limma_results2$logFC),] # ordering data by descending log2foldchange
limma_filtered2
summary(limma_filtered2)
#view(limma_filtered2)

##### extract stat column #####

gene_list <- limma_filtered2$logFC
names(gene_list) <- rownames(limma_filtered2)
as.data.frame(gene_list)
gene_list

######################### gene set enrichment analysis #########################

##### BP #####
gse_BP <- gseGO(gene_list,
                ont = "BP",
                keyType = "SYMBOL",
                OrgDb = "org.Mm.eg.db")
head(gse_BP)
#view(gse_BP)

##### CC #####
gse_CC <- gseGO(gene_list,
                ont = "CC",
                keyType = "SYMBOL",
                OrgDb = "org.Mm.eg.db")
head(gse_CC)
#view(gse_CC)

##### MF #####
gse_MF <- gseGO(gene_list,
                ont = "MF",
                keyType = "SYMBOL",
                OrgDb = "org.Mm.eg.db")
head(gse_MF)
#view(gse_MF)

############################## PLOTTING ########################################

############################## GSAE plots ######################################

# These may vary slightly for each different run of the code dues to adjusted p-values 

gseaplot(gse_BP, geneSetID = 1)   
gseaplot(gse_BP, geneSetID = 2) 
gseaplot(gse_BP, geneSetID = 3)
gseaplot(gse_BP, geneSetID = 4)
gseaplot(gse_BP, geneSetID = 5)
gseaplot(gse_BP, geneSetID = 6)

gseaplot(gse_CC, geneSetID = 1) 
gseaplot(gse_CC, geneSetID = 2)
gseaplot(gse_CC, geneSetID = 3)
gseaplot(gse_CC, geneSetID = 4) 
gseaplot(gse_CC, geneSetID = 5)
gseaplot(gse_CC, geneSetID = 6)

gseaplot(gse_MF, geneSetID = 1) 
gseaplot(gse_MF, geneSetID = 2)
gseaplot(gse_MF, geneSetID = 3)
gseaplot(gse_MF, geneSetID = 4)
gseaplot(gse_MF, geneSetID = 5)
gseaplot(gse_MF, geneSetID = 6)

gene_list[1]
gene_list[8000] # alter the number to see how long the list is

######################### Plotting with ggplot2 ################################

################################ bar plot ######################################

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

########################## GO ENRICHED PATHWAYS ################################

gene_id <-bitr(rownames(limma_filtered),
               fromType = "SYMBOL", 
               toType = "ENTREZID", 
               OrgDb= "org.Mm.eg.db")

gene_id <- gene_id$ENTREZID

KEGG <- enrichKEGG(gene = gene_id,
                   organism = "mmu")

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

######################## GENE SET PATHWAY ENRICHMENT ###########################

# Create the gene list: log2 fold changes as values, Entrez IDs as names
summary(limma_filtered)
gse_genes <- limma_filtered$log2FoldChange

# Ensure the gene list has Entrez gene IDs as names
names(gse_genes) <- rownames(limma_filtered)
gse_genes

# Sort the gene list in decreasing order (important for GSEA)
gse_genes <- sort(gse_genes, decreasing = TRUE)
summary(gse_genes)
gse_genes

gene_id <-bitr(rownames(limma_filtered),
               fromType = "SYMBOL", 
               toType = "ENTREZID",  
               OrgDb= "org.Mm.eg.db")

limma_filtered

limma_filtered <- rownames_to_column(limma_filtered, "SYMBOL")  # Convert rownames to a column
limma_filtered <- left_join(limma_filtered, gene_id, by = "SYMBOL")  # Join to match ENTREZID

# Remove rows with missing ENTREZID
limma_filtered <- limma_filtered %>% filter(!is.na(ENTREZID))

# Create the named vector for GSEA
gse_genes <- limma_filtered$log2FoldChange
names(gse_genes) <- limma_filtered$ENTREZID
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
                    organism = "mmu")

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


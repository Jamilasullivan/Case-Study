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

## loading packages

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
library(stringr)
library(ggupset)

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
#view(genes_to_test)

################################################################################
########################## GO ENRICHMENT ANALYSIS ##############################
################################################################################

## enrichGO looks for overrepresented terms in significant genes.
## enrich uses a filtered list of significant genes
## Good for DEG analysis

################ biological processes (overrepresented GO terms) ###############

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

##### all enrichment terms together ######
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
print(combined_ego_plot) # combines the top 20 results of all ontology terms into one plot object

##### bar plot #####

barplot(ego_all, showCategory=20)

##### dotplot #####

dotplot(ego_all, showCategory=20)

##### GOplot #####

goplot(ego_all)

##### GO network visualisation #####

## the following lines produce the gene_list again to be able to add colour

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

#view(ego_all)

cnet_ego_all <- cnetplot(ego_all, 
              showCategory = 10, 
              circular = TRUE, # Makes the plot radial and cleaner
              size_item = 0.8,
              color_category = "#E9D37A",
              color_edge = "grey",
              size_edge = 0.4,
              node_label = "category",
              categorySize="p.adjust", 
              foldChange=gene_list) # Colors edges by category

cnet_ego_all + 
  labs(title = "GO Term – Gene Network", subtitle = "Top 10 Enriched GO Terms", size = "Number of\ngenes involved") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
  )

#### or 

ego_all2 <- simplify(ego_all)
cnetplot(ego_all2, foldChange=gene_list) # use this to see all gene names but it can get messy.

##### upset plot #####

upsetplot(ego_all)
 
######################## ggplot of less data ###################################

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

######################## Pathway enrichment analysis ###########################

gene_id <-bitr(rownames(limma_filtered),
               fromType = "SYMBOL", 
               toType = "ENTREZID", 
               OrgDb= "org.Mm.eg.db")

gene_id <- gene_id$ENTREZID

KEGG_all <- enrichKEGG(gene = gene_id,
                       organism = "mmu")

KEGG_plot <- KEGG_all[, c("Description", "Count", "p.adjust")]

head(KEGG_all)
nrow(KEGG_all)

ggplot(KEGG_plot, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "GO Term Enrichment",x = NULL, y = "Gene Count") +
  theme_minimal_grid() +
  scale_fill_gradient(low = "#003efc", high = "#C1CFFB") +
  coord_flip()+
  theme()

##### bar plot #####

barplot(KEGG_all, showCategory=20)

##### dotplot #####

dotplot(KEGG_all, showCategory=20)

##### GOplot #####

kegg_res_readable <- setReadable(KEGG_all, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
goplot(kegg_res_readable) # not working

##### GO network visualisation #####

## the following lines produce the gene_list again to be able to add colour

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

#view(KEGG_all)

cnet_KEGG_all <- cnetplot(KEGG_all, 
                         showCategory = 10, 
                         circular = TRUE, # Makes the plot radial and cleaner
                         size_item = 0.8,
                         color_category = "#E9D37A",
                         color_edge = "grey",
                         size_edge = 0.4,
                         node_label = "category",
                         categorySize="p.adjust", 
                         foldChange=gene_list) # Colors edges by category

cnet_KEGG_all + 
  labs(title = "GO Term – Gene Network", subtitle = "Top 10 Enriched KEGG Pathways", size = "Number of\ngenes involved") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
  )

##### upset plot #####

upsetplot(KEGG_all)

################################################################################
################## GO GENE SET ENRICHMENT ANALYSIS (GSEA) ######################
################################################################################

## GSEA takes ranked list of all genes
## gse commands are good for overall trends

## rank-based enrichment analysis on all genes involved

limma_results2 <- read.csv("limma_differential_expression_results.csv", row.names = 1) # read in results file

#data organisation
limma_filtered2 <- limma_results2[order(-limma_results2$logFC),] # ordering data by descending log2foldchange
limma_filtered2
summary(limma_filtered2)
view(limma_filtered2)

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

gse_all <- gseGO(gene_list,
                 keyType = "SYMBOL",
                 OrgDb = "org.Mm.eg.db")
head(gse_all)
#view(gse_all)

############################## PLOTTING ########################################

######################## GSAE plots for specific terms #########################

# These may vary slightly for each different run of the code dues to adjusted p-values. Change the number following 'geneSetID' to the number of the term of interest in the object. 

gseaplot2(gse_BP, geneSetID = 1, title = gse_all@result$Description[1])  

gseaplot2(gse_CC, geneSetID = 1, title = gse_all@result$Description[1])

gseaplot2(gse_MF, geneSetID = 1, title = gse_all@result$Description[1])

gseaplot2(gse_all, geneSetID = 1, title = gse_all@result$Description[1])

gene_list[1]
gene_list[8000] # alter the number to see how long the list is

######################## GSAE plots for the full set ###########################

##### ridge plot #####

ridge_all <- ridgeplot(gse_all, showCategory = 20)

ridge_all + 
  scale_fill_gradient(low = "#FF4F4F", high = "#570000") +  # Blue gradient
  labs(title = "GO Enrichment Ridge Plot",
       x = "Enrichment Score",
       y = "GO Term") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_blank(),  # Clean up background
    panel.grid.minor = element_blank()
  )

##### dot plot #####

dot_all <- dotplot(gse_all, showCategory = 20)$data

ggplot(dot_all, aes(x = NES, 
                     y = reorder(Description, NES), 
                     size = setSize, 
                     color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#56B1F7", high = "#132B43", name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Set Size") +
  labs(
    title = "Top Enriched GO Terms (GSEA)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_blank()
  )

##### heatmap #####

heatplot_gse_all <- heatplot(gse_all, showCategory = 10)

heatplot_gse_all + 
  labs(title = "Top 10 Enriched GO Terms")  # Title for the plot

##### enrichment map #####

emapplot(gse_all, showCategory = 30) # doesn't work
emapplot(pairwise_termsim(gse_all), showCategory = 30) # works

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

################################################################################
######################## GENE SET PATHWAY ENRICHMENT ###########################
################################################################################



# Create the gene list: log2 fold changes as values, Entrez IDs as names
summary(limma_filtered)
#view(limma_filtered)
gse_genes <- limma_filtered$logFC

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
gse_genes <- limma_filtered$logFC
names(gse_genes) <- limma_filtered$ENTREZID
limma_filtered
view(limma_filtered)

limma_filtered <- as.data.frame(limma_filtered)

limma_filtered <- rownames_to_column(limma_filtered, "SYMBOL")  # Convert rownames to a column
limma_filtered <- left_join(limma_filtered, gene_id, by = "SYMBOL")  # Join to match ENTREZID

# Remove rows with missing ENTREZID
limma_filtered <- limma_filtered %>% filter(!is.na(ENTREZID))

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


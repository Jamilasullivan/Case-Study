## LOAD PACKAGES ##############################################################

#BiocManager::install("Gviz")
#BiocManager::install("dzhang32/dasper")

library(ggplot2)
library(dplyr)
library(ggforce)  # for arcs
library(tidyr)
library(gridExtra)

## SET WORKING DIRECTORY ######################################################

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data")

###############################################################################
##################### READ IN AS DATA FROM RMATS ##############################
###############################################################################

## All data here is from the AIRP-CTRL 
## IncLevelDifference goes from condition 1 to condition 2 (<0 - higher exon skipping in condition 1. 0 = minimal difference in splicing). 

##### skipped exon results #####
se_rmats <- read.table("SE.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

##### mutually exclusive exon results #####
mxe_rmats <- read.table("MXE.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

##### alternative 3' splice site results #####
a3ss_rmats <- read.table("A3SS.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

##### alternative 5' splice site results #####
a5ss_rmats <- read.table("A5SS.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

##### retained intron results #####
ri_rmats <- read.table("RI.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## PREVIEW OF ALL RESULTS FILES ###############################################

head(se_rmats)
#View(se_rmats)
head(mxe_rmats)
#view(mxe_rmats)
head(a3ss_rmats)
#view(a3ss_rmats)
head(a5ss_rmats)
#view(a5ss_rmats)
head(ri_rmats)
#view(ri_rmats)

## FILTERING RESULTS ##########################################################

# in rMATS FDR is the adjusted p value
# IncLevelDifference is difference in exon inclusion between case and control (0 to 1)

## filter by FDR < 0.05 and IncLevelDifference > 0.1
sig_se <- subset(se_rmats, FDR < 0.05 & abs(IncLevelDifference) > 0.1)
sig_mxe <- subset(mxe_rmats, FDR < 0.05 & abs(IncLevelDifference) > 0.1)
sig_a3ss <- subset(a3ss_rmats, FDR < 0.05 & abs(IncLevelDifference) > 0.1)
sig_a5ss <- subset(a5ss_rmats, FDR < 0.05 & abs(IncLevelDifference) > 0.1)
sig_ri <- subset(ri_rmats, FDR < 0.05 & abs(IncLevelDifference) > 0.1)

## filter by FDR < 0.1 and IncLevelDifference > 0.05
sig2_se <- subset(se_rmats, FDR < 0.1 & abs(IncLevelDifference) > 0.05)
sig2_mxe <- subset(mxe_rmats, FDR < 0.1 & abs(IncLevelDifference) > 0.05)
sig2_a3ss <- subset(a3ss_rmats, FDR < 0.1 & abs(IncLevelDifference) > 0.05)
sig2_a5ss <- subset(a5ss_rmats, FDR < 0.1 & abs(IncLevelDifference) > 0.05)
sig2_ri <- subset(ri_rmats, FDR < 0.1 & abs(IncLevelDifference) > 0.05)

###############################################################################
############################# PLOTTING ########################################
###############################################################################

############################ HISTOGRAMS #######################################

##### IncLevelDifference histogram for skipped exons #####
ggplot(sig2_se, aes(x = IncLevelDifference)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "ΔPSI Distribution for Significant SE Events",
    x = "Inclusion Level Difference (ΔPSI)",
    y = "Number of Events"
  ) # more exon inclusion was seen in air pollution because more values were close to 1

##### IncLevelDifference histogram for alternative 3' splice sites #####
ggplot(sig2_a3ss, aes(x = IncLevelDifference)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "ΔPSI Distribution for Significant SE Events",
    x = "Inclusion Level Difference (ΔPSI)",
    y = "Number of Events"
  ) + 
  xlim(-1,1) # more exon inclusion was seen in air pollution because more values were close to 1

############################## VOLCANO PLOTS ###################################

## volcano plots don't work massively well for this data

se_rmats$logFDR <- -log10(se_rmats$FDR)
se_rmats$IncLevelDifference <- as.numeric(se_rmats$IncLevelDifference)

ggplot(se_rmats, aes(x = IncLevelDifference, y = logFDR)) +
  geom_point(alpha = 0.5, color = "darkred") +
  theme_minimal() +
  labs(
    title = "Volcano-Like Plot of SE Events",
    x = "ΔPSI (Inclusion Level Difference)",
    y = "-log10(FDR)"
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

########################################################################################################### SASHIMI-STYLE PLOTS ################################
################################################################################

######################## SKIPPED EXON EVENTS ###################################

se_rmats_ordered <- se_rmats %>%
  arrange(FDR)
#View(se_rmats_ordered)
head(se_rmats_ordered)

# Subset the top rows of the ordered se_rmats data
se_rmats_subset <- se_rmats_ordered[1:2, ] # Adjust for actual data frame

# List to store individual plots
plot_list <- list()

# Loop through each row (event) and generate individual plots
for (i in 1:nrow(se_rmats_subset)) {
  event <- se_rmats_subset[i, ] # Extract each event
  
  # Coordinates for the current event
  chrom <- event$chr
  strand <- event$strand
  exonStart <- event$exonStart_0base
  exonEnd <- event$exonEnd
  upstreamExonEnd <- event$upstreamES
  upstreamExonStart <- event$upstreamEE
  downstreamExonStart <- event$downstreamES
  downstreamExonEnd <- event$downstreamEE
  
  # Junction counts 
  IJC1 <- event$IJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  SJC1 <- event$SJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  
  # Sum junction counts 
  IJC_sum <- sum(IJC1)
  SJC_sum <- sum(SJC1)
  
  # Exon boxes
  exons <- data.frame(
    xmin = c(upstreamExonStart, exonStart, downstreamExonStart),
    xmax = c(upstreamExonEnd, exonEnd, downstreamExonEnd),
    ymin = -0.1,
    ymax = 0.1,
    label = c("Upstream Exon", "Skipped Exon", "Downstream Exon")
  )
  
  # Define arcs
  junctions <- data.frame(
    x = c(upstreamExonEnd, upstreamExonEnd),
    xend = c(exonStart, downstreamExonStart),
    y = 0.1,
    yend = 0.1,
    count = c(IJC_sum, SJC_sum),
    label = c("Inclusion", "Skipping")
  )
  
  # Plot for the current event
  plot <- ggplot() +
    # Exon boxes with custom fill color
    geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
              fill = "#0473C7", color = "black", size = 0.5) +
    
    # Arcs for splice junctions with customization (color and line thickness)
    geom_curve(data = junctions, 
               aes(x = x, xend = xend, y = y, yend = yend), 
               curvature = 0.3, 
               linewidth = 1, 
               color = "#B20101") +
    
    # Junction count labels (position, size, color)
    geom_text(data = junctions, 
              aes(x = (x + xend) / 2, y = y + 0.15, label = count), 
              size = 5, color = "black") +
    
    # Title and axis labels with gene name
    labs(title = paste("Sashimi-style Plot (SE Event)\nGene:", event$geneSymbol),
         subtitle = paste("Chromosome:", chrom),
         x = "Genomic Position", y = "") +
    
    # Customizing the plot theme
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(), 
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      legend.position = "none"  
    )
  
  # Store plot in the list
  plot_list[[i]] <- plot
}

grid.arrange(grobs = plot_list, ncol = 2)

###################### MUTUALLY EXCLUSIVE EXON EVENTS ##########################

mxe_rmats_ordered <- mxe_rmats %>%
  arrange(FDR)
#View(mxe_rmats_ordered)
head(mxe_rmats_ordered)

# Subset the top rows of the ordered mxe_rmats data
mxe_rmats_subset <- mxe_rmats_ordered[1:2, ]

# List to store plots
plot_list <- list()

# Loop over each MXE event
for (i in 1:nrow(mxe_rmats_subset)) {
  event <- mxe_rmats_subset[i, ]
  
  # Coordinates
  chrom <- event$chr
  strand <- event$strand
  exon1_start <- event$X1stExonStart_0base
  exon1_end <- event$X1stExonEnd
  exon2_start <- event$X2ndExonStart_0base
  exon2_end <- event$X2ndExonEnd
  upstream_start <- event$upstreamEE
  upstream_end <- event$upstreamES
  downstream_start <- event$downstreamES
  downstream_end <- event$downstreamEE
  
  # Junction counts
  IJC1 <- event$IJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  SJC1 <- event$SJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  
  # Check for missing values in coordinates
  if (any(is.na(c(upstream_end, exon1_start, exon2_start, downstream_start)))) {
    message(paste("Skipping event", i, "due to missing coordinates"))
    next
  }
  
  # Sum junction counts
  IJC_sum <- sum(IJC1)
  SJC_sum <- sum(SJC1)
  
  # Exon boxes
  exons <- data.frame(
    xmin = c(upstream_start, exon1_start, exon2_start, downstream_start),
    xmax = c(upstream_end, exon1_end, exon2_end, downstream_end),
    ymin = -0.1,
    ymax = 0.1,
    label = c("Upstream", "Exon 1", "Exon 2", "Downstream")
  )
  
  # Junction arcs
  junctions <- data.frame(
    x = c(upstream_end, upstream_end),
    xend = c(exon1_start, exon2_start),
    y = 0.1,
    yend = 0.1,
    count = c(IJC_sum, SJC_sum),
    label = c("Inclusion", "Skipping")
  )
  
  # Plot
  plot <- ggplot() +
    geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#01B204", color = "black", size = 0.5) +
    geom_curve(data = junctions,
               aes(x = x, xend = xend, y = y, yend = yend),
               curvature = 0.3, linewidth = 1, color = "#B20101") +
    geom_text(data = junctions,
              aes(x = (x + xend) / 2, y = y + 0.15, label = count),
              size = 5, color = "black") +
    labs(
      title = paste("Sashimi-style Plot (MXE Event)\nGene:", event$geneSymbol),
      subtitle = paste("Chromosome:", chrom),
      x = "Genomic Position", y = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "none"
    )
  
  plot_list[[i]] <- plot
}

grid.arrange(grobs = plot_list, ncol = 2)

#################### ALTERNATIVE 3' SPLICE SITE EVENTS #########################

# Order the A3SS data by FDR
a3ss_rmats_ordered <- a3ss_rmats %>%
  arrange(FDR)

# Subset the top rows of the ordered A3SS data (adjust the range as needed)
a3ss_rmats_subset <- a3ss_rmats_ordered[1:2, ]  # Adjust based on your data

# List to store plots
plot_list <- list()

# Loop over each A3SS event
for (i in 1:nrow(a3ss_rmats_subset)) {
  event <- a3ss_rmats_subset[i, ]
  
  # Coordinates
  chrom <- event$chr
  strand <- event$strand
  exonStart <- event$longExonStart_0base
  exonEnd <- event$longExonEnd
  upstream_start <- event$flankingES
  upstream_end <- event$flankingEE
  downstream_start <- event$shortES
  downstream_end <- event$shortEE
  
  # Junction counts (from your data)
  IJC1 <- event$IJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  SJC1 <- event$SJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  
  # Check for missing values in coordinates
  if (any(is.na(c(upstream_end, exonStart, downstream_start)))) {
    message(paste("Skipping event", i, "due to missing coordinates"))
    next
  }
  
  # Sum counts
  IJC_sum <- sum(IJC1)
  SJC_sum <- sum(SJC1)
  
  # Define exon boxes
  exons <- data.frame(
    xmin = c(upstream_start, exonStart, downstream_start),
    xmax = c(upstream_end, exonEnd, downstream_end),
    ymin = -0.1,
    ymax = 0.1,
    label = c("Upstream", "Exon", "Downstream")
  )
  
  # Define junction arcs
  junctions <- data.frame(
    x = c(upstream_end, upstream_end),
    xend = c(exonStart, downstream_start),
    y = 0.1,
    yend = 0.1,
    count = c(IJC_sum, SJC_sum),
    label = c("Inclusion", "Skipping")
  )
  
  # Plot
  plot <- ggplot() +
    geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#8C019E", color = "black", size = 0.5) +
    geom_curve(data = junctions,
               aes(x = x, xend = xend, y = y, yend = yend),
               curvature = 0.3, linewidth = 1, color = "#B20101") +
    geom_text(data = junctions,
              aes(x = (x + xend) / 2, y = y + 0.15, label = count),
              size = 5, color = "black") +
    labs(
      title = paste("Sashimi-style Plot (A3SS Event)\nGene:", event$geneSymbol),
      subtitle = paste("Chromosome:", chrom),
      x = "Genomic Position", y = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "none"
    )
  
  plot_list[[i]] <- plot
}

grid.arrange(grobs = plot_list, ncol = 2)

#################### ALTERNATIVE 5' SPLICE SITE EVENTS #########################

# Order the A5SS data by FDR
a5ss_rmats_ordered <- a5ss_rmats %>%
  arrange(FDR)

# Subset the top rows of the ordered A5SS data (adjust the range as needed)
a5ss_rmats_subset <- a5ss_rmats_ordered[1:2, ]  # Adjust based on your data

# List to store plots
plot_list <- list()

# Loop over each A5SS event
for (i in 1:nrow(a5ss_rmats_subset)) {
  event <- a5ss_rmats_subset[i, ]
  
  # Coordinates
  chrom <- event$chr
  strand <- event$strand
  exonStart <- event$longExonStart_0base
  exonEnd <- event$longExonEnd
  upstream_start <- event$flankingES
  upstream_end <- event$flankingEE
  downstream_start <- event$shortES
  downstream_end <- event$shortEE
  
  # Junction counts 
  IJC1 <- event$IJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  SJC1 <- event$SJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  
  # Check for missing values in coordinates
  if (any(is.na(c(upstream_end, exonStart, downstream_start)))) {
    message(paste("Skipping event", i, "due to missing coordinates"))
    next
  }
  
  # Sum junction counts
  IJC_sum <- sum(IJC1)
  SJC_sum <- sum(SJC1)
  
  # Exon boxes
  exons <- data.frame(
    xmin = c(upstream_start, exonStart, downstream_start),
    xmax = c(upstream_end, exonEnd, downstream_end),
    ymin = -0.1,
    ymax = 0.1,
    label = c("Upstream", "Exon", "Downstream")
  )
  
  # Junction arcs
  junctions <- data.frame(
    x = c(upstream_end, upstream_end),
    xend = c(exonStart, downstream_start),
    y = 0.1,
    yend = 0.1,
    count = c(IJC_sum, SJC_sum),
    label = c("Inclusion", "Skipping")
  )
  
  # Plot
  plot <- ggplot() +
    geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#FDD91E", color = "black", size = 0.5) +
    geom_curve(data = junctions,
               aes(x = x, xend = xend, y = y, yend = yend),
               curvature = 0.3, linewidth = 1, color = "#B20101") +
    geom_text(data = junctions,
              aes(x = (x + xend) / 2, y = y + 0.15, label = count),
              size = 5, color = "black") +
    labs(
      title = paste("Sashimi-style Plot (A5SS Event)\nGene:", event$geneSymbol),
      subtitle = paste("Chromosome:", chrom),
      x = "Genomic Position", y = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "none"
    )
  
  plot_list[[i]] <- plot
}

library(gridExtra)
grid.arrange(grobs = plot_list, ncol = 2)


######################## RETAINED INTRON EVENTS ################################

# Order the RI data by FDR
ri_rmats_ordered <- ri_rmats %>%
  arrange(FDR)

# Subset the top rows of the ordered RI data
ri_rmats_subset <- ri_rmats_ordered[1:2, ]  # Adjust based on your data

# List to store plots
plot_list <- list()

# Loop over each RI event
for (i in 1:nrow(ri_rmats_subset)) {
  event <- ri_rmats_subset[i, ]
  
  # Coordinates
  chrom <- event$chr
  strand <- event$strand
  exonStart <- event$riExonStart_0base
  exonEnd <- event$riExonEnd
  upstream_start <- event$upstreamES
  upstream_end <- event$upstreamEE
  downstream_start <- event$downstreamES
  downstream_end <- event$downstreamEE
  
  # Junction counts
  IJC1 <- event$IJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  SJC1 <- event$SJC_SAMPLE_1 %>% strsplit(",") %>% unlist() %>% as.numeric()
  
  # Check for missing values in coordinates
  if (any(is.na(c(upstream_end, exonStart, downstream_start)))) {
    message(paste("Skipping event", i, "due to missing coordinates"))
    next
  }
  
  # Sum junction counts
  IJC_sum <- sum(IJC1)
  SJC_sum <- sum(SJC1)
  
  # Exon boxes
  exons <- data.frame(
    xmin = c(upstream_start, exonStart, downstream_start),
    xmax = c(upstream_end, exonEnd, downstream_end),
    ymin = -0.1,
    ymax = 0.1,
    label = c("Upstream", "Retention Intron", "Downstream")
  )
  
  # Junction arcs
  junctions <- data.frame(
    x = c(upstream_end, upstream_end),
    xend = c(exonStart, downstream_start),
    y = 0.1,
    yend = 0.1,
    count = c(IJC_sum, SJC_sum),
    label = c("Inclusion", "Retention")
  )
  
  # Plot
  plot <- ggplot() +
    geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#00FAB4", color = "black", size = 0.5) +
    geom_curve(data = junctions,
               aes(x = x, xend = xend, y = y, yend = yend),
               curvature = 0.3, linewidth = 1, color = "#B20101") +
    geom_text(data = junctions,
              aes(x = (x + xend) / 2, y = y + 0.15, label = count),
              size = 5, color = "black") +
    labs(
      title = paste("Sashimi-style Plot (RI Event)\nGene:", event$geneSymbol),
      subtitle = paste("Chromosome:", chrom),
      x = "Genomic Position", y = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "none"
    )
  
  plot_list[[i]] <- plot
}

grid.arrange(grobs = plot_list, ncol = 2)



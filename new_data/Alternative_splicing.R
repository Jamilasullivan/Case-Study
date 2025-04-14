## LOAD PACKAGES ##############################################################

#BiocManager::install("Gviz")

library(ggplot2)
library(dplyr)
library(ggforce)  # for arcs
library(tidyr)

## SET WORKING DIRECTORY ######################################################

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data")

## READ IN AS DATA FROM RMATS #################################################

## All data here is from the AIRP-CTRL 
## IncLevelDifference goes from condition 1 to condition 2 (<0 - higher exon skipping in condition 1. 0 = minimal difference in splicing). 

## skipped exon results
se_rmats <- read.table("SE.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## mutually exclusive exon results
mxe_rmats <- read.table("MXE.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## alternative 3' splice site results 
a3ss_rmats <- read.table("A3SS.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## alternative 5' splice site results
a5ss_rmats <- read.table("A5SS.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## retained intron results
ri_rmats <- read.table("RI.MATS.JC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## PREVIEW OF ALL RESULTS FILES ###############################################

head(se_rmats)
view(se_rmats)
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

############################## VOLCANO PLOTS ##################################

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

###############################################################################
#################### SELECTING AN EVENT TO PLOT ###############################
###############################################################################

event <- se_rmats[1, ] # selecting the first row of the se data to plot as the event (this is an example)

# Coordinates
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

# Sum counts across replicates
IJC_sum <- sum(IJC1)
SJC_sum <- sum(SJC1)

######################################################################################################### CREATING THE PLOT ##################################
###############################################################################

# Define exon rectangles
exons <- data.frame(
  xmin = c(upstreamExonStart, exonStart, downstreamExonStart),
  xmax = c(upstreamExonEnd, exonEnd, downstreamExonEnd),
  ymin = -0.1,
  ymax = 0.1,
  label = c("Upstream", "Skipped", "Downstream")
)

# Define arcs for junctions
junctions <- data.frame(
  x = c(upstreamExonEnd, upstreamExonEnd),
  xend = c(exonStart, downstreamExonStart),
  y = 0.1,
  yend = 0.1,
  count = c(IJC_sum, SJC_sum),
  label = c("Inclusion", "Skipping")
)

# Plot
ggplot() +
  # Exon boxes
  geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "skyblue") +
  
  # Arcs for splice junctions
  geom_curve(data = junctions, 
             aes(x = x, xend = xend, y = y, yend = yend), 
             curvature = 0.3, linewidth = 1, color = "black") +
  
  # Junction count labels
  geom_text(data = junctions, 
            aes(x = (x + xend) / 2, y = y + 0.15, label = count),
            size = 4) +
  
  theme_minimal() +
  labs(title = "Sashimi-style Plot (rMATS SE Event)",
       x = "Genomic Position", y = "") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())


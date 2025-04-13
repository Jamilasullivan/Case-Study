## LOAD PACKAGES ##############################################################

#BiocManager::install("Gviz")

library(ggplot2)
library(dplyr)
library(ggforce)  # for arcs
library(tidyr)

## SET WORKING DIRECTORY ######################################################

setwd("C:/Users/jamsu/OneDrive - Cardiff University/University/Masters/Big Data Biology/Modules/BIT103. Case Study/NEW DATA - WORKING DIRECTORY/Case-Study/new_data/data")

## READ IN AS DATA FROM RMATS #################################################

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
#view(se_rmats)
head(mxe_rmats)
#view(mxe_rmats)
head(a3ss_rmats)
#view(a3ss_rmats)
head(a5ss_rmats)
#view(a5ss_rmats)
head(ri_rmats)
#view(ri_rmats)

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
















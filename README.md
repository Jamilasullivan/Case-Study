# Case Study - SigLecF+ Neutrophil Analysis
**Big Data Biology** - Module BIT103

**Clients:** 
* Dr Kelly Berube
* Phoebe Ross

The task for this respository is to reanalyse the data collected in a previous study by Shin et al. (2022). Analysis was carried out from raw data processing, through to functional analysis.

The folder named *new_data* contains the analysis of the data mentioned above and all relevant scripts needed for future analysis. 

The folder named *old_data* contains the data and related files for the initially provided data that was not possible to analyse effectively (Pan et al. 2021).

## Raw Data Processing

The bash scripts used for this data processing can be found in the branch named *processing*. 

In this branch you will find the following scripts, which should be used in number order: 

1. 1-QC.sh
2. 2-star_index_genome.sh
3. 3-star.sh
4. 4-markduplicates.sh
5. 5-featurecounts.sh

The purpose of each of these scripts is outlined clearly in the report. 

## Differential Gene Expression Analysis

The scripts used for differential gene expression are *DESeq2.R* and *Limma.R*. Which is used is dependent on your data set. The circumstances under which each script is appropriate is outlined in the report.

These scripts are both capable of outputting results for use in the sebsequent scripts for filtering and visualisation of results.

### Tutorials that assisted in the production of the DESeq2.R script: 
* ["How I analyze RNA Seq Gene Expression data using DESeq2"](https://www.youtube.com/watch?v=kOlMcZujHHA)
* ["RNAseq volcano plot of differentially expressed genes"](https://www.youtube.com/watch?v=vRr78s37CI4)
* ["Draw Heatmap with Clusters Using pheatmap Package in R (4 Examples) | k-means, Row & Column Clusters"](https://www.youtube.com/watch?v=IjperDJ8IaI)
* ["deseq tutorial & visualization. how to plot dispersion estimates"](https://www.youtube.com/watch?v=6EiT5GF5rns)
* ["RNAseq tutorial – part 4 – Differential expression analysis with Deseq2"](https://www.youtube.com/watch?v=Ul-9s8YOOSk)

### Tutorials that assisted in the production of the Limma.R script:
* ["Differential Expression with Limma-Voom"](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html)
* ["DEG isolation using limma voom | A Rstudio Tutorial"](https://www.youtube.com/watch?v=z36fu178jIQ)

## Enrichment Analysis

The *Gene_enrichment_lists.R* gives the opportunity to filter data by specific values and output these results as files to be saved locally. However, this script is not specifically linked to the other enrichment analysis script. 

The *Enrichment_analysis.R* script is capable of assessing and providing visualisations for related processes and pathways to a given gene list.

### Tutorials that assisted in the production of the Enrichment_analysis.R script: 
* ["Overview of enrichment analysis"](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html)
* ["Gene Set Enrichment Analysis (+ R tutorial)"](https://www.youtube.com/watch?v=B7F7a9NcGS0)
* ["RNAseq analysis | Gene ontology (GO) in R"](https://www.youtube.com/watch?v=JPwdqdo_tRg)
* ["3 minute GSEA tutorial in R | RNAseq tutorials"](https://www.youtube.com/watch?v=Mi6u4r0lJvo)
* ["Gene Set Enrichment Analysis with ClusterProfiler"](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/)
* ["Visualization of Functional Enrichment Result"](https://bioconductor.riken.jp/packages/3.7/bioc/vignettes/enrichplot/inst/doc/enrichplot.html#references)

## Alternative Splicing Analysis

Alternative splicing analysis was carried out using the nf-core/rnasplice pipeline (Ashmore et al. 2024). 

The branch *alternative_splicing* contains the netflow config file used to run the pipeline and an annotated version to give a better idea of what the file is doing. 

## R Packages and Version Control
Analysis in R was completed using `R Studio 4.4.1`.

Below is a list of all of the R packages used between the .... scripts.

**Packages and Versions:**
* `edgeR 4.4.2`
* `limma 3.62.2`
* `ggplot 3.5.2`
* `EnhancedVolcano 1.24.0`
* `org.Mm.eg.db 3.20.0`
* `AnnotationDbi 1.68.0`
* `pheatmap 1.0.12`
* `ComplexHeatmap 2.22.0`
* `RColorBrewer 1.1-3`
* `DESeq2 1.46.0`
* `dplyr 1.1.4`
* `ggrepel 0.9.6`
* `clusterProfiler 4.14.6`
* `enrichplot 1.27.5`
* `cowplot 1.1.3`
* `pathview 1.46.0`
* `ReactomePA 1.50.0`
* `tibble 3.2.1`
* `GOplot 1.0.2`
* `stringr 1.5.1`
* `ggforce 0.4.2`
* `tidyr 1.3.1`
* `gridExtra 2.3`
* `ggupset 0.4.1`
  
## References

Ashmore, J. et al. 2024. nf-core/rnasplice: nf-core/rnasplice 1.0.4 (1.0.4). Zenodo. doi: https://doi.org/10.5281/zenodo.15194198 

Pan, S. et al. 2021. Enhanced Transcriptomic Resilience following Increased Alternative Splicing and Differential Isoform Production between Air Pollution Conurbations. *Atmosphere* 12(8), p. 959. Available at: https://www.mdpi.com/2073-4433/12/8/959 

Shin, J. W. et al. 2022. A unique population of neutrophils generated by air pollutant-induced lung damage exacerbates airway inflammation. *Journal of Allergy and Clinical Immunology* 149(4), pp. 1253-1269.e1258. Available at: https://doi.org/10.1016/j.jaci.2021.09.031doi: 10.1016/j.jaci.2021.09.031

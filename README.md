# Mites Analysis Code

This repository contains instruction code and R scripts for analyzing the microbial community of oribatid mites.

unmapping_blast.md 
- instructions to obtain unmapped reads, create BLASTn database, run BLASTn against unmapped reads and filter BLASTn results
- main tools: minimap2, samtools, blastn
- in bash 

length_quality.md
- commands for Seqkit and Nanoplot2: to get read length and read quality
- tools: seqkit, nanoplot
- in bash

script_import_blastn_2rstudio.R 
- R script to import the BLASTn results to RStudio and create a dataframe 

merge_all_blastn.R
- R script to merge all dataframes of all samples / BLASTn results to one dataframe

script_ggplots.R
- R script to visualize the bacteria hit counts and the bacteria hit proportion across all samples with ggplot2
- tool: ggplot2 

scrpt_diversity_test.R
- R script for alpha and beta diversity test
- Shannon Index, Bray-Curtis dissimilarity, heatmap, dendrogram, PERMANOVA, NMDS 

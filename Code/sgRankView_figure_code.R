setwd("~/yourpath/MAGeCK_figures/")
setwd("~/Documents/Projects/MAGeCK_figures")
library(ggplot2)
library(ggrepel)
library(xlsx)
library(stringr)
library(MAGeCKFlute)
source("./Code/Mageck_Figures_source.R")

## example data file 
file2 = "./data/Day14_vs_Day0_sample_rra.sgrna_summary.txt"



## Generate sgRNA Rank view plot
sgRankView(file2 = file2,
           neg = 'NonTargeting',
           file = "sgrankview_test.pdf")






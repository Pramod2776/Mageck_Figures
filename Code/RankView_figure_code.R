setwd("~/yourpath/MAGeCK_figures/")
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(xlsx)
source("./Code/Mageck_Figures_source.R")

## path to the gene summary file from MAGeCK package RRA analysis 
file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
                  "testdata/rra.gene_summary.txt")

## Generate RankView plot from file1 from the package
RankView(file1 = file1,
         filename = "./Figures/Rankview_test2.pdf")

## Generate RankView plot from experimental MAGeCK RRA analysis (example provided) 
file2 = readxl::read_excel("./data/example_gene_summary.xlsx") %>%
  as.data.frame()

## Generate RankView plot
RankView(file1 = file2,
         filename = "./Figures/Rankview_test1.pdf")




setwd("~/Documents/Projects/MAGeCK_figures")
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(poolr)
library(reshape2)
library(rstatix)
library(gginnards)
library(ggthemes)
library(ggeasy)
library(readxl)
library(dplyr)

## read example data with beta score difference and rank
example_data = readxl::read_excel("./data/betascore_mageck_figure.xlsx", sheet = 1) %>%
  as.data.frame()
## Genes tp label 
MLL3_4 = c("Utx", "Kmt2xc", "Kmt2d", "Ncoa6", "Pcgf1")
PRC1.1 = c("Ring1", "Phc1", "Kdm2b", "Bcor")

example_data$grouping <- "Other_complex"
example_data[example_data$Gene %in% toupper(MLL3_4), ]$grouping <- "MLL3_4"
example_data[example_data$Gene %in% toupper(PRC1.1), ]$grouping <- "PRC1.1"

## plot
options(ggrepel.max.overlaps = Inf)
g <-
  ggplot(example_data, aes(x = differencebeta, y = desc(Rank), fill = grouping, label = Gene))+
  geom_point(colour="white", shape=21, size = 4, 
             aes(fill = factor(grouping))) + 
  scale_fill_manual(values=c("yellow", "gray", "orange"))+
  theme_classic()+
  xlim(-2, 3.0)

g <- g + geom_point(
  pch = 21,
  colour = "black",
  fill = "grey",
  size = 4,
  alpha = 0.1
)

g <- g + geom_point(
  data = example_data[example_data$grouping == "MLL3_4",],
  aes(x = differencebeta, y = desc(Rank)),
  pch = 21,
  colour = "black",
  fill = "yellow",
  size = 4,
  alpha = 1.0,
  stroke = 1.0
)

g <- g + geom_point(
  data = example_data[example_data$grouping == "PRC1.1",],
  aes(x = differencebeta, y = desc(Rank)),
  pch = 21,
  colour = "black",
  fill = "orange",
  size = 4,
  alpha = 1.0,
  stroke = 1.0
)

g <- g +
  geom_text_repel(
    data = example_data[example_data$grouping == "MLL3_4",],
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.1,
    direction    = "y",
    ##angle        = 90,
    hjust        = 0,
    segment.size = 0.5,
    max.iter = 1e4, max.time = 1,
    min.segment.length = 1,
    ##size = 3.5,
    box.padding = 2.0,
    # point.padding = 0.5
    
  )

g <- g +
  geom_text_repel(
    data = example_data[example_data$grouping == "PRC1.1",],
    force_pull   = 1, # do not pull toward data points
    nudge_y      = 0.05,
    direction    = "x",
    ##angle        = 90,
    hjust        = 0,
    segment.size = 0.5,
    max.iter = 1e4, max.time = 1,
    min.segment.length = 1,
    box.padding = 3
  )

g <- g + xlab("Difference Beta score")
g <- g + ylab("Rank")
g = g+theme(axis.text.x = element_text(face="bold", color="black", 
                                       size=14, angle=0),
            axis.text.y = element_text(face="bold", color="black", 
                                       size=14, angle=0))+
  theme(text = element_text(size = 20)) +
  theme(axis.ticks.length=unit(.30, "cm"))+
  theme(legend.position = c(0.3, 0.7),
        legend.direction = "vertical")
g = g+
  ggtitle("")+
  ggeasy::easy_center_title()


ggsave(plot=g, filename="./Beta_score_gene_ranking.pdf", units = "in", width=8, height= 8, dpi=600)




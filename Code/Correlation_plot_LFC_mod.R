setwd("~/yourpath/MAGeCK_figures/")
library(ggplot2)
library(ggrepel)

## Read input data with LFC or Gene score or Z score for two comparisons or replicates
subset_edl = readxl::read_excel("./data/subset_edl.xlsx", sheet = 1) %>%
  as.data.frame()
rownames(subset_edl) = subset_edl$Gene

## subset data of interest for labeling and highlighting genes
MC38_high = subset_edl %>%
  arrange(desc(MC38_T_v_I_merge_z)) %>%
  slice(1:10) 

MC38_low = subset_edl %>%
  arrange(desc(MC38_T_v_I_merge_z)) %>%
  slice(291:nrow(subset_edl)) 

MC38 = rbind(MC38_high, MC38_low)

invitro_high = subset_edl %>%
  arrange(desc(invitro_c_v_a_merge_z)) %>%
  slice(1:10) 
invitro_low = subset_edl %>%
  arrange(desc(invitro_c_v_a_merge_z)) %>%
  slice(291:nrow(subset_edl)) 

invitro = rbind(invitro_high, invitro_low)

## code to generate the plots in slide 2 and 3
options(ggrepel.max.overlaps = Inf)
p = ggplot(subset_edl, aes(MC38_T_v_I_merge_z, invitro_c_v_a_merge_z, label = Gene)) +
  geom_text_repel(
    data          = MC38,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    ##direction    = "y",
    ##angle        = 90,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    min.segment.length = 0,
    box.padding = 0.5
  ) +
  geom_point(size = ifelse(subset_edl$Gene %in% MC38$Gene, 3, 0.5), color = ifelse(subset_edl$Gene %in% MC38$Gene, "red", "black"),
             alpha = ifelse(subset_edl$Gene %in% MC38$Gene, 1, 0.1), stroke = 1)+
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 3)+
  geom_hline(yintercept = 0, linetype = 3)+
  geom_text_repel(
    data          = invitro,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    ##direction    = "y",
    ##angle        = 90,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    min.segment.length = 0,
    box.padding = 0.5
  ) +
  geom_point(size = ifelse(subset_edl$Gene %in% invitro$Gene, 3, 0.5), color = ifelse(subset_edl$Gene %in% invitro$Gene, "#007575", "black"),
             alpha = ifelse(subset_edl$Gene %in% invitro$Gene, 1, 0.1), stroke = 1)+
  xlab("MC-38 Tumor vs Input LFC z-score") +
  ylab("Chronic vs Acute LFC z-score") +
  ggtitle("Correlation of in vivo z-score and in vitro z-scores for genes in the CRISPR")
  
ggsave(plot=p, filename="LFC_Genescore_Correlation.pdf", units = "in", width=12, height=8)  



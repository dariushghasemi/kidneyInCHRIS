#=========================================#
#  Visualizing Replication vs. Discovery
#=========================================#


# Started on October 2021
# Last update on February 20, 2025
library(tidyverse)

#-------------#
# inputs
egfr_log_res_EA <- "~/projects/kidneyInCHRIS/inputs/11-Dec-23_Suppl._Table_3_163_replicated_SNPs_in_CHRIS_EUA.csv"

# outputs
maf_effect_plot <- "~/projects/kidneyInCHRIS/outputs/21-Feb-25_Figure_3b_maf_vs_effect_ratio.png"

#-------------#
# read supplementary Table 3
repSNPs <- data.table::fread(egfr_log_res_EA)


#-----------------------------------------------------#
#------    Figure 3B: MAF vs. Effect Ratio     -------
#-----------------------------------------------------#

# Beta ratio vs. MAF ratio in CKDGen and CHRIS 
fig_3b <- repSNPs %>%
  # adding missing columns for plot
  dplyr::mutate(
    MAF_CKDGen = ifelse(EAF_CKDGen < .5, EAF_CKDGen, 1 - EAF_CKDGen),
    MAF_CHRIS  = ifelse(EAF_CHRIS  < .5, EAF_CHRIS,  1 - EAF_CHRIS),
    MAF_Ratio  = MAF_CHRIS / MAF_CKDGen,
    Beta_Ratio = Beta_CHRIS_ald / Beta_CKDGen_ald,
    Locus = str_replace(Locus, "PDILT", "UMOD-PDILT")# rename PDILT to UMOD-PDILT
  ) %>%
  ggplot(aes(x = MAF_Ratio, y = Beta_Ratio, shape = Locus, color = Locus)) + 
  geom_point(alpha = 0.9, size = 3.5) +
  #ggrepel::geom_text_repel(data = tagSNPs, aes(label = Locus), check_overlap = TRUE) +
  #ggfittext::geom_fit_text(grow = TRUE)+
  #geom_text(aes(label = Locus), check_overlap = TRUE) +
  #geom_abline(intercept = 0, slope = +1, color = "grey50", size = .9)+
  scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  scale_color_manual(values= c(
    "maroon1", "darkorchid2", "orange2", "green4", "steelblue2", 
     "tomato", "springgreen2", "royalblue2", "gold", "grey50", "darkturquoise"
    )) +
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(breaks = seq(0,6, 0.5)) +
  scale_x_continuous(breaks = seq(0.75,1.20, .05), limits = c(.75,1.20), expand = c(0,0)) +
  annotate("segment", x = 0.755, xend = 0.995, y = 6, yend = 6, colour = "grey50", size = 1.2,
           arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm"))) +
  annotate("segment", x = 1.005, xend = 1.195, y = 6, yend = 6, colour = "grey50", size = 1.2,
           arrow = arrow(ends = "both", angle = 90, length = unit(0.2, "cm"))) +
  annotate("text", x = 0.85, y = 6.3, size = 4,
           label = "Rarer alleles in CHRIS") + 
  annotate("text", x = 1.09, y = 6.3, size = 4, 
           label = "Rarer alleles in CKDGen") +
  # annotate("text", x = 0.23, y = 5.6, size = 4, label = "CASZ1") +
  # annotate("text", x = 0.7, y = 5.1, size = 4, label = "CASZ1") +
  # annotate("text", x = 0.59, y = 4, size = 4, label = "GAB2") +
  # annotate("text", x = 0.65, y = 3, size = 4, label = "DDX1") +
  # annotate("text", x = 1.15, y = 4, size = 4, label = "IGF1R") +
  # annotate("text", x = 1.2, y = 2.8, size = 4, label = "PIPK1B") +
  labs(
    x = "\nCHRIS-to-CKDGen MAF ratio\n",
    y = "\nCHRIS-to-CKDGen effect ratio\n"
    ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
    panel.grid.major.x = element_blank(),
    axis.text  = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.key.size  = unit(0.9, 'cm'),
    legend.key.width = unit(0.7, 'cm'),
    legend.text  = element_text(size = 12, face = "italic"),
    legend.title = element_text(size = 13, face = "bold")
    )

ggsave(maf_effect_plot, plot = fig_3b, width = 8, height = 5.5, dpi = 500, units = "in")

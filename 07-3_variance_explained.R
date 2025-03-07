#=========================================#
#  Comparing Replication with Discovery
#=========================================#

# Created on December 2023

library(tidyverse)
library(cowplot)
source("00-0_configuration.R") # to load significant loci in CHRIS "loci_CHRIS"


#------------#
# inputs
path_st2 <- "~/projects/kidneyInCHRIS/inputs/07-Dec-23_Suppl._Table2_147_CKDGen_Loci.csv"

# outputs
scatter_all <- "03-Dec-23_violin_plot_variance_explained_in_CHRIS_vs_CKDGen.png"
scatter_replicated <- "~/projects/kidneyInCHRIS/outputs/25-Feb-25_Figure_3a_variance_explained_in_147_loci_EA.png"
variance_maf_comparison <- "~/projects/kidneyInCHRIS/outputs/21-Feb-25_Figure_3_chris_vs_ckdgen.jpg"

#------------#
# read Suppl. Table 2
st2 <- data.table::fread(path_st2)


#-----------------------------------------------------#
#-------         Variance explained           --------
#-----------------------------------------------------#

# response to the comment of the Reviewer 1
# scatter plot of the explained variances
# using the formula used in CKDGen study

# This script uses the results attained 
# in `07-1_chris_vs_ckdgen.R` file: repSNPs_EA_Suppl2

#------------#
# compute variance explained in CHRIS and in CKDGen EU ancestry
st2_var_expl <- st2 %>%
  dplyr::rename(Beta_CKDGen = Beta_CKDGen_EA, Beta_CHRIS = Beta_CHRIS_tmp) %>%
  mutate(
    var_CHRIS  = (Beta_CHRIS  ^ 2) * ((2*EAF_CHRIS *(1 - EAF_CHRIS)) / 0.016),
    var_CKDGen = (Beta_CKDGen ^ 2) * ((2*EAF_CKDGen_EA *(1 - EAF_CKDGen_EA))/ 0.016),
    locus_me   = if_else(Locus %in% loci_CHRIS, "Significant", "Non-significant"),
    locus_rep  = if_else(Locus %in% loci_CHRIS, Locus, ""),
    locus_rep  = str_replace(locus_rep, "PDILT", "UMOD-PDILT")
    )

#------------#
# variance explained by each variant -> for slides
plt_varexpl_all <- st2_var_expl %>%
  ggplot(aes(var_CHRIS, var_CKDGen)) +
  geom_abline(slope = 1) +
  geom_point(color = "steelblue2", size = 2) +
  labs(x = "Variance explained by lead variant in CHRIS",
       y = "Variance explained by lead variant in CKDGen") +
  theme_classic() +
  theme(axis.text  = element_text(size = 8,  face = "bold"),
        axis.title = element_text(size = 12, face = "bold"))

ggsave(scatter_all, plot = plt_varexpl_all, width = 8, height = 5.5, dpi = 300, units = "in")

#------------#
# Fig. 3A: variance of replicated vs. not-replicated loci 
fig_3a <- st2_var_expl %>%
  ggplot(aes(var_CHRIS, var_CKDGen, color = locus_me, label = locus_rep)) +
  geom_abline(slope = 1) +
  geom_hline(yintercept = 0, color = "gray10", lty = 1)+
  geom_vline(xintercept = 0, color = "gray10", lty = 1)+
  geom_point(alpha = 0.8, size = 5) +
  #geom_text(aes(label = locus_rep), size=6, color = "gray20", fontface = 3, vjust = -.8, hjust=-.05) +
  ggrepel::geom_text_repel(size=7, fontface = 3, show.legend=F) +
  scale_color_manual(values = c("steelblue2", "royalblue4")) +
  scale_x_continuous(breaks = seq(0,0.0035, 0.0005), limits = c(0,0.0034)) +
  scale_y_continuous(breaks = seq(0,0.0020, 0.0005)) +
  labs(
    color = "Locus",
    x = "\nVariance explained in CHRIS",
    y = "Variance explained in CKDGen\n"
    ) +
  ylim(0, 0.0020) +
  theme(
    panel.background = element_blank(),
    axis.ticks.length = unit(0.2, 'cm'),
    legend.position = c(.2, .9),
    legend.key.size  = unit(0.9, 'cm'),
    legend.key.width = unit(0.7, 'cm'),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 13,  face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(l = 4, r = 10, t = 6, b = 4, unit = "mm")
    )


ggsave(scatter_replicated, plot = fig_3a, width = 9, height = 6.5, dpi = 500, units = "in")


#-----------------------------------------------------#
#-------            Figure 3 A&B              --------
#-----------------------------------------------------#

# join two plots for paper
fig3 <- cowplot::plot_grid(fig_3a, fig_3b, labels = c("A", "B"), ncol = 1, nrow = 2)

# save combined plot for paper
ggsave(variance_maf_comparison, plot = fig3, width = 10, height = 15, dpi = 500, units = "in")


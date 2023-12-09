#=========================================#
#  Comparing Replication with Discovery
#=========================================#

# Created on December 2023

#-----------------------------------------------------#
#-------         Variance explained           --------
#-----------------------------------------------------#

# response to the comment of the Reviewer 1
# scatter plot of the explained variances
# using the formula used in CKDGen study

# This script uses the results attained 
# in `07-1_chris_vs_ckdgen.R` file: repSNPs_EA_Suppl2

# compute variance in CHRIS and in CKDGen
repSNPs_EA_compar <- repSNPs_EA_Suppl2 %>%
  mutate(var_CHRIS  = (Beta_CHRIS_tmp ^ 2) * ((2*EAF_CHRIS *(1 - EAF_CHRIS)) / 0.016),
         var_CKDGen = (Beta_CKDGen_EA ^ 2) * ((2*EAF_CKDGen_EA *(1 - EAF_CKDGen_EA))/ 0.016)) 

#------------#
# variance explained by each variant
repSNPs_EA_compar %>%
  ggplot(aes(var_CHRIS, var_CKDGen)) +
  geom_abline(slope = 1) +
  geom_point(color = "steelblue2", size = 2) +
  labs(x = "Variance explained by lead variant in CHRIS",
       y = "Variance explained by lead variant in CKDGen") +
  theme_classic() +
  theme(axis.text  = element_text(size = 8,  face = "bold"),
        axis.title = element_text(size = 12, face = "bold"))

ggsave("03-Dec-23_violin_plot_variance_explained_in_CHRIS_vs_CKDGen.png", 
       last_plot(), width = 8, height = 5.5, dpi = 300, units = "in")

#------------#
# Fig. 3A: variance of replicated vs. not-replicated loci 
fig_3a <- repSNPs_EA_compar %>%
  mutate(`Locus replication` = if_else(Locus %in% unique(repSNPs_tmp$Locus),
                                       "Replicated", "Not-Replicated"),
         locus_rep = if_else(Locus %in% unique(repSNPs_tmp$Locus), Locus, "")) %>%
  ggplot(aes(var_CHRIS, var_CKDGen)) +
  geom_abline(slope = 1) +
  geom_hline(yintercept = 0, color = "gray40", lty = 2)+
  geom_vline(xintercept = 0, color = "gray40", lty = 2)+
  geom_point(aes(color = `Locus replication`), size = 2) +
  geom_text(aes(label = locus_rep), color = "gray20", fontface = 3, vjust = 1.6) +
  #ggrepel::geom_text_repel(aes(label = locus_rep),  check_overlap = TRUE) +
  scale_color_manual(values = c("steelblue2", "royalblue4")) +
  scale_x_continuous(breaks = seq(0,0.0035, 0.0005), limits = c(0,0.0034)) +
  scale_y_continuous(breaks = seq(0,0.0020, 0.0005)) +
  labs(x = "\nVariance explained by lead variants at 147 loci in CHRIS",
       y = "Variance explained by lead variants at 147 loci in CKDGen\n") +
  ylim(0, 0.0020) +
  theme_classic() +
  theme(legend.position = c(.3, .95),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text  = element_text(size = 11,  face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.margin = margin(l = 6, r = 6, t = 8, b = 6, unit = "mm"))

ggsave("05-Dec-23_violin_plot_variance_explained_in_147_loci_EA.png", 
       last_plot(), width = 9, height = 6.5, dpi = 300, units = "in")

#-----------------------------------------------------#
#-------            Figure 3 A&B              --------
#-----------------------------------------------------#

cowplot::plot_grid(fig_3a, fig_3b, labels = c("A", "B"), ncol = 1, nrow = 2)


ggsave("09-Dec-23_paper_revised_figure_3_chris_vs_ckdgen_hor.png",
       width = 9.5, height = 15, dpi = 400, units = "in")


#------------#
# violin plot: variance explained in CHRIS vs. CKDGen
repSNPs_EA_compar %>%
  pivot_longer(cols = c(var_CHRIS, var_CKDGen),
               names_to = c("quantity", "study"),
               names_pattern = "(var)_(CHRIS|CKDGen)$",
               values_to = "variance") %>%
  ggplot(aes(study, variance)) +
  geom_violin(aes(fill = study), color = NA, trim = F, show.legend = F)+
  geom_boxplot(color = "royalblue4", width = 0.05, outlier.shape = NA) +
  scale_fill_manual(values = c("steelblue2", "red3")) +
  labs(x = NULL, y = "Variance explained by lead variant at each locus") +
  theme_classic() +
  theme(legend.position = c(.9, .8),
        axis.text.x  = element_text(size = 12,  face = "bold"),
        axis.text.y  = element_text(size = 8,  face = "bold"),
        axis.title = element_text(size = 12, face = "bold"))


ggsave("03-Dec-23_violin_plot_variance_explained_in_CHRIS_vs_CKDGen.png", 
       last_plot(), width = 5, height = 5.5, dpi = 300, units = "in")


#------------#
# compare SNPs effect sizes in CHRIS vs. CKDGen
repSNPs_EA_compar %>%
  mutate(`Locus replication` = if_else(Locus %in% unique(repSNPs_tmp$Locus),
                                       "Replicated", "Not-Replicated"),
         locus_rep = if_else(Locus %in% unique(repSNPs_tmp$Locus), Locus, "")) %>%
  ggplot(aes(Beta_CHRIS_tmp, Beta_CKDGen,
             xmin = Beta_CHRIS_tmp - SE_CHRIS,
             xmax = Beta_CHRIS_tmp + SE_CHRIS)) +
  geom_abline(slope = 1) +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_vline(xintercept = 0, color = "gray20") +
  geom_errorbar(aes(ymin = Beta_CKDGen - SE_CKDGen,
                    ymax = Beta_CKDGen + SE_CKDGen), color = "steelblue3") +
  geom_pointrange(color = "royalblue4", alpha = .5) +
  geom_text(aes(label = locus_rep), size = 3.5, color = "gray20", fontface = 4, vjust = 1.6) +
  #geom_point(color = "steelblue4", size = 2) +
  labs(x = "Effect of lead variant in CHRIS",
       y = "Effect of lead variant in CKDGen") +
  theme_classic() +
  theme(axis.text  = element_text(size = 8,  face = "bold"),
        axis.title = element_text(size = 12, face = "bold"))

ggsave("03-Dec-23_linerange_plot_variants_effects_in_CHRIS_vs_CKDGen.png", 
       last_plot(), width = 13, height = 8.5, dpi = 300, units = "in")

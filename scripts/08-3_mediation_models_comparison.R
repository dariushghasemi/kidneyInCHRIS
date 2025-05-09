#=========================================#
# Mediation Analysis: Sensitivity analysis
#=========================================#

# The script created on December 2023
# to address comment of Reviewer 

library(tidyverse)
source("08-1_mediation_data.R") # to load vcfReg data file

#------------#
# It uses the function defined to run step 3
# of mediation analysis pipeline to benchmark
# the sensitivity of the model in step 1 to 
# the presence of PCs in step 3 rather
# than adjusting for kinship via Emmax

#------------#
# inputs
egfr_log_res_EA <- "~/projects/kidneyInCHRIS/inputs/10-Jan-23_replicatedSNPs_eGFRw.log.Res.txt"

# outputs
fig_s2 <- "03-Dec-23_scatter_plot_variant_effect_estimated_by_emmax_vs_lm_model.png"

#------------#
# replication analysis results: 163 SNPs in European Ancestry
repSNPs_tmp <- data.table::fread(egfr_log_res_EA, header = T, stringsAsFactors = FALSE)


#-----------------------------------------------------#
#------         Linear Regression Model        -------
#-----------------------------------------------------#

# Forming principal components terms
PCs <- paste0("PC", 1:10, collapse = " + ")

# Step 3 results: Long format
step1_lm <- getCoefs3lm(
  vcfReg["eGFRw.log.Res"],
  vcfReg[targets],
  paste("trait ~ SNP +", PCs),
  vcfReg
  )

summary(lm(data = vcfReg, eGFRw.log.Res ~ `chr5:177386403` + PC1 + PC2 + PC3 + PC4 + PC5))


#------------#
# Joint the lm model with emmax model results
res_emmax_lm <- repSNPs_tmp %>%
  dplyr::mutate(SNPid = str_c("chr", CHR, ":", POS38)) %>% # to join with lm model results
  dplyr::select(SNPid, Locus, Beta_CHRIS, SE_CHRIS, Pvalue_CHRIS) %>% #Beta_CHRIS_ald
  # Emmax models used in GWAS for step 1  
  dplyr::rename_with(
    .fn = ~ str_replace(.x, "CHRIS.+$|CHRIS", "emmax"),
    .cols = contains("CHRIS")
    ) %>%
  # join with the Linear Model result
  left_join(step1_lm, join_by(SNPid)) %>%
  dplyr::rename(Beta_lm = Estimate, SE_lm = SE, Pvalue_lm = Pvalue)


#-----------------------------------------------------#
#------             Suppl. Figure S2           -------
#-----------------------------------------------------#

# compare the effect size of variants via scatter plot
plt_fig_s2 <- res_emmax_lm %>%
  ggplot(aes(Beta_emmax, Beta_lm, shape = Locus, color = Locus)) +
  geom_abline(slope = 1) +
  geom_vline(xintercept = 0, color = "gray20") +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  scale_color_manual(values = c(
    "maroon1", "darkorchid2", "orange2", "green4", "steelblue2", "darkturquoise",
    "tomato", "springgreen2", "royalblue2", "gold", "grey50")
    ) +
  labs(
    x = "Effect of variant via mixed model (emmax) adjusted for kinship",
    y = "Effect of variant via linear model (lm) adjusted for 10 PCs"
    ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
    panel.grid.major.x = element_blank(),
    axis.text  = element_text(size = 8,  face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    #legend.position = "none",
    legend.key.size  = unit(0.99, 'cm'),
    legend.key.width = unit(0.7, 'cm'),
    legend.text  = element_text(size = 12, face = "italic"),
    legend.title = element_text(size = 13, face = "bold")
    )


ggsave(fig_s2, plot = plt_fig_s2, width = 8, height = 5.5, dpi = 300, units = "in")

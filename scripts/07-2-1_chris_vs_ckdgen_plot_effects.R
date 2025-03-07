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
effects_rangeplot <- "~/projects/kidneyInCHRIS/outputs/25-Feb-25_Figure_2_effect_comparison_leading_variants.png"

#-------------#
# read supplementary Table 3
repSNPs <- data.table::fread(egfr_log_res_EA)



#-----------------------------------------------------#
#------          Figure 2: Line range           ------
#-----------------------------------------------------#

# The difference of Effect size in CHRIS and CKDGen
effect_comparison <- repSNPs %>%
  dplyr::filter(
    SNPid %in% (target_snps %>% str_remove("chr"))
    ) %>%
  dplyr::select(
    SNPid, RSID, Locus, Beta_CHRIS_ald, Beta_CKDGen_ald, SE_CHRIS, SE_CKDGen
    ) %>%
  dplyr::rename(
    Beta_CHRIS  = Beta_CHRIS_ald, 
    Beta_CKDGen = Beta_CKDGen_ald
    ) %>%
  pivot_longer(
    cols = !c(SNPid, RSID, Locus),   
    names_to = c("trait", "study"),
    names_pattern = "(.+)_(CHRIS|CKDGen)$",
    values_to = c("score")
    ) %>%
  pivot_wider(
    names_from = "trait",
    values_from = "score"
    ) %>%
  dplyr::mutate(label_x = paste(RSID, Locus))



# create effect plot
plt_effect_comparison <- effect_comparison %>%
  ggplot(
    aes(
      x = locus_factor_rsid(Locus),
      y = Beta,
      color = study,
      ymin = Beta - 1.95*SE,
      ymax = Beta + 1.95*SE
    )) +
  geom_pointrange(
    aes(color = study),
    size = 1.1,
    fatten = 2.5,
    alpha = 0.9,
    #lineend = "round",
    position = position_dodge(.3)
    ) +
  #scale_color_grey(start=0.55, end=0.25) +
  scale_color_manual(values=c('turquoise3','tomato3'))+
  #scale_size_manual(values=c(2, 4))+
  geom_hline(yintercept = 0, color="black", lty = 1, size = .5) +
  scale_y_continuous(breaks = seq(-.035, .001, .005)) +
  labs(y = "Effect on ln(eGFRcrea)\n", x = NULL) +
  coord_cartesian(ylim = c(-0.035, 0.001)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(.9, .3),
    legend.key.size  = unit(0.9, 'cm'),
    legend.key.width = unit(0.7, 'cm'),
    legend.text  = element_text(size = 12, face = "plain"),
    axis.ticks.length = unit(0.2, 'cm'),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold.italic", angle=90, hjust=0, vjust=.5),
    axis.title  = element_text(size = 13, face = "bold")
    )

ggsave(effects_rangeplot, plot = plt_effect_comparison, width = 7.5, height = 5.5, dpi = 500, units = "in")

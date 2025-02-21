#=========================================#
#    Mediation Analysis: Visualization
#=========================================#

library(tidyverse)
source("00-0_configuration.R") # to load functions and vectors

#------------#

# inputs
path_med_res <- "~/projects/kidneyInCHRIS/inputs/11-Jan-23_SNP-wise_summary_of_mediation_analysis_steps.csv"

# outputs
outliers_plot <- "~/projects/kidneyInCHRIS/outputs/21-Feb-25_Figure_4_linerange_lead_variants.jpg"

#------------#
# read results of mediation analysis steps
res_mediation <- data.table::fread(path_med_res)



#-----------------------------------------------------#
#------          Figure 4: line range           ------
#-----------------------------------------------------#

# Reorder Trait based on frequency of outlier effect
outlier_freq <- res_mediation %>%
  dplyr::filter(SNPid %in% target_snps) %>% # taking lead variants
  group_by(Trait, outlier) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    id_cols = Trait,
    names_from = outlier,
    values_from = count,
    names_prefix = "outlier_"
    ) %>%
  mutate_all(~replace_na(., 0)) %>%
  mutate(outlier_count = fct_reorder(factor(Trait), outlier_Yes, .desc = T)) %>%
  arrange(outlier_count)


#---------#
#Step2: Line range plot for the paper
plt_outlier <- res_mediation %>%
  dplyr::filter(SNPid %in% target_snps) %>% # CKDGen lead variants or the most significant SNP in CHRIS
  right_join(outlier_freq, by = "Trait") %>% # add frequency of trait with outlier effect
  dplyr::mutate(
    Trait = recode(
      Trait,
      "Body_Fat"   = "Body Fat",
      "Visceral_Fat" = "Visceral Fat",
      "Pulse_Rate" = "Pulse Rate",
      "INR_PT"     = "INR PT",
      "APTT"       = "aPTT",
      "APTT_ratio" = "aPTT ratio",
      "AST_GOT"    = "AST GPT",
      "ALT_GPT"    = "ALT GPT",
      "Urine_pH"   = "Urine pH"
      ),
    trait_ordered  = fct_reorder(Trait, outlier_Yes, .desc = F), # y-axis labels
    Locus_reordered = locus_factor_rsid(Locus) # panel headers
    ) %>%
  ggplot(aes(
    x = Estimate_Step2,
    y = trait_ordered,
    xmin = Estimate_Step2 - SE_Step2,
    xmax = Estimate_Step2 + SE_Step2
    )) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50")+
  geom_pointrange(aes(color = outlier), show.legend = F, fatten = 1.2) +
  scale_color_manual(values = c("grey50", "tomato3")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  coord_cartesian() +
  facet_wrap(~ Locus_reordered, nrow = 1, scales = "free_x", shrink = T)+
  labs(x = "Trait-adjusted effect of CKDGen lead SNP at each locus") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size = 6,  face = "bold"),
    axis.text.y = element_text(size = 6,  face = "bold"),
    axis.title  = element_text(size = 10, face = "bold"),
    #axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray", size = .25),
    strip.text.x = element_text(size = 8, color = "Black", face = "bold.italic")
    )


ggsave(outliers_plot, plot = plt_outlier, width = 15, height = 9, dpi = 500, units = "in")

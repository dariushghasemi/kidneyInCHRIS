#=========================================#
#    Mediation Analysis: Visualization
#=========================================#

library(tidyverse)
source("00-0_configuration.R") # to load functions and vectors

#------------#

# inputs

# outputs
out_barplot <- "08-Oct-22_Frequency of the SNPs with outlier trait in step 2.png"
out_frequency <- "08-Oct-2022_outlier traits in Step2 grouped by trait.csv"
outliers_combined <- '12-Jan-23_Linerange plot of the entire 163 replicated SNPs.pdf'


#------------#
#Outlier traits grouped by Locus and SNPid
results_Step2_long %>% 
  #mutate(outlierTrait = ifelse(outlier == "Yes", trait, "")) %>%
  filter(outlier == "Yes") %>%
  group_by(Locus) %>%
  count(Trait) %>%
  # write.csv(out_frequency, row.names = FALSE)
  ggplot(aes(y = forcats::fct_rev(forcats::fct_infreq(Trait)), fill = Locus)) +
  geom_bar(width = .85,  position = "dodge")+
  scale_x_continuous(breaks = seq(1,10, 1), limits = c(0,10.2), expand = c(0,0)) +
  #scale_fill_brewer(palette = "Paired")+#"Accent"
  #scale_fill_viridis_d(option = "magma")+
  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.y  = element_text(size = 10,  face="bold"),
        legend.position = c(.88, .4),
        legend.key.size = unit(0.8, 'cm'),
        panel.grid.major.x = element_line(color = "lightgray", size = .25)) +
  labs(x = "Frequency of the SNPs with outlier trait in step 2", y = NULL)

ggsave(out_barplot, last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")



#-----------------------------------------------------#
#---------      line range plot per SNP      ---------
#-----------------------------------------------------#

# This is for figure 4 in case we need to add GWAS beta to the plot.
sum3steps_long_eGFR <- sum3steps_long %>%
  # Taking effect size of the top SNPs in GWAS
  filter(SNPid %in% tagSNPs$SNPid, Trait == "SCr") %>%
  mutate(Locus_reordered = Locus_factor(Locus))

#---------#
#Line range plot of adjusted effect size (mediation step 2) for the entire 163 variants
pdf(outliers_combined, width=8, height = 8)

map(1:length(targets), function(i){
  results_Step2_long %>%
    filter(SNPid %in% targets[i]) %>%
    mutate(Trait = recode(Trait,
                          "Body_Fat"   = "Body Fat",
                          "Visceral_Fat" = "Visceral Fat",
                          "Pulse_Rate" = "Pulse Rate",
                          "INR_PT"     = "INR PT",
                          "APTT_ratio" = "APTT ratio",
                          "AST_GOT"    = "AST GPT",
                          "ALT_GPT"    = "ALT GPT",
                          "Urine_pH"   = "Urine pH"),
           Locus_reordered = locus_factor_snpid(Locus),
           Feature = fct_reorder(Trait, Estimate, .desc = F)) %>%
    ggplot(aes(Estimate, 
               Feature,
               xmin = Estimate - SE, 
               xmax = Estimate + SE)) +
    geom_vline(data = repSNPs[i,],
               aes(xintercept = Beta_CHRIS),
               color = "steelblue3",
               lty = 1) +
    geom_vline(xintercept = 0, lty = 2, color = "grey50")+
    geom_pointrange(aes(color = outlier), show.legend = F, fatten = 1.2) +
    scale_color_manual(values = c("grey50", "tomato3")) +
    #scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    coord_cartesian() +
    labs(x = paste0("Trait-adjusted effect of the ", targets[i], " SNP in *", repSNPs[i,"Locus"], "* locus on ln(eGFRcreat)")) + # in SHROOM3 Locus
    theme_classic() +
    theme(strip.background = element_blank(), strip.placement = "outside",
          axis.text.x = element_text(size = 6,  face = "bold"),
          axis.text.y = element_text(size = 6,  face = "bold"),
          axis.title  = element_text(size = 10, face = "bold"),
          axis.title.x = ggtext::element_markdown(),
          axis.title.y.left = element_blank(),
          panel.grid.major.y = element_line(color = "lightgray", size = .25),
          strip.text.x = element_text(size = 8, color = "Black", face = "bold.italic"))
}
)

dev.off()
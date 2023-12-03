#=========================================#
#    Mediation Analysis: Visualization
#=========================================#

library(tidyverse)
#------------#

#Outlier traits grouped by Locus and SNPid
results_Step2_long %>% 
  #mutate(outlierTrait = ifelse(outlier == "Yes", trait, "")) %>%
  filter(outlier == "Yes") %>%
  group_by(Locus) %>%
  count(Trait) %>% View()
#write.csv("08-Oct-2022_outlier traits in Step2 grouped by trait.csv", row.names = FALSE)
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

ggsave("08-Oct-22_Frequency of the SNPs with outlier trait in step 2.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")



#-----------------------------------------------------#
#---------------- Figure 4: line range ---------------
#-----------------------------------------------------#

#New line range plot for the paper
#Changing the order and the label of the Locus top SNPs
Locus_factor <- function(x) {
  factor(x,
         levels = c("CASZ1","DDX1","SHROOM3","DAB2","SLC34A1",
                    "TMEM60","STC1","PIP5K1B","GAB2","IGF1R","PDILT"),
         labels = c("CASZ1\n1:10670853",
                    "DDX1\n2:15642347",
                    "SHROOM3\n4:76480299",
                    "DAB2\n5:39393631",
                    "SLC34A1\n5:177386403",
                    "TMEM60\n7:77714744",
                    "STC1\n8:23885208",
                    "PIP5K1B\n9:68817258",
                    "GAB2\n11:78324786",
                    "IGF1R\n15:98733292",
                    "PDILT\n16:20381010"))
}

# with rsID
Locus_factor <- function(x) {
  factor(x,
         levels = c("CASZ1","DDX1","SHROOM3","DAB2","SLC34A1",
                    "TMEM60","STC1","PIP5K1B","GAB2","IGF1R","PDILT"),
         labels = c("CASZ1\nrs74748843",
                    "DDX1\nrs807624",
                    "SHROOM3\nrs28817415",
                    "DAB2\nrs10062079",
                    "SLC34A1\nrs3812036",
                    "TMEM60\nrs57514204",
                    "STC1\nrs819196",
                    "PIP5K1B\nrs2039424",
                    "GAB2\nrs7113042",
                    "IGF1R\nrs59646751",
                    "PDILT\nrs77924615"))
}

#---------#
# Taking effect size of the top SNPs in GWAS
# by subsetting summary of Mediation-A
sum3steps_long_eGFR <-
  sum3steps_long %>%
  filter(SNPid %in% tagSNPs$SNPid,
         Trait == "SCr") %>% 
  mutate(Locus_reordered = Locus_factor(Locus)) #%>% View()

#---------#
# Reorder Trait based on frequency of 
# outlier effect for the lead SNPs
results_Step2_long_freq <- 
  results_Step2_long %>%
  filter(SNPid %in% tagSNPs$SNPid) %>%
  group_by(Trait, outlier) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(id_cols = Trait,
              names_from = outlier,
              values_from = count,
              names_prefix = "outlier_") %>%
  mutate_all(~replace_na(., 0)) %>%
  mutate(Feature_freq = fct_reorder(factor(Trait), outlier_Yes, .desc = T)) %>%
  arrange(Feature_freq) %>% View

#---------#
#Step2: Line range plot of the CKDGen lead variants for the paper

results_Step2_long %>%
  filter(SNPid %in% tagSNPs$SNPid) %>%
  right_join(results_Step2_long_freq, by = "Trait") %>%
  mutate(Trait = recode(Trait,
                        "Body_Fat"   = "Body Fat",
                        "Visceral_Fat" = "Visceral Fat",
                        "Pulse_Rate" = "Pulse Rate",
                        "INR_PT"     = "INR PT",
                        "APTT"       = "aPTT",
                        "APTT_ratio" = "aPTT ratio",
                        "AST_GOT"    = "AST GPT",
                        "ALT_GPT"    = "ALT GPT",
                        "Urine_pH"   = "Urine pH"),
         Feature         = fct_reorder(Trait, Estimate,    .desc = F),
         Feature_freq    = fct_reorder(Trait, outlier_Yes, .desc = F),
         Locus_reordered = Locus_factor(Locus)) %>%
  ggplot(aes(Estimate,
             Feature_freq,
             xmin = Estimate - SE,
             xmax = Estimate + SE)) +
  # geom_vline(data = sum3steps_long_eGFR,
  #            aes(xintercept = Estimate_GWAS), 
  #            color = "steelblue3",
  #            lty = 1) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50")+
  geom_pointrange(aes(color = outlier), show.legend = F, fatten = 1.2) +
  scale_color_manual(values = c("grey50", "tomato3")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  coord_cartesian() +
  facet_wrap(~ Locus_reordered, nrow = 1, scales = "free_x", shrink = T)+
  labs(x = "Trait-adjusted effect of CKDGen lead SNP at each locus") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside",
        axis.text.x = element_text(size = 6,  face = "bold"),
        axis.text.y = element_text(size = 6,  face = "bold"),
        axis.title  = element_text(size = 10, face = "bold"),
        #axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray", size = .25),
        strip.text.x = element_text(size = 8, color = "Black", face = "bold.italic"))

#ggsave("Top_SNPs_LinerangePlot_CASZ1_10-Apr-22.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

ggsave("28-Mar-23_Linerange plot for the leading SNP in each locus.png",
       last_plot(), width = 12, height = 7, pointsize = 2, dpi = 300, units = "in")

#---------#

#Step2: Line range plot of the entire 163 variants
pdf('12-Jan-23_Linerange plot of the entire 163 replicated SNPs.pdf', width=8, height = 8)

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
           Locus_reordered = Locus_factor(Locus),
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


#---------#

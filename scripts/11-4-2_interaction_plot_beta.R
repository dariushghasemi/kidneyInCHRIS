

library(pheatmap)
library(gridExtra)


#---------#

pdf('02-Feb-23_interaction comparisons.pdf', width=18, height = 18)

map2(c("Dosage_estimate_nonCen", 
       "TSH_cen_estimate", 
       "`Dosage:TSH_nonCen_estimate`", 
       "Dosage_p.value_nonCen.log", 
       "TSH_nonCen_p.value.log", 
       "Dosage_TSH_nonCen_p.value.log"),
     c("Dosage_estimate_Cen", 
       "TSH_cen_estimate",
       "`Dosage:TSH_cen_estimate`", 
       "Dosage_p.value_Cen.log", 
       "TSH_cen_p.value.log", 
       "Dosage_TSH_cen_p.value.log"),
     function(myvar1, myvar2) {
       results_Kidney_TSH_ %>%
         rename_with(~ gsub("TSH", "TSH_nonCen", .x, fixed = TRUE))%>%
         right_join(results_Kidney_TSH_cen,
                    by = c("Locus", "SNPid"),
                    suffix = c("_nonCen", "_Cen")) %>%
         write.csv("02-Feb-23_Kideny interaction with and without centralized thyroid (SNP-TSH)_wide format.csv", row.names = FALSE, quote = FALSE)
       mutate(Dosage_p.value_nonCen.log = -log10(Dosage_p.value_nonCen),
              TSH_nonCen_p.value.log    = -log10(TSH_nonCen_p.value),
              Dosage_TSH_nonCen_p.value.log = -log10(`Dosage:TSH_nonCen_p.value`),
              Dosage_p.value_Cen.log    = -log10(Dosage_p.value_Cen),
              TSH_cen_p.value.log       = -log10(TSH_cen_p.value),
              Dosage_TSH_cen_p.value.log= -log10(`Dosage:TSH_cen_p.value`)) %>%
         ggplot(aes_string(x = myvar1, y = myvar2, color = "Locus")) +
         geom_abline(lty = 2) +
         geom_vline(xintercept = 0) +
         geom_hline(yintercept = 0) +
         geom_point(size = 2.5) +
         labs(x = paste(myvar1, "in ln(eGFR) ~ Dosage * TSH + Sex + Age + PC1:10"),
              y = paste(myvar2, "in ln(eGFR) ~ Dosage * TSH_cen + Sex + Age + PC1:10")) +
         theme_classic()
     })

dev.off()

# ggsave("02-Feb-23_Effect of SNP in model with TSH vs centralized TSH.png",
#        width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")




#-----------------------------------------------------#
#-------                   Other                ------
#-----------------------------------------------------#

# retrieve & save a printed but unsaved console output history by: .Last.value
View(.Last.value)  
history(max.show = Inf) # full history
#---------#

#Plotting and comparing Betas and Pvalues of interaction terms
resModel3 %>% 
  as.data.frame() %>% 
  mutate(coef = rownames(resModel3)) %>% 
  rename(Beta = Estimate, SE = `Std. Error`, Pvalue = `Pr(>|t|)`) %>% 
  mutate(SNPid = rep(repSNPs$SNPid, each = 5),
         Locus = rep(repSNPs$Locus, each = 5)) %>% 
  mutate(coef = replace(coef, coef == "TSH_catHyperT",     "HyperThyrodism"),
         coef = replace(coef, coef == "TSH_catHypoT",      "HypoThyrodism"),
         coef = replace(coef, coef == "SNP:TSH_catHyperT", "SNP_HyperThyrodism"),
         coef = replace(coef, coef == "SNP:TSH_catHypoT",  "SNP_HypoThyrodism")) %>%
  mutate(tagHyper = if_else(coef == "SNP_HyperThyrodism", "GAB2", ""),
         tagHypo  = if_else(coef == "SNP_HypoThyrodism",  "GAB2      IGF1R, PIP5K1B", "")) %>% 
  filter(coef %in% c("SNP_HyperThyrodism", "SNP_HypoThyrodism")) %>% 
  ggplot(aes(x = Pvalue)) +
  geom_density(data = data.frame("rnd" = runif(100000)),
               aes(x = rnd),
               fill = "grey70",
               color = "grey50")+
  geom_density(aes(fill = coef),
               alpha = 0.6,
               color = "grey50") + #, show.legend = NULL
  scale_fill_manual(values = c("violetred1", "turquoise2"),
                    #name   = "SNP vs",
                    labels = c("Hyperthyrodism",
                               "Hypothyrodism")) + 
  geom_vline(aes(xintercept = 0.05), lty = 2) +
  geom_text(aes(label = tagHyper),
            x = 0.14,
            y = 7,
            color = "violetred1",
            size = 4,
            fontface = 2) +
  geom_text(aes(label = tagHypo),
            x = 0.68,
            y = 3,
            color = "turquoise2",
            size = 4,
            fontface = 2) +
  labs(title = "Density plot of the p-values for interaction of SNPs vs:",
       x = "Pvalue for the interaction term") + 
  facet_grid(coef ~ .) + #xlim(0, 0.0025) + scale_fill_discrete() +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(), #"",
        legend.position = "top", #"none", "bottom"
        axis.text.y  = element_text(size = 10, face = "bold"),
        axis.text.x  = element_text(size = 10, face = "bold"),
        axis.title   = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.6, 'cm'),
        legend.text  = element_text(size = 10, face = "bold"),
        strip.text = element_blank()
  )

ggsave("07-Jun-222_SNP-Hyperthyroid_vs_SNP-Hypothyroid_interactions_Pvalues_ESHG_poster.png", 
       width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-----------------------------------------------------#

#Plots
ggplot(resAllModels, aes(y = (SNP_Model0 - T3_beta)))+
  geom_boxplot(alpha = 0.6, fill = "Orange", outlier.color = "red")
#geom_text(aes(label = Locus, color = Locus), size = 3, vjust = -1)

#scatter
ggplot(resAllModels, aes(x = TSH_beta , y = TSH.q_beta))+
  geom_point(aes(shape = Locus, color = Locus), size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) +
  geom_abline(slope = 1, intercept = 0, color="gray", linetype="dashed", size = 1.1) +
  theme_classic()
#geom_text(aes(label = Locus, color = Locus), size = 2, vjust = -1)

ggplot(resAllModels, aes(x = (SNP_Model0 - T3_beta) , y = (Effect.ckdgen)))+
  geom_point(aes(shape = Locus, color = Locus), size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) + theme_classic()
#geom_abline(slope = 1, intercept = 0, color="gray", linetype="dashed", size = 1.1) + 

ggsave("SNP_Model0-T3_beta) vs. (Effect.ckdgen).png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#---------#


pdf('pheatmap_MediationAnalysis_TSH_cat_25Jan.pdf', width=18, height = 18)

#Betas
pheatmap(medAna_Merged[c("BETA", "SNP.Beta2", "SNP.Beta5")] / medAna_Merged$BETA,
         cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = T, labels_row = medAna_Merged$SNPid, 
         border_color = 'Black', fontsize_row = 6, fontsize_col = 12, angle_col = "0")
#P-values
pheatmap(medAna_Merged[c("PVALUE", "SNP.Pvalue2", "SNP.Pvalue5")],
         cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = T, labels_row = medAna_Merged$SNPid, 
         border_color = 'Black', fontsize_row = 6, fontsize_col = 12, angle_col = "0")
#---------#

#Scatter plots of Betas
p11 <- medAna_Merged %>% 
  ggplot(aes(BETA, SNP.Beta2)) + 
  geom_point(aes(shape = Locus, color = Locus), size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) + 
  geom_abline(slope = 1, intercept = 0, color="gray") + theme_classic()

p12 <- medAna_Merged %>% 
  ggplot(aes(BETA, SNP.Beta5)) + 
  geom_point(aes(shape = Locus, color = Locus), size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) + 
  geom_abline(slope = 1, intercept = 0, color="gray") + theme_classic()

# join plots
p13 <- gridExtra::grid.arrange(p11, p12, ncol = 1)

dev.off()

#=========================================#
#               Proteomics
#=========================================#


# ------ Package installation ------

# Look for the path where the packages are installed.
# to change the path use this -> .libPaths( "/Users/tex/lib/R" )
.libPaths()

# Then upload chrisProteomicsData package tar file there and install it manually -> 
# "C:/Users/dghasemisemeskandeh/Documents/R/win-library/4.1"
#install.packages("~/R/win-library/4.1/chrisProteomicsData_0.3.0.tar.gz", repos = NULL, type = "source")

library(chrisProteomicsData)
library(chrisProteomicsData, quietly = TRUE, warn.conflicts = FALSE)

browseVignettes("chrisProteomicsData")


#Download required packages
#install.packages("BiocManager")
BiocManager::install("QFeatures")

#-----------------------------------------------------#
# Reading plasma proteins data
data("chris_plasma_proteins")

rowData(chris_plasma_proteins)

# Store protein names
proteinsInfo <- rowData(chris_plasma_proteins)

# contains AID
colData(chris_plasma_proteins) %>% View()
table(chris_plasma_proteins$sampleOR)

# Protein concentration
chrisProteomy <- 
  assay(chris_plasma_proteins) %>%
  as.data.frame %>%
  rownames_to_column() %>%
  pivot_longer(-rowname, names_to = "sample_id") %>% 
  pivot_wider(names_from=rowname, values_from=value) %>% 
  left_join(by = "sample_id",
            tibble("sample_id" = colData(chris_plasma_proteins)$sample_id,
                   "AID"       = colData(chris_plasma_proteins)$AID),
            .) %>% as.data.frame()

#Saving metabolites names for mediation analysis
proteins <- chrisProteomy %>% colnames %>% setdiff(.,c("AID", "sample_id"))

#Merge 148 proteins for 4,087 with CHRIS baseline
vcfProteomy <-
  chris[c("AID", "Age", "Sex", "eGFRw.log.Res")] %>%
  inner_join(PCs_13K,   by = "AID") %>% #dim()
  inner_join(chrisProteomy, by = "AID") %>% #dim()
  inner_join(vcfmod,    by = "AID") #%>% dim()

#-----------#
pool <- chris_plasma_proteins[, chris_plasma_proteins$sampleOR == "POOL"]
ncol(pool)
chris_plasma_proteins$colData %>% View()
assay(pool) %>% dim
boxplot(assay(pool)) 
grid()
chrisProteomy %>% View

#-----------#
rowData(chris_plasma_proteins) %>% 
  as.data.frame()  %>%  tibble::rownames_to_column(var="PROT_NAME") %>% View
  mutate(PROT_NAME=str_replace_all(PROT_NAME, ":", "_"), 
         PP = str_replace_all(PROT_NAME,";", "_"),
         PP = str_replace_all(PP, "-", ":"),
         PP2 = str_replace_all(PP,":", ".")) %>% View
  
#-----------#
library(QFeatures)

data("chris_plasma_proteomics")
rowData(chris_plasma_proteomics)$peptides %>% View

assay(chris_plasma_proteomics) %>% rownames
assay(chris_plasma_proteomics) %>% dim
assay(chris_plasma_proteomics, "peptides") %>% dim
rowData(chris_plasma_proteomics)$peptides  %>% dim
rowData(chris_plasma_proteomics)$proteins %>% dim

#-----------------------------------------------------#
# ---------------- Mediation Analysis ----------------
#-----------------------------------------------------#


#-----------#
#- QC step 1 - 

summary(lm(eGFRw.log.Res ~ `chr1:10599281`, data = vcfMass))

vcfReg %>%
  mutate(Peptited = ifelse(AID %in% vcfProteomy$AID, 1, 0)) %>% 
  count(Peptited)

# Indicator variable showing if the individual has proteins
vcfReg$Peptited <- case_when(vcfReg$AID %in% vcfProteomy$AID ~ "1", TRUE ~ "0")

# Defining LR model for passing to map function
vcf_model <- function(df){lm(eGFRw.log.Res ~ Dosage, data = df)}

vcfModels <- 
  vcfReg %>%
  select(AID, eGFRw.log.Res, Peptited, starts_with("chr")) %>%
  pivot_longer(cols = starts_with("chr"),
               names_to = c("SNPid"),
               values_to  = c("Dosage")) %>%
  group_by(SNPid, Peptited) %>%
  nest()

  
vcfModels %>%
  mutate(model = data %>% map(vcf_model),
         tidy  = model %>% map(broom::tidy)) %>%
  unnest(tidy) %>%
  filter(term != "(Intercept)") %>% 
  left_join(repSNPs[c("Locus", "SNPid")],
            .,
            by = "SNPid") %>%
  select(Locus, SNPid, Peptited, estimate) %>%
  pivot_wider(id_cols = c(Locus, SNPid),
              names_from  = Peptited,
              values_from = estimate,
              names_prefix = "Peptited_") %>%
  ggplot(aes(x = Peptited_0,
             y = Peptited_1,
             shape = Locus,
             color = Locus)) +
  geom_abline(slope = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 0.9, size = 2.5) +
  scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  scale_color_manual(values = c(
    "maroon1", "darkorchid2","orange2", "green4", "steelblue2", "darkturquoise", "tomato", 
             "springgreen2", "royalblue2", "gold", "grey50"))+
  labs(x = "SNPs effect on ln(eGFRcreat) for subset with proteins",
       y = "SNPs effect on ln(eGFRcreat) for subset without proteins")+
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.placement = "outside",
        axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=7, face="bold"),
        axis.title = element_text(size=14, face="bold"),
        legend.key.size  = unit(0.99, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))
  
ggsave("26-Jan-23_SNPs effect for subset with_without proteins.png", last_plot(), width = 10, height = 7, pointsize = 4, dpi = 300, units = "in")

                                          

# nest(data = -Peptited) %>%
# mutate(model   = map(.x = data, ~ lm(eGFRw.log.Res ~ `chr5:177386403`, data = .x)),
#        glance  = model %>% map(broom::glance),
#        rsq     = glance%>% map_dbl("r.squared"),
#        tidy    = model %>% map(broom::tidy),
#        augment = model %>% map(broom::augment)) %>%
# unnest(glance, tidy)


vcfModels %>%
  mutate(
         glance  = vcfModels %>% map(broom::glance),
         rsq     = glance    %>% map_dbl("r.squared"),
         tidy    = vcfModels %>% map(broom::tidy),
         augment = vcfModels %>% map(broom::augment)
         )


# Checking the distribution
vcfReg[c("AID", "eGFRw.log.Res", "Peptited", "chr5:177386403")] %>%
  mutate(Dosage_SLC34A1 = cut(`chr5:177386403`,
                              breaks = c(-Inf, 0.500, 1.500, Inf),
                              labels = c("0", "1", "2"))) %>%
  #ggplot(aes(Dosage_SLC34A1, eGFRw.log.Res))+
  #geom_violin(aes(fill = Peptited), position = "dodge")+
  #geom_boxplot(aes(fill = Peptited), width = 0.3, position = position_dodge(.9))+
  ggplot(aes(x = `chr5:177386403`)) +
  geom_density(aes(fill = Peptited), alpha = 0.6, color = "grey50") + #, show.legend = NULL
  #scale_fill_manual(values = c("violetred1", "turquoise2"), labels = c("No", "Yes")) +
  theme_classic()


ggsave("25-Jan-23_Density plot eGFR vs SLC34A1 grouped by Peptited2.png", last_plot(), width = 10, height = 7, pointsize = 4, dpi = 300, units = "in")

#-----------#
#Step 2: trait as covariate

PMA_Step2_raw <-
  step2_to_table(vcfProteomy[proteins],
                 vcfProteomy[targets],
                 vcfProteomy)

#Turning the results to longer format for merging with step3
PMA_results_Step2_long <- 
  PMA_Step2_raw %>% 
  pivot_longer(cols = -c(SNPid, Locus),
               names_to = c("Trait", "value"),
               names_pattern = "(.+).(Estimate|SE|Pvalue)$",
               values_to = c("score")) %>%
  pivot_wider(names_from = "value",
              values_from = "score") %>% 
  group_by(Locus, SNPid) %>%
  mutate(outlier = ifelse(is_outlier(Estimate),
                          "Yes",
                          "No")) %>% 
  ungroup()

write.csv(PMA_results_Step2_long,
          "24-Jan-23_PMA_Step2_SNPs adjusted for proteins_long format.csv",
          row.names = F,
          quote = F)

#-----------#
#Step 3: protein as outcome

PMA_Step3_raw <-
  step3_to_table(vcfProteomy[proteins],
                 vcfProteomy[targets],
                 paste("trait ~ SNP + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
                 vcfProteomy)

#Appending Locus column to the results
PMA_results_Step3_long <-
  repSNPs %>%                    
  select(SNPid, Locus) %>% 
  inner_join(PMA_Step3_raw,
             by = c("SNPid")) %>% 
  #change the order of the columns
  select(contains(c("SNPid",
                    "Locus",
                    "pheno",
                    "Estimate",
                    "SE",
                    "Pvalue"))) %>%
  #Significant association between variants and metabolites
  mutate(associated = ifelse(Pvalue <= 0.05/1480,
                             "Yes",
                             "No"))

write.csv(results_Step3_long, "29-Nov-2022_MMD_Step3_SNPs associated with metabolites.csv", row.names = F, quote = F)

#-----------------------------------------------------#
# ----------------- Merge 3 Steps --------------------
#-----------------------------------------------------#

PMA_results_Step2_long %>%
  inner_join(PMA_results_Step3_long,
             by = c("SNPid" = "SNPid",
                    "Locus" = "Locus",
                    "Trait" = "pheno"),
             suffix = c("_Step2", "_Step3")) %>%
  inner_join(repSNPs[c("SNPid",
                       "Locus",
                       "RA_CHRIS_disc",
                       "EA_CHRIS_disc",
                       "Beta_CHRIS",
                       "SE_CHRIS",
                       "Pvalue_CHRIS")],
             by = c("SNPid", "Locus")) %>% 
  mutate(outlierAssociated = case_when(
      outlier == "Yes" & associated == "Yes" ~ "Yes/Yes",
      outlier == "Yes" & associated != "Yes" ~ "Yes/No",
      outlier != "Yes" & associated == "Yes" ~ "No/Yes",
      outlier != "Yes" & associated != "Yes" ~ "No/No"),
    mediator = ifelse(outlier == "Yes" & associated == "Yes", "Yes", "No"),
    SNP = factor(SNPid,
                 levels = str_sort(unique(SNPid),
                                   numeric = TRUE,
                                   decreasing = TRUE))) %>% 
  rename(OA = RA_CHRIS_disc, 
         EA = EA_CHRIS_disc,
         Estimate_GWAS = Beta_CHRIS,
         SE_GWAS       = SE_CHRIS,
         Pvalue_GWAS   = Pvalue_CHRIS) %>%
  select(SNPid, Locus, EA, OA, Estimate_GWAS, SE_GWAS, Pvalue_GWAS, everything()) %>% #View()
  #count(Locus, associated) %>% View() #pivot_longer(cols = -Locus) 
  filter(associated == "Yes") %>% View()
  #write.csv(., "29-Nov-22_Heatmap_MMD_Step 1&2&3_outlierRelated traits_Mediatory metabolites.csv", row.names = FALSE)
  ggplot(aes(x = trait, y = SNP, fill = outlierAssociated)) +
  geom_tile() +
  theme_classic() +
  scale_fill_manual(values=c('white', "grey50", "#FF6666"))+ #'#999999', "#E69F00"
  labs(x    = "",
       y    = "")+
  theme(legend.title= element_text(size=7, face="bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text  = element_text(size = 8),
        axis.text.x = element_text(size=4, face="bold", angle=90, vjust=1.05, hjust=1),
        axis.text.y = element_text(size=4, face="bold"),
        axis.title  = element_text(size=12,face="bold"))

ggsave("29-Nov-22_Heatmap_MMD_Step 1&2&3_outlierRelated traits.png", 
       last_plot(), width = 10, height = 9, pointsize = 4, dpi = 300, units = "in")
#-----------#

proteinsInfo %>%
  as_tibble() %>% 
  filter(Protein.Group == "P00748") %>%
    select(First.Protein.Description)

#-----------------------------------------------------#
# ----------------- Linerange Plot -------------------
#-----------------------------------------------------#

#Step2: Line range plot of the entire 163 variants
pdf('24-Jan-23_Linerange plot of protein-adjusted effect of SNPs.pdf', width=8, height = 10)

map(1:length(targets), function(i){
  PMA_results_Step2_long %>%
    filter(SNPid %in% targets[i]) %>%
    # mutate(Trait = recode(Trait,
    #                       "Body_Fat"   = "Body Fat",
    #                       "Visceral_Fat" = "Visceral Fat",
    #                       "Pulse_Rate" = "Pulse Rate",
    #                       "INR_PT"     = "INR PT",
    #                       "APTT_ratio" = "APTT ratio",
    #                       "AST_GOT"    = "AST GPT",
    #                       "ALT_GPT"    = "ALT GPT",
    #                       "Urine_pH"   = "Urine pH"),
    mutate(Locus_reordered = Locus_factor(Locus),
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
    labs(x = paste0("Protein-adjusted effect of the ", targets[i], " SNP in *", repSNPs[i,"Locus"], "* locus on ln(eGFRcreat)")) + # in SHROOM3 Locus
    theme_classic() +
    theme(strip.background = element_blank(), strip.placement = "outside",
          axis.text.x = element_text(size = 6,  face = "bold"),
          axis.text.y = element_text(size = 5,  face = "bold"),
          axis.title  = element_text(size = 10, face = "bold"),
          axis.title.x = ggtext::element_markdown(),
          axis.title.y.left = element_blank(),
          panel.grid.major.y = element_line(color = "lightgray", size = .25),
          strip.text.x = element_text(size = 8, color = "Black", face = "bold.italic"))
}
)

dev.off()

#-----------------------------------------------------#

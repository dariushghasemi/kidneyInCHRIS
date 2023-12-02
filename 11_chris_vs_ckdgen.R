#=========================================#
#  Comparing Replication with Discovery
#=========================================#

# Started on October 2021
# Last update on January 2023 

#---------- template result -----------

library(data.table)
library(readxl)
library(ggplot2)
library(stringr)
library(tidyverse)
#-----------------------------------------------------#

# Reading the initial replication analysis results on ln(eGFRcrea)
tempResult0 <- read_excel("F:/Dariush/PhD/Analysis/1st Paper/Outputs_tables/Replication_eGFRw.log.Res_MultiA/Template_results_211012_eGFR.log.std.xlsx", sheet = 2, na = "NA")#skip = 1, [-1, ]

# Reading the first replication analysis results on the residuals of ln(eGFRcrea) made by indirect method
tempResult1 <- read.table("F:/Dariush/PhD/Analysis/1st Paper/Outputs_tables/Replication_eGFRw.log.Res_MultiA/Meta_replicatedSNPs_eGFRw.log.Res.all.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# adding the locus name to replication results 
tempResult  <- merge(tempResult0[,c("CHR", "BEG", "Locus")],
                     tempResult1, by = c("CHR", "BEG"), all.y = TRUE)

tempResult1 <- as.data.frame(tempResult1)
tempResult$Locus <- as.factor(tempResult$Locus)

# missing SNPs
mergedSupTable <- read.delim("E://Dariush//PhD//Analysis//mergedtable_hg19_SerumCreatLog.txt", header= TRUE, sep = "\t")

missedSNPs <- setdiff(Supptable$BEG, mergedSupTable$BEG)
Supptable[Supptable$BEG == missedSNPs, c("chr","BEG")]

#-----------------------------------------------------#
# 147 CKDGen Loci
repSNPs_EA <- read.csv("F:/Dariush/PhD/Analysis//1st Paper/Outputs_tables/Replication_eGFRw.log.Res_EUR_A/31-Dec-22_Replication_of_147_EA_Loci_in_CHRIS.txt", sep = '\t')

# replication analysis results: 162 SNPs (Multi-Ancestry)
repSNPs_old <- read.csv("F:/Dariush/PhD/Analysis//1st Paper/Outputs_tables/ReplicatedSNPs by eGFRw.log.Res.csv",
                        stringsAsFactors = FALSE)

# replication analysis results: 163 SNPs (European Ancestry) - including GWAS hit at PDILT locus
repSNPs_tmp <- read.table("F:/Dariush/PhD/Analysis//1st Paper/Outputs_tables/Replication_eGFRw.log.Res_EUR_A/10-Jan-23_replicatedSNPs_eGFRw.log.Res.txt",
                          header = T, stringsAsFactors = FALSE)

#-----------------------------------------------------#
#CHRIS
#str_extract(tempResult$MARKER_ID.log.Std, "[A-Z]:[A-Z]")

#------- changing the signs of positive alleles to negative effect -------
propVars <- c("Effect.ckdgen",
              "Effect.CHRIS",
              "CHRIS.CKDGen.Effect.Ratio")

#-----------------------------------------------------#
# ------------- Concordant Effect Allele -------------
#-----------------------------------------------------#

# Making concordant effect allele and changing the sign
# of the positive effect allele in both studies to have
# negative effect for consistency in one single operation

repSNPs <- 
  repSNPs_tmp %>%
  rename(EA_CKDGen_disc  = EA_CKDGen,
         RA_CKDGen_disc  = RA_CKDGen,
         EAF_CKDGen_disc = EAF_CKDGen,
         MAF_CHRIS_disc  = MAF_CHRIS) %>% #Pay attention to SNP order in SLC34A1, DAB2 Loci
         #Locus           = Closest.Gene) %>% #Closest.Gene was manually added later 
  mutate(SNPid           = str_c("chr", CHR, ":", POS38),  #CHR:POS
         RSID            = case_when(POS38 == 10603539 ~ "rs1641044855",
                                     POS38 == 78319103 ~ "rs200047230",
                                     POS38 == 78323023 ~ "rs34351881",
                                     POS38 == 78331260 ~ "rs202186505",
                                     POS38 == 78336491 ~ "rs114142124",
                                     POS38 == 78362081 ~ "rs1376948391",
                                     POS38 == 78381674 ~ "rs369377897",
                                     POS38 == 98741048 ~ "rs199654172", TRUE ~ RSID),
         RA_CHRIS_disc   = as.character(noquote(str_split(MARKER_ID_CHRIS, "\\W", simplify = T)[,5])),  #Reference allele
         EA_CHRIS_disc   = as.character(noquote(str_split(MARKER_ID_CHRIS, "\\W", simplify = T)[,6])),  #Effect allele
         EAF_CHRIS_disc  = round(AC_CHRIS / (NS_CHRIS * 2), 5),
         #size direction when CKDGen EA is discordant to CHRIS
         #Cristian noticed when the CKDGen Effect Allele Frequency is > 0.5, 
         #the CKDGen Effect Allele is always discordant from the CHRIS Effect Allele
         Beta_CHRIS_tmp  = ifelse(EA_CKDGen_disc != EA_CHRIS_disc, (-1) * Beta_CHRIS, Beta_CHRIS),
         EA_CHRIS_tmp    = ifelse(EA_CKDGen_disc != EA_CHRIS_disc, RA_CHRIS_disc, EA_CHRIS_disc),
         RA_CHRIS_tmp    = ifelse(EA_CKDGen_disc != EA_CHRIS_disc, EA_CHRIS_disc, RA_CHRIS_disc), #Changing CHRIS effect
         EAF_CHRIS_tmp   = ifelse(EA_CKDGen_disc != EA_CHRIS_disc, 1 - EAF_CHRIS_disc, EAF_CHRIS_disc),
         Check_EA        = ifelse(EA_CKDGen_disc == EA_CHRIS_tmp, "Yes", "No"),
         Check_Beta      = case_when(Beta_CKDGen > 0 & Beta_CHRIS_tmp > 0 ~ "Yes_Pos",
                                     Beta_CKDGen < 0 & Beta_CHRIS_tmp < 0 ~ "Yes_Neg"),
         #To align the direction of associations in both studies, we need to:
         Beta_CKDGen_ald = ifelse(Beta_CKDGen > 0, (-1) * Beta_CKDGen, Beta_CKDGen),
         Beta_CHRIS_ald  = ifelse(Beta_CKDGen > 0, (-1) * Beta_CHRIS_tmp, Beta_CHRIS_tmp),
         Check_Beta_ald  = ifelse(Beta_CKDGen_ald < 0 && Beta_CHRIS_ald < 0, "Yes", "No"),
         # Check_Beta      = case_when(Beta_CKDGen > 0 & Beta_CHRIS_tmp > 0 ~ "Yes_Pos",
         #                             Beta_CKDGen < 0 & Beta_CHRIS_tmp < 0 ~ "Yes_Neg"),
         EA              = ifelse(Beta_CKDGen > 0, RA_CHRIS_tmp, EA_CHRIS_tmp),
         OA              = ifelse(Beta_CKDGen > 0, EA_CHRIS_tmp, RA_CHRIS_tmp),
         EAF_CKDGen      = ifelse(Beta_CKDGen > 0, 1 - EAF_CKDGen_disc, EAF_CKDGen_disc),
         EAF_CHRIS       = ifelse(Beta_CKDGen > 0, 1 - EAF_CHRIS_tmp, EAF_CHRIS_tmp),
         #RAF_CKDGen     = 1 - EAF_CKDGen,         
         MAF_CKDGen      = ifelse(EAF_CKDGen < .5, EAF_CKDGen, 1 - EAF_CKDGen),
         MAF_CHRIS       = ifelse(EAF_CHRIS  < .5, EAF_CHRIS,  1 - EAF_CHRIS),
         MAF_Diff        = MAF_CHRIS - MAF_CKDGen,
         MAF_Ratio       = MAF_CHRIS / MAF_CKDGen,
         Beta_Diff       = Beta_CHRIS_ald - Beta_CKDGen_ald,
         Beta_Ratio      = Beta_CHRIS_ald / Beta_CKDGen_ald) #%>% View
  # select(MARKER_ID_CHRIS, 
  #        Beta_CHRIS, Beta_CKDGen, Beta_CKDGen_ald, Beta_CHRIS_ald, Check_Beta,
  #        EAF_CKDGen_disc, EAF_CKDGen, EAF_CHRIS_disc, EAF_CHRIS,
  #        RA_CKDGen_disc, EA_CKDGen_disc,
  #        RA_CHRIS_tmp, EA_CHRIS_tmp, Check_EA, EA, OA)  %>% View #%>% count(Check_Beta)

#drop_na(any_of(c("Beta_CHRIS_ult", "Beta_CKDGen_ult"))

#------------#
# Supplementary table 3: 163 replicated SNPs in CKDGen and CHRIS

repSNPs %>%
  mutate(SNPid = str_remove(SNPid, "chr"), #below I tried to fill in properly the missing EA/OA, EAF, beta, and SE
         EA_OA = case_when(
           is.na(Beta_CHRIS_ald) & Beta_CHRIS < 0 ~ paste0(EA_CHRIS_disc, "/", RA_CHRIS_disc),
           is.na(Beta_CHRIS_ald) & Beta_CHRIS > 0 ~ paste0(RA_CHRIS_disc, "/", EA_CHRIS_disc),
           TRUE ~ paste0(EA, "/", OA)),
         EAF_CHRIS = case_when(
           is.na(Beta_CHRIS_ald) & Beta_CHRIS < 0 ~ EAF_CHRIS_disc,
           is.na(Beta_CHRIS_ald) & Beta_CHRIS > 0 ~ 1 - EAF_CHRIS_disc,
           TRUE ~ EAF_CHRIS),
         Beta_CHRIS_ald = case_when(
           is.na(Beta_CHRIS_ald) & Beta_CHRIS < 0 ~ Beta_CHRIS,
           is.na(Beta_CHRIS_ald) & Beta_CHRIS > 0 ~ Beta_CHRIS * (-1),
           TRUE ~ Beta_CHRIS_ald),
         Pvalue_CHRIS_1s = Pvalue_CHRIS/2) %>%
  select(Locus, RSID, SNPid, EA_OA,
         EAF_CKDGen, Beta_CKDGen_ald, SE_CKDGen, Pvalue_CKDGen,
         EAF_CHRIS,  Beta_CHRIS_ald,  SE_CHRIS,  Pvalue_CHRIS,  Pvalue_CHRIS_1s) %>% #View()
  write.csv("19-Jan-23_Suppl. Table 3_163 replicated SNPs in CHRIS.csv", row.names = F, quote = F)
  
  
  

#------------#
# Supplementary table 2: 147 Loci in CKDGen and CHRIS

repSNPs_EA %>%
  rename(EA_CKDGen_disc  = EA_CKDGen,
         RA_CKDGen_disc  = RA_CKDGen,
         #EAF_CKDGen_disc = EAF_CKDGen,
         MAF_CHRIS_disc  = MAF_CHRIS,
         Locus           = Closest.Gene) %>%
  mutate(SNPid           = str_c(CHR, ":", POS38),                                                      #CHR:POS
         RA_CHRIS_disc   = as.character(noquote(str_split(MARKER_ID_CHRIS, "\\W", simplify = T)[,5])),  #Reference allele
         EA_CHRIS_disc   = as.character(noquote(str_split(MARKER_ID_CHRIS, "\\W", simplify = T)[,6])),  #Effect allele
         EAF_CHRIS_disc  = round(AC_CHRIS / (NS_CHRIS * 2), 5),                                         #Computing CHRIS effect allele frequency
         Beta_CHRIS_tmp  = ifelse(EA_CKDGen_disc != EA_CHRIS_disc, (-1) * Beta_CHRIS, Beta_CHRIS),      #Changing CHRIS effect size
         EA_CHRIS_tmp    = ifelse(EA_CKDGen_disc != EA_CHRIS_disc, RA_CHRIS_disc, EA_CHRIS_disc),       #Changing CHRIS reference allele
         RA_CHRIS_tmp    = ifelse(EA_CKDGen_disc != EA_CHRIS_disc, EA_CHRIS_disc, RA_CHRIS_disc),       #Changing CHRIS effect allele
         EAF_CHRIS       = ifelse(EA_CKDGen_disc != EA_CHRIS_disc, 1 - EAF_CHRIS_disc, EAF_CHRIS_disc), #Correcting EAF
         Check_EA        = ifelse(EA_CKDGen_disc == EA_CHRIS_tmp, "Yes", "No"),                         #Check if the effect alleles are concordant
         EA              = ifelse(EA_CKDGen_disc == EA_CHRIS_tmp, EA_CKDGen_disc, NA),                  #Extracting the concordant effect allele
         OA              = ifelse(RA_CKDGen_disc == RA_CHRIS_tmp, RA_CKDGen_disc, NA),                  #Extracting the concordant reference allele
         Check_Beta      = case_when(Beta_CKDGen > 0 & Beta_CHRIS_tmp > 0 ~ "Yes_Pos",                  #Check if association signs are in one direction
                                     Beta_CKDGen < 0 & Beta_CHRIS_tmp < 0 ~ "Yes_Neg"),
         Pvalue_CHRIS_1s = case_when(Check_Beta == "Yes_Pos" | Check_Beta == "Yes_Neg" ~ Pvalue_CHRIS/2,#Compute 1-sided P-value for concordant alleles
                                     is.na(Check_Beta) ~ 1-(Pvalue_CHRIS/2)),                           #Compute complementary 1-sided P-value for discordant alleles
         EA_OA           = paste0(EA_CKDGen_disc, "/", RA_CKDGen_disc),                                 #Compute consistent EA/OA in bothe studies
         Beta_SE_CKDGen  = paste0(Beta_CKDGen, '(', round(SE_CKDGen,6), ')'),                           #Compute b(SE) in CKDGen
         Beta_SE_CHRIS   = paste0(Beta_CHRIS_tmp,  '(', round(SE_CHRIS, 6), ')')) %>%                   #Compute b(SE) in CHRIS
  select(Locus, RSID, SNPid, EA_OA,
         #EA_CKDGen_disc, RA_CKDGen_disc, EA_CHRIS_tmp, RA_CHRIS_tmp, Check_EA, EA, OA, EA_OA,
         EAF_CKDGen, Beta_CKDGen,    SE_CKDGen, Pvalue_CKDGen,
         EAF_CHRIS,  Beta_CHRIS_tmp, SE_CHRIS,  Pvalue_CHRIS, Pvalue_CHRIS_1s) %>% #View
  write.csv("19-Jan-23_Suppl. Table 2_147 CKDGen Loci.csv", row.names = FALSE)

# Alt+- is a shortcut for <-
#------------#
# Subsetting the leading or most significant SNPs

tagSNPs <- 
  repSNPs %>%
  select(SNPid, Locus, RSID, 
         Beta_CKDGen_ald, Pvalue_CKDGen,
         Beta_CHRIS_ald,  Pvalue_CHRIS) %>%
  group_by(Locus) %>% 
  top_n(n = -1, wt = Pvalue_CKDGen) #%>% View

#-------------#
# Supplementary tables of the paper: full list of 162 replicated SNPs

repSNPs %>%
  filter(SNPid %in% tagSNPs$SNPid) %>%
  #filter(SNPid %in% CKDGenSNPs) %>%
  #filter(SNPid %in% leadingSNPs) %>%
  #filter(Pvalue_CHRIS <= 0.00068) %>%     #Checking UMOD in all 147 CKDGen Loci
  mutate(EA_OA          = paste0(EA, "/", OA),
         Beta_SE_CKDGen = paste0(Beta_CKDGen_ald, '(', round(SE_CKDGen,6), ')'),
         Beta_SE_CHRIS  = paste0(Beta_CHRIS_ald,  '(', round(SE_CHRIS, 6), ')'),
         Pvalue_CHRIS_1s = Pvalue_CHRIS/2) %>%
  select(SNPid, Locus, RSID, EA_OA,
         EAF_CKDGen, Beta_SE_CKDGen, Pvalue_CKDGen,
         EAF_CHRIS,  Beta_SE_CHRIS,  Pvalue_CHRIS, Pvalue_CHRIS_1s) %>%
  # Find the smallest and the largest p-values among the loci
  #top_n(n = -1, wt = Pvalue_CHRIS_1s) %>% View
  write.csv("10-Jan-23_Table 2_Leading SNPs of the replicated loci.csv", row.names = FALSE) #"10-Jan-23_Supp Table 3_162 Replicated SNPs.csv"
  
#-----------------------------------------------------#
#-------------- Testing the correlations -------------
#-----------------------------------------------------#

cor.test(repSNPs$MAF.Diff,
         repSNPs$Effect.CHRIS - repSNPs$Effect.ckdgen,
         use = "complete.obs")

cor(repSNPs[c("MAF.Diff", "Effect.CHRIS", "Effect.ckdgen")], 
    use = "complete.obs")

#Visualizing the correlation matrix
library(ggcorrplot)

#Select MAF and effect size columns
cor_data <- 
  repSNPs %>% 
  select(#Locus,
    MAF_CHRIS, MAF_CKDGen, MAF_Diff, MAF_Ratio,
    Beta_CHRIS_ald, Beta_CKDGen_ald, Beta_Diff, Beta_Ratio) %>% 
  drop_na(MAF_CKDGen)

#Compute correlation Matrix
corr <- round(cor(cor_data, use = "complete.obs"), 2)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(cor_data)

# Add correlation significance level
ggcorrplot(corr, hc.order = TRUE, type = "lower", p.mat = p.mat)

# Add correlation coefficients
ggcorrplot(corr,
           hc.order = TRUE,
           type = "lower", # get the lower triangle
           outline.col = "white",
           ggtheme = ggplot2::theme_bw, # change theme
           colors = c("#6D9EC1", "white", "#E46726"), # change color palette
           p.mat = p.mat, #Barring the non-significant coefficient
           lab = TRUE)# Add correlation coefficients

ggsave("10-Jan-23_Correlation matrix for MAF_Ratio and Effect_Ratio in CHRIS and CKDGen.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#-------------#
# Quantify Effect difference by MAF.diff

summary(lm(Beta_Diff ~ MAF_Diff, data = repSNPs))

# Model the interaction id different levels of TSH
repSNPs %>%
  select(Locus, SNPid, Beta_Diff, MAF_Diff) %>%
  group_by(Locus) %>%
  nest() %>%
  mutate(model = data %>% map(~lm(Beta_Diff ~ MAF_Diff, data = .)),
         tidy  = model %>% map(broom::tidy)) %>%
  unnest(tidy) %>%
  filter(term != "(Intercept)") %>% 
  select(-data, -model) %>%
  write.csv("31-Jan-23_Quantifying MAF difference by regressing on effect difference.csv", 
            quote = F, row.names = F)



#-----------------------------------------------------#
# ----------- Figure 3: MAF vs Effect ratio ----------
#-----------------------------------------------------#

library(ggfittext)

#-------------#
# Beta ratio vs. MAF ratio in CKDGen and CHRIS 
repSNPs %>% 
  ggplot(aes(x = MAF_Ratio, y = Beta_Ratio, shape = Locus, color = Locus)) + 
  geom_point(alpha = 0.9, size = 3.5) +
  #ggrepel::geom_text_repel(data = tagSNPs, aes(label = Locus), check_overlap = TRUE) +
  #ggfittext::geom_fit_text(grow = TRUE)+
  #geom_text(aes(label = Locus), check_overlap = TRUE) +
  #geom_abline(intercept = 0, slope = +1, color = "grey50", size = .9)+
  scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  scale_color_manual(values= c("maroon1", "darkorchid2", "orange2", "green4", "steelblue2", "darkturquoise",
                               "tomato", "springgreen2", "royalblue2", "gold", "grey50")) +
  #scale_color_brewer(palette = "Dark2") +
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(breaks = seq(0,6, 0.5)) +
  scale_x_continuous(breaks = seq(0.75,1.20, .05), limits = c(.75,1.20), expand = c(0,0)) +
  annotate("segment", x = 0.755, xend = 0.995, y = 6, yend = 6, colour = "grey50", size = 1.2,
           arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm"))) +
  annotate("segment", x = 1.005, xend = 1.195, y = 6, yend = 6, colour = "grey50", size = 1.2,
           arrow = arrow(ends = "both", angle = 90, length = unit(0.2, "cm"))) +
  annotate("text", x = 0.85, y = 6.3, size = 4,
           label = "Rarer alleles in CHRIS") + #\nLarger effect in CHRIS
  annotate("text", x = 1.09, y = 6.3, size = 4, 
           label = "Rarer alleles in CKDGen") +#\nLarger effect in CHRIS
  # annotate("text", x = 0.23, y = 5.6, size = 4, label = "CASZ1") +
  # annotate("text", x = 0.7, y = 5.1, size = 4, label = "CASZ1") +
  # annotate("text", x = 0.59, y = 4, size = 4, label = "GAB2") +
  # annotate("text", x = 0.65, y = 3, size = 4, label = "DDX1") +
  # annotate("text", x = 1.15, y = 4, size = 4, label = "IGF1R") +
  # annotate("text", x = 1.2, y = 2.8, size = 4, label = "PIPK1B") +
  labs(x = "CHRIS-to-CKDGen MAF ratio", 
       y = "CHRIS-to-CKDGen effect ratio") + 
  #theme_minimal()+ #theme_light(base_size = 10)+ #lims(y = c(1, 6)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
        panel.grid.major.x = element_blank(),
        axis.text  = element_text(size = 8,  face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        #legend.position = "none",
        legend.key.size  = unit(0.99, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.text  = element_text(size = 12, face = "italic"),
        legend.title = element_text(size = 13, face = "bold"))

ggsave("22-Mar-23_MAF Ratio vs Effect Size Ratio in CHRIS and CKDGen.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")



#-----------------------------------------------------#
# ----------------- Figure 2: Line range -------------
#-----------------------------------------------------#

# The difference of Effect size in CHRIS and CKDGen

repSNPs %>%
  filter(SNPid %in% tagSNPs$SNPid) %>% #%in% leadingSNPs
  select(SNPid, Locus, Beta_CHRIS_ald, Beta_CKDGen_ald, SE_CHRIS, SE_CKDGen) %>% 
  rename(Beta_CHRIS  = Beta_CHRIS_ald, 
         Beta_CKDGen = Beta_CKDGen_ald) %>% 
  pivot_longer(cols = !c(SNPid, Locus),   
               names_to = c("trait", "study"),
               names_pattern = "(.+)_(CHRIS|CKDGen)$", #([A-Za-z0-9_]+)
               values_to = c("score")) %>%
  pivot_wider(names_from = "trait",
              values_from = "score") %>% #View()
  ggplot(aes(x = Locus,
             y = Beta,
             color = study,
             ymin = Beta - 1.95*SE,
             ymax = Beta + 1.95*SE)) + 
  geom_pointrange(aes(color = study),
                  size = 1.1, fatten = 2.5, alpha = 0.9,
                  #lineend = "round",
                  position = position_dodge(.3)) +
  #scale_color_grey(start=0.55, end=0.25) +
  scale_color_manual(values=c('turquoise3','tomato3'))+
  #scale_size_manual(values=c(2, 4))+
  geom_hline(yintercept = 0, color="black", lty = 1, size = .85) +
  scale_x_discrete(limits = tagSNPs$Locus,
                   breaks = tagSNPs$Locus,
                   labels = paste0(tagSNPs$Locus,
                                   "\n",
                                   tagSNPs$RSID)) +
  scale_y_continuous(breaks = seq(-.035, .001, .005)) +
  #theme_light(base_size = 10) +
  labs(y = "Effect on ln(eGFRcrea)", x = NULL) + # x = "Locus"
  coord_cartesian(ylim = c(-0.035, 0.001)) +
  theme(legend.title = element_blank(),
        legend.position = c(.9, .365),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
        axis.text.y  = element_text(size = 8,  face = "bold"),
        axis.text.x  = element_text(size = 8,  face = "bold.italic"),
        axis.title   = element_text(size = 14, face = "bold"),
        legend.key.size = unit(0.6, 'cm'),
        legend.text  = element_text(size = 12))

ggsave("22-Mar-23_Comparison of leading SNPs Effect size_wot Locus.png",
       last_plot(), width = 8.5, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-----------------------------------------------------#

#filling the locus information missed for some of the replicated variants

tempResult <- tempResult %>%
  mutate(Locus = replace(Locus, which(CHR == '5' & BEG == '39385539' ), "DAB2"))%>%
  mutate(Locus = replace(Locus, which(CHR == '5' & BEG == '39393631' ), "DAB2"))%>%
  mutate(Locus = replace(Locus, which(CHR == '5' & BEG == '39397030' ), "DAB2"))%>%
  mutate(Locus = replace(Locus, which(CHR == '7' & BEG == '77714744' ), "RSBN1L"))%>%
  mutate(Locus = replace(Locus, which(CHR == '7' & BEG == '77736048' ), "RSBN1L"))%>%
  mutate(Locus = replace(Locus, which(CHR == '8' & BEG == '23885208' ), "STC1"))
#-----------------------------------------------------#

p1 <- 
  repSNPs %>% 
  #tempResult %>% 
  #drop_na(any_of(propVars))
  drop_na(any_of(c("Effect.ckdgen",
                   "Effect.CHRIS",
                   "CHRIS.CKDGen.Effect.Ratio"))) %>%
  ggplot(aes(x = Effect.ckdgen,
             y = Effect.CHRIS)) + 
  geom_point(aes(shape = Locus,
                 color = Locus),
             size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) +
  #scale_color_manual(values=c( ' #999999 ' , ' #E69F00 ' , ' #56B4E9 ' ))+
  #scale_size_manual(values=c(1.5, 2, 3))+
  #theme(legend.position="top")
  #geom_quantile()+
  #xlim(-0.0074, -0.002)+
  geom_abline(intercept = 0,
              slope = +1,
              color = "grey50",
              size = 1.5,
              lty = 2)+
  geom_vline(xintercept=0)+#, linetype="dashed", color = "darkgray", size=1.5)+
  geom_hline(yintercept=0)+#, linetype="dashed", color = "darkgray", size=1.5)+
  xlab("Effect CKDGen")+#"Association coefficient for log(eGFRcrea) in CKDGen"
  ylab("Effect CHRIS")+#"Association coefficient for log(eGFRcrea) in CHRIS"
  theme_classic()+
  theme(axis.text  = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.key.size = unit(0.65, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16))

ggsave("effectDiff_eGFRw_slides.png",last_plot(), width = 8, height = 5.5, 
       pointsize = 5, dpi = 300, units = "in")
#-------------#

p2 <- 
  repSNPs %>% 
  #tempResult %>% 
  #drop_na(any_of(propVars))
  drop_na(any_of(c("Effect.ckdgen",
                   "Effect.CHRIS",
                   "CHRIS.CKDGen.Effect.Ratio"))) %>%
  ggplot(aes(x = Effect.ckdgen,
             y = MAF.Diff)) +
  geom_point(aes(color = Locus),
             position  = position_jitter(height = 0L, seed = 1L)) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) + 
  geom_linerange(aes(x = Effect.ckdgen,
                     ymax = MAF.Diff,
                     ymin = 0,
                     color= Locus),
                 position = position_jitter(height = 0L, seed = 1L))+
  xlab("Effect CKDGen") + #"Association coefficient for log(eGFRcrea) in CKDGen"
  ylab("MAF CHRIS - MAF CKDGen") + #"Diffrence of Alleles frequency between CHRIS and CKDGen"
  theme_classic()+
  theme(axis.text  = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.key.size = unit(0.65, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16))
#theme(panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
#      axis.text.x = element_text(#face="bold",color="steelblue",
#  size=8, angle=90))#+
#theme(axis.title.x = element_text(size=12))#,face="italic"
#-------------#

library(gridExtra)
p3 <- grid.arrange(p1, p2, ncol=1)

ggsave("21-Apr-22_CKDGen Betas versus MAF.Diff & CHRIS Betas_4Poster.png", 
       p3, width = 10, height = 6.5, pointsize = 8, dpi = 600, units = "in")
#-----------------------------------------------------#
library(reshape2)
#melt((tempResult[, c("Locus", "MAF.ckdgen", "MAF.log.Std" ,"MAF.Diff")]), id.vars = )
tempResult2 <- data.frame("Study"=rep(c("CKDGen", "CHRIS"), each=nrow(tempResult)), 
                          "MAF"=rbind(as.matrix(tempResult$MAF.ckdgen), as.matrix(tempResult$MAF)),
                          "MAF.Diff"=rep(tempResult$MAF.Diff, each=2))
#-------------#
ggplot(tempResult2, aes(x=MAF.Diff, y=MAF))+ geom_point(aes(color=Study), alpha=0.6, size=3)+
  geom_vline(xintercept=0,linetype="dashed") + geom_hline(yintercept=0, linetype="dashed")+
  theme(panel.background = element_blank(), panel.border = element_blank())
#-------------#

dim(data.frame("Study"=rep(c("CKDGen", "CHRIS"), each=nrow(tempResult))))
#-----------------------------------------------------#
library(tidyverse)

tempResult %>% 
  drop_na(any_of(propVars)) %>%
  ggplot(aes(x=Effect.ckdgen, y=CHRIS.CKDGen.Effect.Ratio))+ 
  geom_point(aes(color=Locus, shape=Locus), size=2)+
  geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
  #geom_text(aes(label= Locus), size = 1.5)+
  scale_shape_manual(values=c(0, 1, 2,3, 7, 8, 9, 15, 16, 17))+ theme_bw()#, 15, 16, 17

ggsave("CKDGen Beta vs CHRIStoCKDGen.Effects Ratio.png",
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-------------#
tempResult %>% 
  drop_na(CHRIS.CKDGen.Effect.Ratio) %>%
  ggplot(aes(MAF.Diff, CHRIS.CKDGen.Effect.Ratio))+ 
  geom_point(aes(shape=Locus, color=Locus), size=2)+
  geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
  scale_shape_manual(values=c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17))+ theme_bw()

ggsave("MAF.Diff vs CHRIStoCKDGen.Effects Ratio.png",
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

write.csv(tempResult,"ReplicatedSNPs by eGFRw.log.Res.csv")

#-----------------------------------------------------#
# ----------- Multi-Anc vs. European-Anc -------------
#-----------------------------------------------------#

#Comapring CKDGen Multi-Ancestry Meta-GWAS vs. European Ancestry

repSNPs %>%
  mutate(SNPid = as.character(noquote(str_split(MARKER_ID_CHRIS,
                                                "_",
                                                simplify=TRUE)[,1]))) %>%
  #filter(Locus == "CASZ1") %>% View
  select(SNPid, Locus, EAF_CKDGen, Beta_CKDGen, SE_CKDGen, Pvalue_CKDGen) %>%
  inner_join(repSNPs_old[c("SNPid", "Locus", "Freq1", "Effect", "StdErr", "P.value")],
             by = c("SNPid", "Locus")) %>%
  mutate(Beta_to_SE_MultiAncestry = Effect / StdErr,
         Beta_to_SE_EuropAncestry = Beta_CKDGen / SE_CKDGen) %>%
  ggplot(aes(x = Beta_to_SE_MultiAncestry, y = Beta_to_SE_EuropAncestry, color = Locus)) +
  #ggplot(aes(x = P.value, y = Pvalue_CKDGen, color = Locus)) +
  #ggplot(aes(x = Freq1, y = EAF_CKDGen, color = Locus)) +
  geom_point() +
  geom_abline(slope = 1) + #geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  theme_classic()

#-------------#
#-------------#


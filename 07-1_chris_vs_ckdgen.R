#=========================================#
#  Comparing Replication with Discovery
#=========================================#

# Started on October 2021
# Last update on January 2023 

#------------#

library(data.table)
library(readxl)
library(ggplot2)
library(stringr)
library(tidyverse)


#-----------------------------------------------------#
#----------        template result         -----------
#-----------------------------------------------------#

# input files directories
base.dir <- "F:/Dariush/PhD/Analysis/1st Paper/Outputs_tables"
egfr_log_std <- paste0(base.dir, "/Replication_eGFRw.log.Res_MultiA/Template_results_211012_eGFR.log.std.xlsx")
egfr_log_res <- paste0(base.dir, "/Replication_eGFRw.log.Res_MultiA/Meta_replicatedSNPs_eGFRw.log.Res.all.txt")
egfr_log_res_MA <- paste0(base.dir, "/Replication_eGFRw.log.Res_MultiA/ReplicatedSNPs by eGFRw.log.Res.csv")
egfr_log_res_EA <- paste0(base.dir, "/Replication_eGFRw.log.Res_EUR_A/10-Jan-23_replicatedSNPs_eGFRw.log.Res.txt")
ckdgen_EA    <- paste0(base.dir, "/Replication_eGFRw.log.Res_EUR_A/31-Dec-22_Replication_of_147_EA_Loci_in_CHRIS.txt")
ckdgen_MA    <- paste0(base.dir, "/Replication_eGFRw.log.Res_MultiA/07-Dec-23_Replication_of_147_Trans-A_loci_in_CHRIS.txt")

#------------#
# Reading the initial replication analysis results on ln(eGFRcrea)
tempResult0 <- read_excel(egfr_log_std, sheet = 2, na = "NA")#skip = 1, [-1, ]

# Reading the first replication analysis results on the residuals of ln(eGFRcrea) made by indirect method
tempResult1 <- read.table(egfr_log_res, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# adding the locus name to replication results 
tempResult  <- merge(tempResult0[,c("CHR", "BEG", "Locus")],
                     tempResult1, by = c("CHR", "BEG"), all.y = TRUE)

tempResult$Locus <- as.factor(tempResult$Locus)
tempResult1 <- as.data.frame(tempResult1)

# missing SNPs
mergedSupTable <- read.delim("E://Dariush//PhD//Analysis//mergedtable_hg19_SerumCreatLog.txt", header= TRUE, sep = "\t")

missedSNPs <- setdiff(Supptable$BEG, mergedSupTable$BEG)
Supptable[Supptable$BEG == missedSNPs, c("chr","BEG")]

#------------#
# 147 CKDGen Loci European Ancestry
repSNPs_EA <- read.csv(ckdgen_EA, sep = '\t')

# 147 CKDGen Loci Multi Ancestries
repSNPs_MA <- read.delim(ckdgen_MA, sep = '\t')

# replication analysis results: 162 SNPs (Multi-Ancestry)
repSNPs_old <- read.csv(egfr_log_res_MA, stringsAsFactors = FALSE)

# replication analysis results: 163 SNPs (European Ancestry) - including GWAS hit at PDILT locus
repSNPs_tmp <- read.table(egfr_log_res_EA, header = T, stringsAsFactors = FALSE)

#------------#
# CHRIS
#str_extract(tempResult$MARKER_ID.log.Std, "[A-Z]:[A-Z]")


#-----------------------------------------------------#
#-------      Concordant Effect Allele        --------
#-----------------------------------------------------#

# Making concordant effect allele and changing the sign
# of the positive effect allele in both studies to have
# negative effect for consistency in one single operation

repSNPs <- repSNPs_tmp %>%
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
         EAF_CHRIS,  Beta_CHRIS_ald,  SE_CHRIS,  Pvalue_CHRIS,  Pvalue_CHRIS_1s) %>%
  write.csv("19-Jan-23_Suppl. Table 3_163 replicated SNPs in CHRIS.csv", row.names = F, quote = F)
  

#------------#
# Align alleles and association direction in CHRIA & CKDGen
concordinating_alleles <- function(df) {
  
  df %>%
    as_tibble() %>%
    arrange(CHR, POS38) %>%
    rename(EA_CKDGen_disc  = EA_CKDGen,
           RA_CKDGen_disc  = RA_CKDGen,
           #EAF_CKDGen_disc = EAF_CKDGen,
           MAF_CHRIS_disc  = MAF_CHRIS,
           Locus           = Closest.Gene,
           n_CKDGen        = n_total_sum_CKDGen) %>%
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
           EAF_CKDGen, Beta_CKDGen,    SE_CKDGen, Pvalue_CKDGen, n_CKDGen,
           EAF_CHRIS,  Beta_CHRIS_tmp, SE_CHRIS,  Pvalue_CHRIS, Pvalue_CHRIS_1s)
}

#------------#
# Supplementary table 2: 147 Loci in CKDGen EA and CHRIS

# CHRIS merged with CKDGen EA
repSNPs_EA_Suppl2 <- repSNPs_EA %>%
  concordinating_alleles() %>%
  rename_with(~paste0(.x, "_EA"), ends_with("CKDGen"))

#write.csv(repSNPs_EA_Suppl2, "19-Jan-23_Suppl. Table 2_147 CKDGen Loci.csv", row.names = FALSE)

# merging CKDGen EA-merged CHRIS with CKDGen Trans-Ancestry
repSNPs_MA %>%
  concordinating_alleles() %>%
  rename_with(~paste0(.x, "_MA"), ends_with("CKDGen")) %>%
  inner_join(repSNPs_EA_Suppl2) %>%
  select(Locus, RSID, SNPid, EA_OA, ends_with("CKDGen_MA"), ends_with("CKDGen_EA"), contains("CHRIS")) %>%
  write.csv("07-Dec-23_Suppl._Table2_147_CKDGen_Loci.csv", row.names = FALSE)

#------------#
# Alt+- is a shortcut for <-
#------------#
# Subsetting the leading or most significant SNPs

tagSNPs <- repSNPs %>%
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
  write.csv("10-Jan-23_Table 2_Leading SNPs of the replicated loci.csv", 
            #"10-Jan-23_Supp Table 3_162 Replicated SNPs.csv"
            row.names = FALSE)


#-----------------------------------------------------#
#-------       Testing the correlations       --------
#-----------------------------------------------------#

cor.test(repSNPs$MAF.Diff,
         repSNPs$Effect.CHRIS - repSNPs$Effect.ckdgen,
         use = "complete.obs")

# changing the signs of positive alleles to negative effect
propVars <- c("Effect.ckdgen",
              "Effect.CHRIS",
              "CHRIS.CKDGen.Effect.Ratio")

cor(repSNPs[c("MAF.Diff", "Effect.CHRIS", "Effect.ckdgen")], 
    use = "complete.obs")

#Visualizing the correlation matrix
library(ggcorrplot)

#Select MAF and effect size columns
cor_data <- repSNPs %>% 
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


#-------------#


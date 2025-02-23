#=========================================#
#    Mediation Analysis: Running Steps
#=========================================#

# Last edit on 10/12/2023

#------------#

library(tidyverse)
library(stringr)
library(purrr)
source("08-1_mediation_data.R") # to load vcfReg data file


# outputs
mediation_table <- "11-Jan-23_SNP-wise summary of mediation analysis steps.csv"


#-----------------------------------------------------#
#-------           Analysis Pipeline           -------
#-----------------------------------------------------#

# Mediation analysis steps are as follows:

# step 1: ln(eGFR) ~ SNP + kinship_matrix via emmax model in EPACTS
# step 2: ln(eGFR) ~ SNP + trait    + PC1:10
# step 3:    trait ~ SNP + PC1:10
# step 3:    trait ~ SNP + ln(eGFR) + PC1:10


#-----------------------------------------------------#
#------      Step 1: association with eGFR     -------
#-----------------------------------------------------#

# Step 1 results is indeed the results of GWAS and
# the replication of CKDGen EA in CHRIS 10K done on the servers.
# Commands are available here:
# /home/dghasemisemeskandeh/projects/gwas/replicationAnalysis/scripts/
# The summary or final results of replication analysis can be found in: 
# 


#-----------------------------------------------------#
#------        Step 2: trait as covariate      -------
#-----------------------------------------------------#


# Forming principal components terms
PCs <- paste0("PC", 1:10, collapse = " + ")

#------------#
# preparing the results of step 2 (covariate)
# for merging with step3 (direct)


# Step 2: long format
results_Step2_long <-  getCoefs2table(
  quantVars,
  targets,
  paste("eGFRw.log.Res ~ SNP + trait +", PCs),
  vcfReg
  ) %>%
  group_by(SNPid) %>%
  mutate(outlier = ifelse(is_outlier(Estimate), "Yes", "No")) %>% 
  ungroup()



#-----------------------------------------------------#
#------        Step 3a: trait as outcome       -------
#-----------------------------------------------------#

# Mediation Analysis -> Step 3: trait as the outcome
# The direct association between the SNP and quantitative traits

# Step 3 results: Long format
results_Step3_long <- getCoefs3table(
  vcfReg[quantVars],
  vcfReg[targets],
  paste("trait ~ SNP +", PCs),
  vcfReg) %>%
  mutate(related = ifelse(Pvalue <= 0.05/770, "Yes", "No"))

# Total yes_outlier in Step2 = 1,135; excluding SCr = 973



#-----------------------------------------------------#
#------        Step4: eGFR as mediator         -------
#-----------------------------------------------------#

#Step4: testing eGFR as a mediator for each trait
results_Step4_long <- getCoefs3table(
  vcfReg[quantVars],
  vcfReg[targets],
  paste("trait ~ SNP + eGFRw.log.Res +", PCs),
  vcfReg)



#-----------------------------------------------------#
#------           Mediation Analysis           -------
#-----------------------------------------------------#

# Collecting the 3 Steps of the Mediation Analysis

#step 1: ln(eGFR) ~ SNP + PC1:10
#step 2: ln(eGFR) ~ SNP + trait    + PC1:10
#step 3:    trait ~ SNP + PC1:10
#step 3:    trait ~ SNP + ln(eGFR) + PC1:10


sum3steps_long <- repSNPs %>%
  mutate(EA_OA = paste0(EA_CHRIS_disc, "/", RA_CHRIS_disc)) %>% 
  select(SNPid, Locus, RSID, EA_OA, Beta_CHRIS, SE_CHRIS, Pvalue_CHRIS) %>% #Step1
  inner_join(results_Step2_long, by = c("SNPid", "Locus")) %>%              #Step2
  inner_join(results_Step3_long, by = c("SNPid" = "SNPid",                  #Step3
                                        "Locus" = "Locus",
                                        "Trait" = "Trait"),
             suffix = c("_Step2", "_Step3")) %>%
  inner_join(results_Step4_long, by = c("SNPid" = "SNPid",                  #Step4
                                        "Locus" = "Locus",
                                        "Trait" = "pheno"),
             suffix = c("", "_Step4")) %>%
  mutate(Mediator = ifelse(outlier == "Yes" & related == "Yes", "Yes", "No"),
         Change_P = round((Beta_CHRIS - Estimate_Step2) / Beta_CHRIS * -100, 2)) %>% 
  rename(Estimate_GWAS  = Beta_CHRIS,
         SE_GWAS        = SE_CHRIS,
         Pvalue_GWAS    = Pvalue_CHRIS,
         Estimate_Step4 = Estimate,
         SE_Step4       = SE,
         Pvalue_Step4   = Pvalue) %>%
  select(Locus, SNPid, RSID, EA_OA, ends_with("GWAS"), everything()) %>%
  filter(Mediator == "Yes",
         Trait    != "SCr") %>% 
  as_tibble()


# save mediation results
write.csv(sum3steps_long, file = mediation_table, row.names = FALSE, quote = FALSE)

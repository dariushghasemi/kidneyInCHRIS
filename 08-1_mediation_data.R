#=========================================#
#   Mediation Analysis: Data Preparation
#=========================================#

library(tidyverse)
library(stringr)
library(purrr)


#------------------------------------------------------#
#------------ Merging Genotypes & Phenotypes ----------
#------------------------------------------------------#


# Genotypes file containing the Dosage levels of the alleles of the SNPs
vcf <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\163Targets.txt",
                  col.names = c("AID", "MARKER_ID", "Dosage"), sep = "\t", stringsAsFactors = FALSE)
#------------#

# Principle components
PCs_raw <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\CHRIS.10K.genotypes.evecs", 
                      sep = "\t", stringsAsFactors = FALSE)

PCs_13K <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\CHRIS13K.GT.evecs", sep = "\t", 
                      stringsAsFactors = FALSE)

PCs_13K <- PCs_13K %>% rename(AID = X.IND_ID)

#------------#
# Genotypes file in wide format
vcfmod <- vcf %>%
  mutate_at("MARKER_ID", 
            str_replace_all, 
            ":[A-T-C-G]+:[A-T-C-G]+", "")%>%
  rename(SNPid = MARKER_ID) %>%
  pivot_wider(names_from = "SNPid",
              values_from= "Dosage")%>% 
  select(-c("chr11:78335892", "chr11:78392251"))
#------------#

# Combining dosages + kidney-thyroid traits + PCs
# This dataset is being used for mediation analysis.
# Later we will used this data to filter cases to be 
# used in interaction analyses with TSH and thyroid disease.

vcfReg <- chris[c(
  "AID",
  "Sex",
  "Age",
  "Operator",
  "Municipality",
  "Cancer",
  "Thyroid_cancer",
  "Alteration_thyroid_pregnancy",
  "Thyroid_DrugName",
  "kidneyCancer",
  "Goiter",
  "Operation_thyroid_gland",
  "TSH_cat",
  "eGFRw.log.Res",
  "eGFRw.log",
  "eGFR",
  quantVars)] %>% 
  inner_join(vcfmod,  by = "AID") %>%
  inner_join(PCs_13K, by = "AID", all = FALSE)
#------------#

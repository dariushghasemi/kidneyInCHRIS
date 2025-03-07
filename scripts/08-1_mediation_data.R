#=========================================#
#   Mediation Analysis: Data Preparation
#=========================================#

library(tidyverse)
library(stringr)
library(purrr)
source("~/projects/kidneyInCHRIS/scripts/00-0_configuration.R") #to load quantVars

#------------#
# inputs
path_dosage <- "~/projects/kidneyInCHRIS/inputs/163Targets.txt"
path_pcs    <- "~/projects/kidneyInCHRIS/inputs/CHRIS13K.GT.evecs"
path_clinicals <- "~/projects/kidneyInCHRIS/inputs/chris_q-norm.csv"

#------------------------------------------------------#
#------               Reading data              -------
#------------------------------------------------------#


# Genotypes file containing the Dosage levels of the alleles of the SNPs
vcf <- data.table::fread(path_dosage, col.names = c("AID", "MARKER_ID", "Dosage"))

# Principle components
pcs <- data.table::fread(path_pcs) %>% rename(AID = "#IND_ID")

# load clinical traits
chris <- read.csv(path_clinicals)


#------------#
# Genotypes file in wide format
vcfmod <- vcf %>%
  dplyr::mutate_at("MARKER_ID", str_replace_all, ":[A-T-C-G]+:[A-T-C-G]+", "") %>%
  dplyr::rename(SNPid = MARKER_ID) %>%
  pivot_wider(names_from = "SNPid", values_from= "Dosage") %>% 
  dplyr::select(-c("chr11:78335892", "chr11:78392251")) # remove variants that are not among 163 signif variants


#------------------------------------------------------#
#------      Merging Genotypes & Phenotypes      ------
#------------------------------------------------------#

# Combining dosages + kidney-thyroid traits + PCs
# This dataset is being used for mediation analysis.
# Later we will used this data to filter cases to be 
# used in interaction analyses with TSH and thyroid disease.

vcfReg <- chris %>%
  dplyr::select(
    AID, Sex, Age, Operator, Municipality,
    Cancer, Thyroid_cancer, Alteration_thyroid_pregnancy,
    Thyroid_DrugName, kidneyCancer, Goiter, Operation_thyroid_gland,
    TSH_cat, eGFRw.log.Res, eGFRw.log, eGFR,
    all_of(quantVars)
  ) %>%
  right_join(vcfmod, by = "AID") %>%
  inner_join(pcs, by = "AID") # n=10,758

#------------#
# Storing the variants name (SNPid)
targets <- vcfReg %>% select(starts_with("chr")) %>% colnames()


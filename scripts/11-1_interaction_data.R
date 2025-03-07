
#=========================================#
#   Prepare data to model interaction
#=========================================#

library(tidyverse)
source("~/projects/kidneyInCHRIS/scripts/08-1_mediation_data.R") #to load vcfReg


#-----------------------------------------------------#
#------      Preparing data for interaction    -------
#-----------------------------------------------------#

#Excluding participants with serious 
#Thyroid problems for interaction model

vcfReg_TSHmod <- vcfReg %>%
  dplyr::mutate(Thyroid_DrugName = na_if(Thyroid_DrugName, "")) %>%
  # below filters will remove 416 cases
  dplyr::filter(
    !if_all(c(Thyroid_DrugName, TSH), is.na),                        # Removed are:   4 with NAs
    #Cancer                 != 1 | is.na(Cancer),                    # Removed are: 356 with Cancer
    Thyroid_cancer          != 1 | is.na(Thyroid_cancer),            # Removed are:  16 with Thyroid_cancer (only 3 more removed after removing Cancer cases)
    kidneyCancer            != 1 | is.na(kidneyCancer),              # Removed are:   1 with Kidney cancer
    Goiter                  != 1 | is.na(Goiter),                    # Removed are: 277 with Goitre
    Operation_thyroid_gland != 1 | is.na(Operation_thyroid_gland),   # Removed are: 312 with operation on thyroid gland
    #Alteration_thyroid_pregnancy != 1 | is.na(Alteration_thyroid_pregnancy), # Removed are: 64 removed (but 62 should be eliminated)
    ) %>%
  dplyr::mutate(
    TSH_cat = replace(TSH_cat, Thyroid_DrugName == "Iodine therapy",       "HyperT"), # no cases taking iodine therapy
    TSH_cat = replace(TSH_cat, Thyroid_DrugName == "Levothyroxine sodium", "HypoT"),  # 512 cases taking Levothyn
    TSH_cat = replace(TSH_cat, Thyroid_DrugName == "Propylthiouracil",     "HyperT"), # 1 case  taking Propylthiouracil
    TSH_cat = replace(TSH_cat, Thyroid_DrugName == "Thiamazole",           "HyperT"), # 3 cases taking thiamazole
    #Changing the reference level of TSH_cat to TSH_cat = "2" or "NormT"
    TSH_cat = relevel(as.factor(TSH_cat), ref = 2)
    ) 

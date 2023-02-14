
#=========================================#
#             Data Cleansing
#=========================================#

library(tidyverse)
library(dplyr)
library(data.table)

#------------------------------------------------------#
#------------------- CHRIS data -----------------------
#------------------------------------------------------#

#CHRIS<-read.table("E:...\\MostTraits01\\CHRIS.xlsx", header=T, colClasses="character") # encoding="UTF-8", fileEncoding="utf-16"
#CHRIS<-read.delim2("E:...\\MostTraits01\\traits.txt", colClasses="character")
#CHRISkidney<-read.table("E:...\\MostTraits01+KidneyQue+BMI\\traits.txt", sep='\t')

#importing CHRIS Baseline data on 8th Nov 2021 and merging with Anthropometry features
#On 21th January, 2022, we added drug summary data to CHRIS baseline data for investigating thyroidic patients.

#------------------------------------------------------#

# Importing lab traits CHRIS baseline containing 10,388
CHRISlab  <- fread("D:\\Dariush\\PhD\\Analysis\\Data\\CHRIS_TEMP\\CHRISbaseline_clinical_traits.txt", 
                   sep = "\t", header = T, colClasses = "character", stringsAsFactors = FALSE) # quote="" #%>% writeLines()

# Importing drug questionnaire
CHRISdrug <- fread("D:\\Dariush\\PhD\\Analysis\\Data\\CHRIS_TEMP\\thyroid_drugs_baseline.txt", 
                   sep = "\t", header = T, colClasses = "character", stringsAsFactors = FALSE)

# Importing anthropomorphic traits
CHRISdemo <- fread("D:\\Dariush\\PhD\\Analysis\\Data\\CHRISbaseTotal\\traits.txt", 
                   sep = "\t", header = T, colClasses = "character", stringsAsFactors = FALSE)
#------------------------------------------------------#

#Merging data sets

CHRISbase <-   
  CHRISlab %>% 
  select(c("AID","x0an01","x0an02","x0an03q","x0an04","x0an05","x0bp01","x0bp02","x0bp03","x0bp04a","x0bp05")) %>%
  inner_join(CHRISdrug, by = "AID") %>% #CHRISdrug[, -c("x0dd02a", "x0dd24") #group_by(x0dd06) %>% summarise(n())
  inner_join(CHRISdemo, by = "AID")

#CHRISbase <- merge(data.frame(CHRISlab[,c("AID",    "x0an01", "x0an02", "x0an03q", "x0an04", "x0an05", "x0bp01", "x0bp02", "x0bp03", "x0bp04a", "x0bp05")],
#                               stringsAsFactors = FALSE), CHRISdemo, by = "AID", all.x = TRUE)
#------------------------------------------------------#

# Properly naming the variables
chrisMessy <- data.frame( "AID"                = CHRISbase$AID,
                          "Sex"                = CHRISbase$x0_sex,
                          "Age"                = CHRISbase$x0_age,
                          "Participation"      = CHRISbase$x0_examd,
                          "Operator"           = CHRISbase$x0_opintc,
                          "Municipality"       = CHRISbase$x0_residp,
                          "Height"             = CHRISbase$x0an01,
                          "Weight"             = CHRISbase$x0an02,
                          "BMI"                = CHRISbase$x0an03q,
                          "Body_Fat"           = CHRISbase$x0an04,
                          "Visceral_Fat"       = CHRISbase$x0an05,
                          "SBP"                = CHRISbase$x0bp01, #Systolic Blood Pressure
                          "DBP"                = CHRISbase$x0bp02, #Diastolic BP
                          "Pulse_Rate"         = CHRISbase$x0bp03,
                          "BP_Operator"        = CHRISbase$x0bp04a,
                          "BP_Device"          = CHRISbase$x0bp05,
                          "Pregnant"                 = CHRISbase$x0wo05,
                          "Pregnancy_week"           = CHRISbase$x0wo05a,
                          "Birth_Weight1"            = CHRISbase$x0bi01,
                          "Birth_Weight2"            = CHRISbase$x0bi01a,
                          "Pregnancy1"               = CHRISbase$x0bi02,
                          "Pregnancy2"               = CHRISbase$x0bi02a,
                          "PrePostMature"            = CHRISbase$x0bi03,
                          "BirthType"                = CHRISbase$x0bi04,
                          "Breast_Fed1"              = CHRISbase$x0bi05,
                          "Breast_Fed2"              = CHRISbase$x0bi05a,
                          "CongenitalMalformations1" = CHRISbase$x0bi06,
                          #"CongenitalMalformations2" = CHRISbase$x0bi06a,
                          "HbA1c"             = CHRISbase$x0lp06a,
                          "HbA1c.Ins"         = CHRISbase$x0lp06ac,#0:6700-1:6688
                          "UCr"               = CHRISbase$x0lp01, #UrinaryCreatinine
                          "UCr.Ins"           = CHRISbase$x0lp01c, #0:4420-1:8968
                          "UAlb"              = CHRISbase$x0lp64, #UrinaryAlbumin
                          "UAlb.Ins"          = CHRISbase$x0lp64c, #0:2015-1:11373
                          "SCr"               = CHRISbase$x0lp07, #SerumCreatinine
                          "SCr.Ins"           = CHRISbase$x0lp07c, #0:4420-1:8968
                          "SAlb"              = CHRISbase$x0lp09, #SerumAlbumin
                          "SAlb.Ins"          = CHRISbase$x0lp09c, #0:4420-1:8968
                          "Cancer"            = CHRISbase$x0ca00,
                          "Diabetes"               = CHRISbase$x0dm00,
                          "Diabetes.Type"          = CHRISbase$x0dm02,
                          "Diabetes.Treat"         = CHRISbase$x0dm03,
                          "kidneybydoctor"     = CHRISbase$x0ki00,
                          "glomerulonephritis" = CHRISbase$x0ki01,
                          "glomerTreatment"    = CHRISbase$x0ki01c,
                          "pyelonephritis"     = CHRISbase$x0ki02,
                          "pyelonTreatment"    = CHRISbase$x0ki02c,
                          "renalArteries"      = CHRISbase$x0ki04,
                          "renalArtTreatment"  = CHRISbase$x0ki04c,
                          "congenitalKidney"   = CHRISbase$x0ki05,
                          "kidneyCancer"       = CHRISbase$x0ki06,
                          "kidneyStones"       = CHRISbase$x0ki07,
                          "anotherKidDis"      = CHRISbase$x0ki08,
                          "renalFailure"       = CHRISbase$x0ki09,
                          "kidneyTransplant"   = CHRISbase$x0ki10,
                          "renalSurgery"       = CHRISbase$x0ki16,
                          "kidneyStonesOperat" = CHRISbase$x0ki17,
                          "medProblemOperat"   = CHRISbase$x0ki18,
                          "donatedKidney"      = CHRISbase$x0ki19,
                          "angioplasty"        = CHRISbase$x0ki20,
                          "AnotherRenalSurg"   = CHRISbase$x0ki21,
                          "EverDialysis"       = CHRISbase$x0ki22,
                          "stillDialysis"      = CHRISbase$x0ki23, 
                          "Thyroid_disease"              = CHRISbase$x0th00,
                          "Hyperthyroidism"              = CHRISbase$x0th01,
                          "Hyperthyroidism_dy"           = CHRISbase$x0th01a,
                          "Hypothyroidism"               = CHRISbase$x0th02,
                          "Hypothyroidism_dy"            = CHRISbase$x0th02a,
                          "Goiter"                       = CHRISbase$x0th03,
                          "Goiter_dy"                    = CHRISbase$x0th03a,
                          "Nodule"                       = CHRISbase$x0th04,
                          "Nodule_dy"                    = CHRISbase$x0th04a,
                          "GravesBasedow_disease"        = CHRISbase$x0th05,
                          "GravesBasedow_disease_dy"     = CHRISbase$x0th05a,
                          "Thyroid_cancer"               = CHRISbase$x0th06,
                          "Thyroid_cancer_dy"            = CHRISbase$x0th06a,
                          "Hashimoto_disease"            = CHRISbase$x0th07,
                          "Hashimoto_disease_dy"         = CHRISbase$x0th07a,
                          "Alteration_thyroid_pregnancy" = CHRISbase$x0th08,
                          "Other_thyroid_disease"        = CHRISbase$x0th09,
                          "Other_thyroid_disease_dy"     = CHRISbase$x0th09a,
                          "Other_thyroid_disease_ft"     = CHRISbase$x0th09b,
                          "Radioiodine_therapy"          = CHRISbase$x0th11,
                          "Radioiodine_therapy_year"     = CHRISbase$x0th11a,
                          "Operation_thyroid_gland"      = CHRISbase$x0th12,
                          "Operation_year"               = CHRISbase$x0th12a,
                          "Operation_type_list"          = CHRISbase$x0th12b,
                          "Operation_type_free_text"     = CHRISbase$x0th12c,
                          "Therapy_thyroid_list"         = CHRISbase$x0th13, #1:Yes, radioiodine therapy, 2:Yes, medical or pharmacological therapy, 3: No.
                          "Therapy_thyroid_from_year"    = CHRISbase$x0th13a,
                          "Therapy_thyroid_to_year"      = CHRISbase$x0th13b,
                          "Therapy_thyroid_until_today"  = CHRISbase$x0th13c,
                          "Family_thyroid_diseases"      = CHRISbase$x0th14,
                          "Family_thyroid_Mother"        = CHRISbase$x0th14a,
                          "Family_thyroid_Father"        = CHRISbase$x0th14b,
                          "Family_thyroid_Brothers"      = CHRISbase$x0th14c,
                          "Family_thyroid_Sisters"       = CHRISbase$x0th14d,
                          "Family_thyroid_Sons"          = CHRISbase$x0th14e,
                          "Family_thyroid_Daughters"     = CHRISbase$x0th14f,
                          "Family_thyroid_ft"            = CHRISbase$x0th14g,
                          "Notes_diseases"               = CHRISbase$x0thn1,
                          "Notes_surgery"                = CHRISbase$x0thn2,
                          "Notes_therapy"                = CHRISbase$x0thn3,
                          "Notes_familiarity"            = CHRISbase$x0thn4,
                          "Notes"                        = CHRISbase$x0thnote,
                          "Version"                      = CHRISbase$x0thver, #TBG & AntiTPO
                          "Thyroid_Medication"           = CHRISbase$x0dd24.y, #Yes/No Answers
                          "Thyroid_DrugName"             = CHRISbase$x0dd06,
                          "PTT"                = CHRISbase$x0lp02a, #Prothrombin_Time
                          "INR_PT"             = CHRISbase$x0lp02b, #INR_PT_INR
                          "APTT_ratio"         = CHRISbase$x0lp02c,
                          "APTT"               = CHRISbase$x0lp02d, #Previously labeled "PTT_sec"
                          "Fibrinogen"         = CHRISbase$x0lp03,
                          "AT"                 = CHRISbase$x0lp04, #Antithrombin
                          "BGlucose"           = CHRISbase$x0lp05, #Blood_Glucose
                          "Urate"              = CHRISbase$x0lp08, #Uric_Acid
                          "AST_GOT"            = CHRISbase$x0lp10,
                          "ALT_GPT"            = CHRISbase$x0lp11,
                          "GGT"                = CHRISbase$x0lp12, #Gamma_GT
                          "ALP"                = CHRISbase$x0lp13,
                          "TB"                 = CHRISbase$x0lp14a, #Total_Bilirubin
                          "DB"                 = CHRISbase$x0lp14b, #Direct_Bilirubin
                          "Lipase"             = CHRISbase$x0lp15,
                          "TC"                 = CHRISbase$x0lp16, #Total_Cholesterol
                          "HDL"                = CHRISbase$x0lp17,
                          "LDL"                = CHRISbase$x0lp18,
                          "TG"                 = CHRISbase$x0lp19, #Triglycerides
                          "Sodium"             = CHRISbase$x0lp20,
                          "Potassium"          = CHRISbase$x0lp21,
                          "Chlorine"           = CHRISbase$x0lp22,
                          #"Calcium_mg"         = CHRISbase$x0lp23,
                          #"Calcium_mmol.L"     = CHRISbase$x0lp23a,
                          "Calcium"            = CHRISbase$x0lp24, #Calcium_Corrected_mg
                          "Phosphorus"         = CHRISbase$x0lp25, #Phosphorus_mg
                          #"Phosphorus_mmol"    = CHRISbase$x0lp25a,
                          "Magnesium"          = CHRISbase$x0lp26, #Magnesium_mg
                          #"Magnesium_mmol"     = CHRISbase$x0lp26a,
                          "Iron"               = CHRISbase$x0lp27,
                          "Ferritin"           = CHRISbase$x0lp28,
                          "Transferrin"        = CHRISbase$x0lp29,
                          "TIBC"               = CHRISbase$x0lp30, #T_Iron_Binding
                          "TS"                 = CHRISbase$x0lp31, #Transferrin_saturation
                          "Homocyst"           = CHRISbase$x0lp32, #Homocysteine
                          "CRP"                = CHRISbase$x0lp33, #C_Reactive_Protein
                          "TSH"                = CHRISbase$x0lp35,
                          "FT4"                = CHRISbase$x0lp36, #"T4ng.dL"
                          #'T4pg.mL'            = CHRISbase$x0lp36a,
                          "FT3"                = CHRISbase$x0lp37,
                          "Cortisol"           = CHRISbase$x0lp38,
                          "WBC"                = CHRISbase$x0lp39,
                          "RBC"                = CHRISbase$x0lp40,
                          "HGB"                = CHRISbase$x0lp41,
                          "HCT"                = CHRISbase$x0lp42,
                          "MCV"                = CHRISbase$x0lp43,
                          "MCH"                = CHRISbase$x0lp44,
                          "MCHC"               = CHRISbase$x0lp45,
                          "RDW"                = CHRISbase$x0lp46,
                          "PLT"                = CHRISbase$x0lp47,
                          "MPV"                = CHRISbase$x0lp48,
                          "Neutrophils"        = CHRISbase$x0lp50a, #Neutrophils1
                          "Lymphocytes"        = CHRISbase$x0lp50b, #Lymphocytes1
                          "Monocytes"          = CHRISbase$x0lp50c, #Monocytes1
                          "Eosinophils"        = CHRISbase$x0lp50d, #Eosinophils1
                          "Basophils"          = CHRISbase$x0lp50e, #Basophils1
                          "Neutrophils2"       = CHRISbase$x0lp50f,
                          "Lymphocytes2"       = CHRISbase$x0lp50g,
                          "Monocytes2"         = CHRISbase$x0lp50h,
                          "Eosinophils2"       = CHRISbase$x0lp50i,
                          "Basophils2"         = CHRISbase$x0lp50j,
                          ##"ANA"                = CHRISbase$x0lp51a,
                          #"ANA_pattern_flouroscopico"  = CHRISbase$x0lp51b,
                          #"ANA_titolo1"        = CHRISbase$x0lp51c,
                          #"ANA_pattern"        = CHRISbase$x0lp51d,
                          #"ANA_titolo2"        = CHRISbase$x0lp51e,
                          #"ANA_cromosomi_metafase"     = CHRISbase$x0lp51f,
                          #"ANA_titolo_cromosomi"       = CHRISbase$x0lp51g,
                          #"ANA_observation"            = CHRISbase$x0lp51h,
                          #"ANA_pattern_flouroscopico2" = CHRISbase$x0lp51i,
                          #"ANA_pattern_flouroscopico3" = CHRISbase$x0lp51j,
                          "AntiTPO"             = CHRISbase$x0lp52,
                          "Urine_pH"            = CHRISbase$x0lp53, #pH
                          "UGlucose"            = CHRISbase$x0lp54, #UrinaryGlucose
                          "UProteins"           = CHRISbase$x0lp55, #Proteins
                          "UHGB"                = CHRISbase$x0lp56, #Hemoglobin
                          "Ketone_bodies"       = CHRISbase$x0lp58,
                          ##"Bilirubin"           = CHRISbase$x0lp59,
                          ##"Urobilinogen"        = CHRISbase$x0lp60,
                          "Specific_weight"     = CHRISbase$x0lp61a,
                          ##"Color"               = CHRISbase$x0lp61b,
                          ##"Appearance"          = CHRISbase$x0lp61c,
                          ##"Nitrites"            = CHRISbase$x0lp61d,
                          ##"Leukocytes_esterase" = CHRISbase$x0lp61e,
                          ##"Erythrocytes"        = CHRISbase$x0lp62b,
                          ##"Leukocytes"          = CHRISbase$x0lp62c,
                          ##"Epithelial_cells"    = CHRISbase$x0lp62d,
                          ##"Bacteria"            = CHRISbase$x0lp62e,
                          #"LIPAEMIC_INDEX"      = CHRISbase$x0lp82,
                          #"hemolitic_index"     = CHRISbase$x0lp83,
                          #"ICTERIC_INDEX"       = CHRISbase$x0lp84,
                          #"HIL_ABBOTT"          = CHRISbase$x0lp85,
                          "PTT.Ins"                = CHRISbase$x0lp02ac, #0:3975-1:9413 #Prothrombin_Time.Ins
                          "INR_PT.Ins"             = CHRISbase$x0lp02bc, #0:3975-1:9413 #INR_PT_INR.Ins
                          "APTT_ratio.Ins"         = CHRISbase$x0lp02cc, #0:3975-1:9413
                          "APTT.Ins"               = CHRISbase$x0lp02dc, #0:3975-1:9413 #APTT_sec.Ins
                          "Fibrinogen.Ins"         = CHRISbase$x0lp03c,  #0:3975-1:9413
                          "AT.Ins"                 = CHRISbase$x0lp04c,  #0:3975-1:9413 #Antithrombin.Ins
                          "BGlucose.Ins"           = CHRISbase$x0lp05c, #0:4420-1:8968 #Blood_Glucose.Ins
                          "Urate.Ins"              = CHRISbase$x0lp08c, #0:4420-1:8968 #Uric_Acid.Ins
                          "AST_GOT.Ins"            = CHRISbase$x0lp10c, #0:4420-1:8968
                          "ALT_GPT.Ins"            = CHRISbase$x0lp11c, #0:4420-1:8968
                          "GGT.Ins"                = CHRISbase$x0lp12c, #0:4420-1:8968 #Gamma_GT.Ins
                          "ALP.Ins"                = CHRISbase$x0lp13c, #0:4420-1:8968
                          "TB.Ins"                 = CHRISbase$x0lp14ac,#0:4420-1:8968 #Total_bilirubin.Ins
                          "DB.Ins"                 = CHRISbase$x0lp14bc,#0:4420-1:8968 #Direct_bilirubin.Ins
                          "Lipase.Ins"             = CHRISbase$x0lp15c, #0:4420-1:8968
                          "TC.Ins"                 = CHRISbase$x0lp16c, #0:4420-1:8968 #Total_Cholesterol.Ins
                          "HDL.Ins"                = CHRISbase$x0lp17c, #0:4420-1:8968
                          "LDL.Ins"                = CHRISbase$x0lp18c, #0:4420-1:8968
                          "TG.Ins"                 = CHRISbase$x0lp19c, #0:4420-1:8968 #Triglycerides.Ins
                          "Sodium.Ins"             = CHRISbase$x0lp20c, #0:4420-1:8968
                          "Potassium.Ins"          = CHRISbase$x0lp21c, #0:4420-1:8968
                          "Chlorine.Ins"           = CHRISbase$x0lp22c, #0:4420-1:8968
                          #"Calcium_mg.Ins"         = CHRISbase$x0lp23c, #0:4420-1:8968
                          "Calcium.Ins"            = CHRISbase$x0lp24c,#0:4420-1:8968 #Calcium_Corrected_mg.Ins
                          "Phosphorus.Ins"         = CHRISbase$x0lp25c, #0:4420-1:8968
                          "Magnesium.Ins"          = CHRISbase$x0lp26c, #0:4420-1:8968
                          "Iron.Ins"               = CHRISbase$x0lp27c, #0:4420-1:8968
                          "Ferritin.Ins"           = CHRISbase$x0lp28c, #0:4420-1:8968
                          "Transferrin.Ins"        = CHRISbase$x0lp29c, #0:4420-1:8968
                          "TIBC.Ins"               = CHRISbase$x0lp30c, #0:4420-1:8968 #T_Iron_Binding.Ins
                          "TS.Ins"                 = CHRISbase$x0lp31c,#0:4420-1:8968 #Transferrin_saturation.Ins
                          "Homocyst.Ins"           = CHRISbase$x0lp32c, #0:10750-1:555-"":2083 #Homocysteine.Ins
                          "CRP.Ins"                = CHRISbase$x0lp33c, #0:4420-1:8968 #C_Reactive_Protein.Ins
                          "TSH.Ins"                = CHRISbase$x0lp35c, #0:4420-1:8968
                          "FT4.Ins"                = CHRISbase$x0lp36c, #0:4420-1:8968 #FT4.Ins
                          "FT3.Ins"                = CHRISbase$x0lp37c, #0:4420-1:8968
                          "Cortisol.Ins"           = CHRISbase$x0lp38c, #0:4420-1:8968
                          "WBC.Ins"                = CHRISbase$x0lp39c, ##0:4150-1:9238
                          "RBC.Ins"                = CHRISbase$x0lp40c, ##0:4150-1:9238
                          "HGB.Ins"                = CHRISbase$x0lp41c, ##0:4150-1:9238
                          "HCT.Ins"                = CHRISbase$x0lp42c, ##0:4150-1:9238
                          "MCV.Ins"                = CHRISbase$x0lp43c, ##0:4150-1:9238
                          "MCH.Ins"                = CHRISbase$x0lp44c, ##0:4150-1:9238
                          "MCHC.Ins"               = CHRISbase$x0lp45c, ##0:4150-1:9238
                          "RDW.Ins"                = CHRISbase$x0lp46c, ##0:4150-1:9238
                          "PLT.Ins"                = CHRISbase$x0lp47c, ##0:4150-1:9238
                          "MPV.Ins"                = CHRISbase$x0lp48c, ##0:4150-1:9238
                          "Neutrophils.Ins"        = CHRISbase$x0lp50ac, #0:4150-1:9238 Neutrophils1.Ins
                          "Lymphocytes.Ins"        = CHRISbase$x0lp50bc, #0:4150-1:9238
                          "Monocytes.Ins"          = CHRISbase$x0lp50cc, #0:4150-1:9238
                          "Eosinophils.Ins"        = CHRISbase$x0lp50dc, #0:4150-1:9238
                          "Basophils.Ins"          = CHRISbase$x0lp50ec, #0:4150-1:9238 Basophils1.Ins
                          "Neutrophils2.Ins"       = CHRISbase$x0lp50fc, #0:4150-1:9238
                          "Lymphocytes2.Ins"       = CHRISbase$x0lp50gc, #0:4150-1:9238
                          "Monocytes2.Ins"         = CHRISbase$x0lp50hc, #0:4150-1:9238
                          "Eosinophils2.Ins"       = CHRISbase$x0lp50ic, #0:4150-1:9238
                          "Basophils2.Ins"         = CHRISbase$x0lp50jc, #0:4150-1:9238
                          ##"ANA.Ins"                = CHRISbase$x0lp51ac,   #0:4600-1:8788
                          "AntiTPO.Ins"            = CHRISbase$x0lp52c,    #0:6786-1:6602
                          "Urine_pH.Ins"           = CHRISbase$x0lp53c, ##0:7358-6030 #pH.Ins
                          "UGlucose.Ins"           = CHRISbase$x0lp54c, ##0:7358-6030 #UrinaryGlucose.Ins
                          "UProteins.Ins"          = CHRISbase$x0lp55c, ##0:7358-6030 #Proteins.Ins
                          "UHGB.Ins"               = CHRISbase$x0lp56c, ##0:7358-6030 #Hemoglobin.Ins
                          "Ketone_bodies.Ins"      = CHRISbase$x0lp58c, ##0:7358-6030
                          "Bilirubin.Ins"          = CHRISbase$x0lp59c, ##0:7358-6030
                          "Urobilinogen.Ins"       = CHRISbase$x0lp60c, ##0:7358-6030
                          "Specific_weight.Ins"    = CHRISbase$x0lp61ac, #0:7358-6030
                          "Color.Ins"              = CHRISbase$x0lp61bc, #0:7358-6030
                          "Appearance.Ins"         = CHRISbase$x0lp61cc, #0:7358-6030
                          "Nitrites.Ins"           = CHRISbase$x0lp61dc, #0:7358-6030
                          "Leukocytes_esterase.Ins"= CHRISbase$x0lp61ec, #0:7358-6030
                          #"Erythrocytes.Ins"       = CHRISbase$x0lp62bc, #0:7358-6030
                          #"Leukocytes.Ins"         = CHRISbase$x0lp62cc, #0:7358-6030
                          #"Epithelial_cells.Ins"   = CHRISbase$x0lp62dc, #0:7358-6030
                          "Bacteria.Ins"           = CHRISbase$x0lp62ec, #0:7358-6030
                          stringsAsFactors = FALSE
)

#--------------------------------------------------------------------------#
#------------------------------ data wrangling ----------------------------#
#--------------------------------------------------------------------------#

#-------- Descriptive Stat -------- 

str(chrisMessy)
summary(chrisMessy)
lapply(chris, sd, na.rm = TRUE)
#---------#

lapply(
  list("UAlb", "UCr", "SAlb", "SCr", "AntiTPO"), function (i) {
    sort(unique(chrisMessy[,i]))[1:10]
  }
)#levels(chris[,i])[1:10]
#---------#

#Table 1 for paper
varsTable1 <- list("eGFR", "Age", "TSH", "FT3", "FT4", "Magnesium", "APTT")

#---------#
#Counting Instruments for Suppl. Table 1 of the paper
sd(chris[,paste0(varsTable1[4],".Ins")], na.rm = T)
#---------#

lapply(varsTable1, #quantVars
       function (i) {
         Mean     = mean(vcfReg[,i],           na.rm = T) %>% round(2)
         SD       = sd(vcfReg[,i],             na.rm = T) %>% round(2)
         Median   = median(vcfReg[,i],         na.rm = T) %>% round(2)
         Q1       = quantile(vcfReg[,i], 0.25, na.rm = T) %>% round(2)
         Q3       = quantile(vcfReg[,i], 0.75, na.rm = T) %>% round(2)
         IQR      = IQR(vcfReg[,i],            na.rm = T) %>% round(2)
         Trait    = unique(i)
         Missed   = sum(is.na(vcfReg[,i]))
         Missed.p = (sum(is.na(vcfReg[,i])) / length(vcfReg[,i])) %>% round(2)
         data.frame(Trait, Mean, SD, Median, Q1, Q3, Missed, Missed.p) %>% 
           unlist()
       }
) %>%
  bind_rows(.) %>% View()
write.csv("29-Dec-2022_Supplementary Table 1 of paper.csv", quote = F)
#------------------------------------------------------#

# Finding Loss of Detection (LOC) values
lapply(chris[,108:171], function (i) {sort(unique(i))[1:10]})
lapply(seq(4,12,2),     function (i) sort(unique(chris[,i])))
lapply(chris[,1:10],    function (i) {sum(is.na(i))})
#---------#

sort(unique(chrisMessy$HbA1c))
sort(unique(chrisMessy$UrinaryCreatinine))[1:10]
sort(unique(chrisMessy$UrinaryAlbumin))[1:10]
sort(unique(chrisMessy$SerumCreatinine))
sort(unique(chrisMessy$SerumAlbumin))
#---------#

sum(is.na(chris$HbA1c))
sum(chris == -9989)
sum(chris == -9998)
sum(chris == -99)
sum(chris == -98)
sum(chris == -89)
sum(chris == -88)
sum(chris$UrinaryAlbumin == -89)
#---------#

# Loss of Detection values
as.data.frame(table(chris$UrinaryAlbumin))[1:8,]
length(which(chrisMessy$UrinaryAlbumin == "< 2.0"))
length(which(chrisMessy$UrinaryAlbumin == "< 2.26"))
sum(is.na(chris$UrinaryAlbumin))
#---------#

#contingacy table
chrisMessy %>% 
  group_by(Sex)%>% 
  count(TSH.Ins) %>% 
  mutate(prop = prop.table(n)) %>% 
  spread(TSH.Ins, n, fill = 0) 
#---------#

chris %>%
  group_by(UrinaryAlbumin.Ins) %>%
  summarise(UrinAlb.M  = mean(UrinaryAlbumin, na.rm= TRUE), 
            UrinAlb.sd =   sd(UrinaryAlbumin, na.rm= TRUE))
#---------#

chris %>% group_by(TSH.Ins) %>% summarise(TSH, mean)
#---------#

library(janitor)
mtcars %>% tabyl(cyl, gear)
tabyl(c("hi", "med", "med", "lo"))
#---------#

# The old way for changing column names
#names(chris)[which(names(chris)=="UrinAlbumin")] <- "UrinaryAlbumin"


#------------------------------------------------------#
#------------------- Missing Values -------------------
#------------------------------------------------------#

library(naniar)
library(VIM)

#dt$Age[dt$Age == 99] <- NA
#chris %>% replace_with_na(replace = list(
#  HbA1c = c(-89),
#  UrinCreatinine= c(-89, -99),
#  UrinAlbumin= c(-99), 
#  SerumCreatinine = c(-89), 
#  Diabetes.Type = c(-99),
#  Diabetes.Treat = c(-99)
#  ))


#chris<- chris %>% replace_with_na_all(condition = ~.x == c(-99))
#chris<- chris %>% replace_with_na_all(condition = ~.x == c(-89))
#chris<- chris %>% replace_with_na_all(condition = ~.x == c(-88))


# Plotting MVs with naniar
gg_miss_var(chris[, -1],show_pct = TRUE)

# Plotting MVs with VIM
res <- summary(aggr(chris[,c(2,3,4,6,8,10,12,14,15,16,17)],
                    sortVar=TRUE, prop=FALSE, numbers=TRUE))$combinations
matrixplot(chris[, c(-1)])
matrixplot(chris[, 2:10])



#------------------------------------------------------#
#------------------ Missing imputation ----------------
#------------------------------------------------------#

# Nt used in this analysis

#library(qdap)
#lookup(chris$UrinCreatinine, c(-88, -89, -99), c(NA, NA, NA), missing = NULL)

#chris$UrinAlbumin[chris$UrinaryAlbumin == "< 2.0"] <- "2"
#chris$UrinaryAlbumin == "< 2.26"] = "2"

#UAlblev<- list('-89'="NA", '< 2.0'="2", '< 2.26'="2")
#chris$UrinaryAlbumin <- recode(chris$UrinaryAlbumin, !!!UAlblev)

#library(plyr) #revalue()
#mapvalues(chris$UrinaryAlbumin, c("-89", "< 2.0", "< 2.26"), c(NA, 2, 2))

#library(dplyr)
#chris$UrinAlbumin<-replace(chris$UrinAlbumin, "-89", NA)
#chris <- chris %>% mutate(UrinAlbumin = replace(UrinAlbumin, UrinAlbumin == -89, NA))



#------------------------------------------------------#
#-------------------- Data cleansing ------------------
#------------------------------------------------------#

library(dplyr)
library(naniar)
library(stringr)
#---------#

# Recoding Sex
# Replacing LOD signs (e.g <1.5) with LOD (e.g 1.5)
# Replacing missing signs (e.g. -88) with NA
# Filling empty info for lab traits with NA

chris <- 
  chrisMessy %>%
  mutate(Sex           = replace(Sex, Sex == "2", "0"))%>%
  mutate(BP_Operator   = replace(BP_Operator, BP_Operator == "-89", 0))%>%
  #mutate(Homocyst.Ins  = replace(Homocyst.Ins, Homocyst.Ins == "", "0"))%>%
  mutate(APTT.Ins      = case_when(APTT.Ins == "0" ~ "1",
                                   APTT.Ins == ""  ~ "0"))%>%
  #mutate(UrinaryAlbumin = replace(UrinaryAlbumin, UrinaryAlbumin == "< 2.0",  2))%>%
  #mutate_all(HbA1c         = replace(HbA1c, HbA1c == "-89", NA))%>%
  mutate(across(where(is.character), ~na_if(., "-88")))%>%
  mutate(across(where(is.character), ~na_if(., "-89")))%>%
  mutate(across(where(is.character), ~na_if(., "-98")))%>%
  mutate(across(where(is.character), ~na_if(., "-99")))%>%
  mutate(across(where(is.character), ~na_if(., "-9998")))%>%
  mutate(across(where(is.character), ~na_if(., "-9989")))%>%
  mutate_if(is.character, str_replace_all, pattern = "< ", "")%>%
  mutate_if(is.character, str_replace_all, pattern = "> ", "")%>%
  mutate_if(is.character, str_replace_all, pattern = ">", "")%>%
  mutate_if(is.character, str_replace_all, pattern = "<", "")%>%
  mutate(UProteins = replace(UProteins, UProteins == "negative", 0))%>%
  mutate(Ketone_bodies = replace(Ketone_bodies, Ketone_bodies == "negative", 0))%>%
  mutate(Operator      = replace(Operator, is.na(Operator), 0))
#---------#

# Changing the format of non-numeric variables
#chris$Participation         <- as.Date(chris$Participation, tryFormats = c("%d/%m/%Y"))
chris$Diabetes              <- as.factor(chris$Diabetes)
chris$Diabetes.Type         <- as.factor(chris$Diabetes.Type)
chris$Operator              <- as.factor(chris$Operator)
chris$Municipality          <- as.factor(chris$Municipality)
#---------#

# Defining quantitative traits for phenome-wide scan
quantVars <- c(#"Age",
  "Height","Weight","BMI","Body_Fat","Visceral_Fat",
  "SBP", "DBP", "Pulse_Rate", "HbA1c", #"Birth_Weight1", "Pregnancy1",
  "SAlb", "SCr", "UAlb", "UCr", "UACR",
  "PTT", "INR_PT", "APTT_ratio", "APTT", "Fibrinogen", "AT", "BGlucose",
  "Urate", "AST_GOT", "ALT_GPT", "GGT", "ALP", "TB", "DB", "Lipase", "TC",
  "HDL", "LDL", "TG", "Sodium", "Potassium", "Chlorine", #"Calcium_mg",
  "Calcium", "Phosphorus", "Magnesium", "Iron", "Ferritin", 
  "Transferrin", "TIBC" , "TS", "Homocyst", "CRP", "TSH", "FT3", "FT4", 
  "Cortisol", "WBC","RBC","HGB","HCT", "MCV","MCH","MCHC","RDW","PLT","MPV",
  "Neutrophils","Lymphocytes","Monocytes","Eosinophils","Basophils",
  #"Neutrophils2","Lymphocytes2","Monocytes2","Eosinophils2","Basophils2",
  "AntiTPO","Urine_pH","UGlucose","UProteins","UHGB"
)
#---------#

#changing the class of quantitative variables into numeric
chris[quantVars] <- lapply(chris[quantVars], as.numeric)
#---------#

#checking replacement
merge(CHRISbase[1:1100, c("AID", "x0lp64", "x0lp52")], 
      chris[1:1100, c("AID", "UAlb",   "AntiTPO")], by = "AID")[1000:1030,]

#Eliminating unnecessary data-sets
remove(CHRISbase, CHRISdemo, CHRISlab)
#---------#

#intersect(colnames(chris), quantVars)
#setdiff(colnames(chris), quantVars)

#"|> [0-9]+|>[0-9]+|<[0-9]+" | "< \\d+.\\d+"
#chris$UrinaryAlbumin <- as.character(chris$UrinaryAlbumin)
#na_strings <- c("-89", "-89", "-98", "-99")
#------------------------------------------------------#

#changing the class of variables into integer

#changeClass(chris, c("chr", rep("int",1:10)))
#mydf[,2:3] <- lapply(mydf[,2:3], as.factor)
#chris[quantVars] <- gsub(",", "", chris[quantVars])
#sapply(mydata, as.numeric)
#---------#

# Installing a package from Github
#install.packages("remotes")
#remotes::install_github("hstojic/hfunk")
#devtools::install_github("hstojic/hfunk")
#chris %>% mutate_at(vars(c(52:112, 113:)), as.numeric)




#------------------------------------------------------#
#----------------- Making new Variables ---------------
#------------------------------------------------------#

library(tidyverse)
library(lubridate)

#Ethnicity for eGFR formula: 0 for Europeans, and 1 for Africans
chris$ethnicity <- 0

#Checking the format of Participation
cbind(chris$Participation,
      as.Date(chris$Participation, tryFormats = c("%d/%m/%Y")),
      CHRISbase$x0_examd) %>% View()

#Participation period for entering CHRIS study in weeks
chris$week <- interval(
  as.Date("15/08/2011", tryFormats = c("%d/%m/%Y")),
  as.Date(chris$Participation, tryFormats = c("%d/%m/%Y"))
) %/% weeks(1)

chris$week <- as.factor(chris$week)
#------------------------------------------------------#

#checking week calculation
head(data.frame("date" = chris$Participation,
                "week" = difftime(chris$Participation, 
                                  as.Date("15/08/2011", tryFormats = c("%d/%m/%Y")),
                                  units = c("weeks")),
                "weekr"= chris$week))

chris[chris$Participation=="2011-08-24" | 
        chris$Participation=="2011-08-29" | 
        chris$Participation=="2011-08-30" |
        chris$Participation=="2011-09-01" |
        chris$Participation=="2011-09-06" |
        chris$Participation=="2018-12-21", 
      c("Participation", "week")]
#------------------------------------------------------#


# Recoding Sex for eGFR formula
#chris$Sex <- dplyr::recode(chris$Sex, `1` = 1L, `2` = 0L) #1:Males 0:Females

# Converting date participation to week
# Sequential week order since the beginning of the study
#chris$week <- difftime(as.Date(chris$Participation, tryFormats = c("%d/%m/%Y")),
#                       as.Date("01/08/2011",        tryFormats = c("%d/%m/%Y")),
#                       units = c("weeks")) %>% floor_date(as.numeric())#floor_date()


#round_date(chris$week, unit = "week", week_start = getOption("lubridate.week.start", 7))
#chris$week <- as.numeric(strftime(as.POSIXlt(chris$Participation, tryFormats=c("%d/%m/%Y")), format="%W"))
#head(data.frame(CHRISbase$x0_examd, "col1"=as.numeric(strftime(CHRISbase$x0_examd, format = "%V")), "col2"=strftime(as.POSIXlt(CHRISbase$x0_examd), format="%W"), "week"=chris$week))

# or using >>> library(lubridate)
#week(ymd(CHRIS$x0_examd)) 
#data.frame(CHRIS$x0_examd, as.numeric(strftime(CHRIS$x0_examd, format = "%V")))

#chris$UrinaryAlbumin.Sensor <- ifelse(chris$UrinaryAlbumin == '< 2.0' | chris$UrinAlbumin == '< 2.26', 1, 0)

#library(kimisc)
#chris$UACR.Ins <- with(chris, coalesce.na(UrinAlbumin.Ins, UrinCreatinine.Ins))
#chris$UACR.Ins <- ifelse(chris$UrinAlbumin.Ins==0 & chris$UrinCreatinine.Ins==0, 0, 1)
#------------------------------------------------------#



#------------------------------------------------------#
#------------------------------------------------------#
#------------------------------------------------------#






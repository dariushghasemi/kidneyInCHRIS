#=========================================#
#    Mediation Analysis: Running Steps
#=========================================#

# Last edit on 10/12/2023

#------------#

library(tidyverse)
library(stringr)
library(purrr)

#-----------------------------------------------------#
#-------           Analysis Pipeline           -------
#-----------------------------------------------------#

# Mediation analysis steps are as follows:

# step 1: ln(eGFR) ~ SNP + kinship_matrix via emmax model in EPACTS
# step 2: ln(eGFR) ~ SNP + trait    + PC1:10
# step 3:    trait ~ SNP + PC1:10
# step 3:    trait ~ SNP + ln(eGFR) + PC1:10

#------------#
# Storing the variants name (SNPid)
targets <- vcfReg %>% select(starts_with("chr")) %>% colnames()

# most CHRIS significant variant at 11 replicated loci
leadingSNPs <- c("chr1:10670853",  "chr2:15648568", 
                 "chr4:76492991",  "chr5:39385539",
                 "chr5:177386403", "chr7:77733187", 
                 "chr8:23894869",  "chr9:68819791",
                 "chr11:78339803", "chr15:98729803")


# CKDGen lead variant at 11 replicated loci
CKDGenSNPs <- c("chr1:10670853",  "chr2:15642347",
                "chr4:76480299",  "chr5:39377763",
                "chr5:177386403", "chr7:77793266",
                "chr8:23890907",  "chr9:68817258",
                "chr11:78312060", "chr15:98733292")



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


# Step 2 function: trait-adjusted effect of variants on eGFRcrea
getCoefs2table <- function(mytrait, mytarget, myformula, data)
  {
  map2_df(
    .x = mytrait,
    .y = myformula,
    function(mytrait, myformula)
      {
      trait <- data[, mytrait]
      res0  <- map_df(
        mytarget,
        function(mySNP)
          {
          SNP <- data[, mySNP]
          myformula <- as.formula(myformula)
          m <- lm(myformula, data = data)
          p <- summary(m)$coefficients[2, c(1,2,4)] #Beta, se, Pvalue
          names(p) <- str_replace_all(
            names(p),
            c("Estimate"   = "Estimate",
              "Std. Error" = "SE",
              "Pr\\([^\\(]*\\)" = "Pvalue")
            )
          p$SNPid <- mySNP
          return(p)
        }
      )
      res0$Trait <- mytrait
      res1 <- repSNPs %>%
        select(SNPid, Locus) %>%
        inner_join(res0, by = c("SNPid")) %>%
        select(SNPid, Locus, Trait, everything())
      return(res1)
      }
  )
}

#------------#
#outlier detection function for SNPs effects of eGFRw.res.ln
is_outlier <- function(x) { 
  return(x < quantile(x, 0.10) - 1.5 * IQR(x) | x > quantile(x, 0.90) + 1.5 * IQR(x))
}

#------------#
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

#------------#

# Step 2 function: check SNP association with a trait
getCoefs3table <- function(mytrait, mytarget, myformula, mydata){
    res3 <- map2_dfr(
      .x = mytrait,
      .y = myformula,
      .f = function(trait, formula){
        res1 <- map_dfr(
          .x = mytarget,
          .f = function(SNP)
          {
            m <- lm(as.formula(formula), data = mydata)
            p <- summary(m)$coefficients[2, c(1,2,4)]
            return(p)
          }
        )
        res2 <- as.data.frame(res1)
        colnames(res2) <- c("Estimate",
                            "SE",
                            "Pvalue")
        Nsnp       <- length(mytarget)
        res2$SNPid <- rep(colnames(mytarget)[1:Nsnp], 1)
        return(res2)
      }
    )
    res3$pheno <- rep(colnames(mytrait), each = length(unique(res3$SNPid)))
    res4 <- repSNPs %>% 
      select(SNPid, Locus) %>% 
      inner_join(res3, by = "SNPid")
    return(res4)
}

#------------#

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


write.csv(sum3steps_long,
          "11-Jan-23_SNP-wise summary of mediation analysis steps.csv",
          row.names = FALSE, quote = FALSE)

#-----------------------------------------------------#
#                       The end                       #
#-----------------------------------------------------#

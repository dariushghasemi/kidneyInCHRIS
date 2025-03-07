

library(tidyverse)

#------------#
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

#------------#
# lead SNPs in Table 2 of paper
# CHRIS: the most significant variant at 11 replicated loci -> equivalent to tagSNPs
target_snps <- c(
  "chr1:10670853",  "chr2:15642347", # "chr2:15648568", 
  "chr4:76480299",  "chr5:39393631", # "chr4:76492991", "chr5:39385539",
  "chr5:177386403", "chr7:77714744", # "chr7:77733187", 
  "chr8:23885208",  "chr9:68817258", # "chr8:23894869", "chr9:68819791",
  "chr11:78324786", "chr15:98733292",# "chr11:78339803", "chr15:98729803"
  "chr16:20381010"
)



# CKDGen: lead variant at 11 replicated loci
CKDGenSNPs <- c(
  "chr1:10670853",  "chr2:15642347",
  "chr4:76480299",  "chr5:39377763",
  "chr5:177386403", "chr7:77793266",
  "chr8:23890907",  "chr9:68817258",
  "chr11:78312060", "chr15:98733292"
  )


# Significant loci in CHRIS
loci_CHRIS <- c(
  "CASZ1", "DDX1", "SHROOM3", "DAB2", "SLC34A1", 
  "TMEM60","STC1", "PIP5K1B", "GAB2", "IGF1R", "PDILT"
  )


#------------#
#Changing the order and the label of the Locus top SNPs with rsID
locus_factor_rsid <- function(x) {
  factor(
    x,
    levels = c(
      "CASZ1", "DDX1", "SHROOM3", "DAB2", "SLC34A1",
      "TMEM60", "STC1", "PIP5K1B", "GAB2", "IGF1R", "PDILT"
      ),
    labels = c(
      "CASZ1\nrs74748843",
      "DDX1\nrs807624",
      "SHROOM3\nrs28817415",
      "DAB2\nrs10062079",
      "SLC34A1\nrs3812036",
      "TMEM60\nrs57514204",
      "STC1\nrs819196",
      "PIP5K1B\nrs2039424",
      "GAB2\nrs7113042",
      "IGF1R\nrs59646751",
      "UMOD-PDILT\nrs77924615"
      )
    )
}

#Changing the order and the label of the Locus top SNPs
locus_factor_snpid <- function(x) {
  factor(
    x,
    levels = c(
      "CASZ1", "DDX1", "SHROOM3", "DAB2", "SLC34A1",
      "TMEM60", "STC1", "PIP5K1B", "GAB2", "IGF1R", "PDILT"
      ),
    labels = c(
      "CASZ1\n1:10670853",
      "DDX1\n2:15642347",
      "SHROOM3\n4:76480299",
      "DAB2\n5:39393631",
      "SLC34A1\n5:177386403",
      "TMEM60\n7:77714744",
      "STC1\n8:23885208",        
      "PIP5K1B\n9:68817258",
      "GAB2\n11:78324786",
      "IGF1R\n15:98733292",
      "PDILT\n16:20381010"
      )
    )
}

#----------#
# recode variants
gene_coder <- function(x){
  factor(
    x,
    levels = c(
      "rs28394165", "rs10025351", "rs28817415", "rs4859682", "rs13146355", "rs3812036"
    ),
    labels = c(
      "rs28394165\nSHROOM3", "rs10025351\nSHROOM3", "rs28817415\nSHROOM3", "rs4859682\nSHROOM3", "rs13146355\nSHROOM3", "rs3812036\nSLC34A1"
    )
  )
}


#------------#
#outlier detection function for SNPs effects of eGFRw.res.ln
is_outlier <- function(x) { 
  return(x < quantile(x, 0.10) - 1.5 * IQR(x) | x > quantile(x, 0.90) + 1.5 * IQR(x))
}

#------------#
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
# Step 3 function: check SNP association with a trait
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
# get variants effect in lm(outcome ~ variant + ...)
getCoefs3lm <- function(mytrait, mytarget, myformula, mydata){
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
      colnames(res2) <- c("Estimate", "SE", "Pvalue")
      Nsnp       <- length(mytarget)
      res2$SNPid <- rep(colnames(mytarget)[1:Nsnp], 1)
      return(res2)
    }
  )
  return(res3)
}

#------------#

#------------#

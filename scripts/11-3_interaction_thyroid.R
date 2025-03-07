

out_tbl_quantitativ <- "13-Feb-23_Interaction of kidney with TSH (eGFRw.log-TSH) adjusted for sex age PCs.csv"
out_tbl_categorized <- "13-Feb-23_Interaction of kidney with thyroid disease (eGFRw.log-TSH_cat) adjusted for sex age PCs.csv"
out_tbl_centralized <- "01-Feb-23_Kideny interaction with thyroid adj for sex and age (SNP-TSH)_wide format.csv"

#------------------#
# variables and constants

# Forming principal components term of the model 
PCs <- paste0("PC", 1:10, collapse = " + ")


#------------------#
# function to run interaction models
regresModel <- function(data, mySNP, formula, vecTraits){
  
  SNP <- data[, mySNP]
  #SNPid <- as.character(noquote(mySNP))
  m <- lm(as.formula(formula), data = data)
  p <- summary(m)$coefficients[vecTraits, c(1,2,4)] #c(2:3,16) applies to TSH; c(2:4,17:18)applies to TSH_cat
  r <- as.data.frame(p) %>%
    rename(Beta = Estimate, SE = `Std. Error`, Pvalue = `Pr(>|t|)`) %>%
    rownames_to_column(var = "Term") %>%
    mutate(associated = ifelse(Pvalue < 0.05/11, "Yes", "No")) %>%
    pivot_wider(#id_cols     = SNPid,
      names_from  = Term,
      values_from = c(Beta, SE , Pvalue, associated),
      names_glue  = "{Term}_{.value}") %>%
    #as.data.frame() %>%
    mutate("SNPid" = mySNP) %>%
    inner_join(repSNPs[c("SNPid", "Locus")],
               by = "SNPid") %>% 
    select(SNPid, Locus, starts_with("SNP_"), starts_with("TSH_"), everything())
  return(r)
}

#-----------------------------------------------------#
#------            TSH-SNP interaction          ------
#-----------------------------------------------------#

# Test interaction model
# summary(
#   lm(
#     eGFRw.log.Res ~ `chr1:10599281` * TSH + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
#     data = vcfReg_TSHmod)
# )$coefficients[c(2:3,16), c(1,2,4)]

#---------#
# Testing regression model function
#regresModel(vcfReg_TSHmod,
#            targets[1],  #vcfReg_TSHmod[targets[1:3]]
#            paste("eGFRw.log.Res ~ SNP * TSH + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))


#---------#
# Linear interaction: TSH-SNP

# Iterating intearction model
results_Kidney_TSH <- map_dfr(
  targets, function(SNP) {
    regresModel(
      vcfReg_TSHmod,
      SNP, #vcfReg_TSHmod[,SNP],
      paste("eGFRw.log ~ SNP * TSH + Sex + Age +", PCs),
      c(2:3,16)
      )
  }
  )

#---------#
# Save the results of linear interaction model
write.csv(results_Kidney_TSH, file = out_tbl_quantitativ, row.names = FALSE, quote = FALSE)

# #---------#
# # Table XX of the paper 
# results_Kidney_TSH %>%
#   filter(SNPid %in% tagSNPs$SNPid) %>% View



#-----------------------------------------------------#
#-------------- SNP:TSH_cat interaction ---------------
#-----------------------------------------------------#


# Non-linear interaction: TSH_cat-SNP

results_Kidney_TSH_cat <- map_dfr(
  targets, 
  function(SNP) {
    regresModel(
      vcfReg_TSHmod,
      SNP,
      paste("eGFRw.log ~ SNP * TSH_cat + Sex + Age+", PCs),
      c(2:4,17:18)
    )
    }
  ) %>%
  dplyr::rename_with(~ gsub("TSH_catHyperT", "Hyperthyrodism", .x, fixed = TRUE)) %>%
  dplyr::rename_with(~ gsub("TSH_catHypoT",  "Hypothyrodism",  .x, fixed = TRUE)) %>%
  dplyr::select(
    SNPid, 
    Locus, 
    starts_with(
      c(
        "SNP_",
        "Hyperthyrodism_",
        "Hypothyrodism_",
        "SNP:Hyperthyrodism_",
        "SNP:Hypothyrodism_"
        )
      )
    )

#---------#

# Save the results of non-linear interaction model
write.csv(results_Kidney_TSH_cat, file = out_tbl_categorized, row.names = FALSE, quote = FALSE)
#---------#

# Table XX of the paper 
results_Kidney_TSH_cat %>% 
  filter(SNPid %in% tagSNPs$SNPid) %>% View



#-----------------------------------------------------#
#--------------- TSH_cen-SNP interaction -------------
#-----------------------------------------------------#

# Model the interaction in different levels of TSH for STC1 locus

# By default, scale() function will standardize the 
# data (mean zero, unit variance). To indicate that
# we just want to subtract the mean, we need to turn 
# off the argument scale = FALSE.

#results_Kidney_TSH_cen <- 
vcfReg_TSHmod %>%
  dplyr::select(AID, eGFR, eGFRw.log, TSH, TSH_cat, Age, Sex, starts_with(c("chr8", "PC"))) %>% 
  # mutate(eGFR.log = log(eGFR),
  #        #TSH_cen  = scale(TSH, scale = F)
  #        ) %>% 
  pivot_longer(cols      = starts_with("chr"),
               names_to  = c("SNPid"),
               values_to = c("Dosage")) %>%
  group_by(SNPid) %>%
  nest() %>%
  dplyr::mutate(
    model_tsh = data  %>% map(~lm(eGFRw.log ~ Dosage * TSH     + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)),
    model_cat = data  %>% map(~lm(eGFRw.log ~ Dosage * TSH_cat + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)),
    tidy_tsh  = model_tsh %>% map(broom::tidy),
    tidy_cat  = model_cat %>% map(broom::tidy)) %>%
  unnest(tidy_tsh) %>%
  #unnest(tidy_cat) %>%
  dplyr::filter(!str_detect(term, "(Intercept)|PC|Age|Sex")) %>%
  dplyr::select(-data, -model_tsh, -statistic) %>%
  inner_join(repSNPs[c("Locus", "SNPid")],
             .,
             by = "SNPid") %>%
  mutate(associated = ifelse(p.value < 0.05/11, "Yes", "No")) %>% View
#filter(associated == "Yes") %>% View
filter(str_detect(term, "Dosage:")) %>% View
#write.csv("01-Feb-23_Kideny interaction with centralized thyroid (SNP-TSH)_long format.csv", row.names = FALSE, quote = FALSE)
pivot_wider(id_cols = c(Locus, SNPid),
            names_from = term,
            names_glue = "{term}_{.value}",
            values_from = c(estimate,
                            std.error,
                            #statistic,
                            p.value,
                            associated))


#---------#
write.csv(results_Kidney_TSH_cen, file = out_tbl_centralized, row.names = FALSE, quote = FALSE)




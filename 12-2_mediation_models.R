#=========================================#
#    Mediation Analysis: Running Steps
#=========================================#

library(tidyverse)
library(stringr)
library(purrr)


#-----------------------------------------------------#
#---------         phenom-wide scan        -----------
#-----------------------------------------------------#

# Mediation analysis steps are as follows:

# step 1: ln(eGFR) ~ SNP + kinship_matrix via emmax model in EPACTS
# step 2: ln(eGFR) ~ SNP + trait    + PC1:10
# step 3:    trait ~ SNP + PC1:10
# step 3:    trait ~ SNP + ln(eGFR) + PC1:10

#------------#
# Storing the variants name (SNPid)
targets <- vcfReg %>% select(starts_with("chr")) %>% colnames()


leadingSNPs <- c("chr1:10670853",  "chr2:15648568", 
                 "chr4:76492991",  "chr5:39385539",
                 "chr5:177386403", "chr7:77733187", 
                 "chr8:23894869",  "chr9:68819791",
                 "chr11:78339803", "chr15:98729803")


CKDGenSNPs <- c("chr1:10670853",  "chr2:15642347",
                "chr4:76480299",  "chr5:39377763",
                "chr5:177386403", "chr7:77793266",
                "chr8:23890907",  "chr9:68817258",
                "chr11:78312060", "chr15:98733292")



#-----------------------------------------------------#
#-------------- Step 1: association with eGFR --------
#-----------------------------------------------------#

# Step 1 results is indeed the results of GWAS and
# the replication of CKDGen EA in CHRIS 10K done on the servers.
# Commands are available here:
# /home/dghasemisemeskandeh/projects/gwas/replicationAnalysis/scripts/
# The summary or final results of replication analysis can be found in: 
# 



#-----------------------------------------------------#
#-------------- Step 2: trait as covariate -----------
#-----------------------------------------------------#


#------------#
getCoefs2table <- function(mytrait, mytarget, myformula, data){
    map2_df(.x = mytrait,
            .y = myformula,
            function(mytrait, myformula)
            {
              trait <- data[, mytrait]
              res0  <- map_df(mytarget,
                              function(mySNP)
                              {
                                SNP <- data[, mySNP]
                                myformula <- as.formula(myformula)
                                m <- lm(myformula, data = data)
                                p <- summary(m)$coefficients[2, c(1,2,4)] #Beta, se, Pvalue
                                names(p) <- str_replace_all(names(p),
                                                            c("Estimate"   = "Estimate",
                                                              "Std. Error" = "SE",
                                                              "Pr\\([^\\(]*\\)" = "Pvalue"))
                                p$SNPid <- mySNP
                                return(p)
                              })
              
              res0$Trait <- mytrait
              res1 <- 
                repSNPs %>%
                select(SNPid, Locus) %>% 
                inner_join(res0,
                           by = c("SNPid")) %>%
                select(SNPid, Locus, Trait, everything())
              return(res1)
            })
  }

#------------#
#outlier detection function for SNPs effects of eGFRw.res.ln
is_outlier <- function(x) { 
  return(x < quantile(x, 0.10) - 1.5 * IQR(x) | x > quantile(x, 0.90) + 1.5 * IQR(x))#x > quantile(x, 0.75) + 1.5 * IQR(x)
}
#------------#
# preparing the results of step 2 (covariate)
# for merging with step3 (direct)

# Forming principal components terms
PCs <- paste0("PC", 1:10, collapse = " + ")

# Step 2: long format
results_Step2_long <- 
  getCoefs2table(quantVars,
                 targets,
                 paste("eGFRw.log.Res ~ SNP + trait +", PCs),
                 vcfReg) %>%
  group_by(SNPid) %>%
  mutate(outlier = ifelse(is_outlier(Estimate), "Yes", "No")) %>% 
  ungroup()



#-----------------------------------------------------#
#--------------- Step 3a: trait as outcome -----------
#-----------------------------------------------------#

# Mediation Analysis -> Step 3: trait as the outcome
# The direct association between the SNP and quantitative traits

#------------#

#Example code: function for iterating
map_dfr(.x = vcfReg_TSHmod[targets[1:3]], 
        function(SNP){
          m <- lm(TSH.q ~ SNP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = vcfReg_TSHmod)
          p <- summary(m)$coefficients[2,c(1,4)]
          return(p) })
#------------#

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
    res4 <- repSNPs %>% select(SNPid, Locus) %>% inner_join(res3, by = "SNPid")
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
#------------#

# Total yes_outlier in Step2 = 1,135; excluding SCr = 973

#------------#
#Step 3 results: Wide format

# 1st way: Reshape the raw results
#results_Step3_wide <- 
results_Step3_long %>%
  pivot_wider(id_cols     = c(SNPid, Locus),
              names_from  = pheno,
              values_from = c(Estimate, SE , Pvalue, related),
              names_glue  = "{pheno}_{.value}") #%>% #View()
#write.csv("11-Jan-23_Step3_SNPs association with health traits.csv", quote = F, row.names = F)

#------------#
# 2st way: Use the latest version of function
results_Step3_wide <- step3_to_table(
  vcfReg[quantVars],
  vcfReg[targets],
  paste("trait ~ SNP +", PCs),
  vcfReg)

#------------#
results_Step3_long %>% 
  #count(Related)
  filter(related == "Yes") %>% 
  View()
#write.csv(., "12-Feb-22_Step2_significantly associated SNPs with health traits.csv", row.names = FALSE)



#-----------------------------------------------------#
#-------------- Step4: eGFR as mediator --------------
#-----------------------------------------------------#

#Step4: testing eGFR as a mediator for each trait
results_Step4_long <- getCoefs3table(
  vcfReg[quantVars],
  vcfReg[targets],
  paste("trait ~ SNP + eGFRw.log.Res +", PCs),
  vcfReg)



#-----------------------------------------------------#
#----------------- Mediation Analysis ----------------
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
  select(Locus, SNPid, RSID, EA_OA, ends_with("GWAS"), everything()) %>% #View()
  filter(Mediator == "Yes",
         Trait    != "SCr") %>% as_tibble() %>% View()
# write.csv("11-Jan-23_Table 3_Step1&2&3&4_6 SNPs successfully passed both steps.csv", row.names = F, quote = F)


write.csv(sum3steps_long,
          "11-Jan-23_SNP-wise summary of mediation analysis steps.csv",
          row.names = FALSE, quote = FALSE)

#---------#
#---------#
# Finding the frequency of SNPs of which 
# trait led to an outlier SNP effect size

#2nd try
sum3steps_long %>%
  group_by(Trait, related) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(id_cols = Trait,
              names_from = related,
              values_from = count,
              names_prefix = "related_") %>%
  mutate_all(~replace_na(., 0)) %>%
  arrange(factor(Trait, levels = quantVars)) %>%
  write.csv("11-Jan-23_Related traits in step 3 grouped by trait.csv", row.names = FALSE)

# The same table for step 2
sum3steps_long %>%
  group_by(Trait, outlier) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(id_cols = Trait,
              names_from = outlier,
              values_from = count,
              names_prefix = "outlier_") %>%
  mutate_all(~replace_na(., 0)) %>%
  arrange(fct_reorder(factor(Trait), outlier_Yes,.desc = T)) %>% View

#---------#
#1st try
sum3steps_outlier %>%
  select(#-ends_with("Step2"),
    -SNPid,
    -Locus) %>% 
  map(.x = ., .f = function(x) table(x)) %>% 
  bind_rows() %>% 
  mutate(variable = names(sum3steps_outlier)[-c(1,2)], #69:134, 3:68
         trait    = str_split(variable, "_Step[2-3]$", simplify = TRUE)[,1],
         step     = str_extract(variable, "Step[2-3]$")) %>% 
  #trait    = str_split(trait, "\\W", simplify = TRUE)[,1],
  #trait    = str_replace_all(trait, "\\_$", ""),
  #step     = replace(step, is.na(step), "Step3"),
  # step     = str_replace(step, "_Pvalue", "Step2")
  #mutate(No_Step2 = replace(No_Step2, is.na(No_Step2), 0))
  mutate_all(~replace_na(., 0)) %>% 
  select(-variable) %>%
  pivot_wider(names_from  = step, values_from = c(Mediator, No)) %>%  View()
#write.csv("31-May-22_Summary_of_3_steps_Mediation_Analysis_oulier beta P10-P90_step3.csv", row.names = FALSE, quote = FALSE)
ggplot(aes(x = trait, y = Mediator)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           mapping = aes(x = reorder(trait, Mediator), y = Mediator),
           nshow.legend = FALSE,
           width = 0.7,
           fill = "steelblue2",
           color = "grey50") +
  geom_text(aes(label = Mediator),
            hjust = -.3,
            color = "Black", 
            size = 2) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face="bold"),
        axis.text.y  = element_text(size = 5, face="bold")) +
  labs(x = NULL,
       y = "Number of SNPs whom traits led to an outlier effect size for Dosage (P10-P90)")+
  #y = "Number of SNPs whom traits led to a nonsignificant effect size for Dosage",
  #y = "Number of SNPs with significant effect size for Dosage associated with traits") +
  coord_flip()

ggsave("31-May-22_Outlier Beta (P10-P90) for Step3.png", last_plot(), width = 8, height = 5.5, pointsize = 0.5, dpi = 600, units = "in")

#---------#
pdf('27-May-22_DenstiyPlot_4features.pdf', width=15)

ggplot(vcfReg, aes(Hemoglobin)) + #pH, UrinaryGlucose, Proteins, Hemoglobin
  geom_density(color = "steelblue", fill = "gold2", alpha=0.9) +
  theme_minimal()

dev.off()


#-----------------------------------------------------------------#

#---------#
resAllModels2 <- 
  resAllModels %>%
  select(!contains(c("Model",
                     "GWAS",
                     "eGFR",
                     "eGFRw", 
                     "Age", 
                     "Proteins",
                     "Pvalue"))) #%>% View()
#%>% rename(eGFR.Beta = Model0.Beta, eGFR.SE = Model0.SE, eGFR.Pvalue = Model0.Pvalue)

#---------#

write.csv(resAllModels2, "15-Apr-2022_resAllModels_phase3_Coefs+Pvalues.csv", row.names = FALSE, quote = FALSE)

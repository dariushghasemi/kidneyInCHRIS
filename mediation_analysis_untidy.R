#=========================================#
#         Mediation Analysis: Full Script
#=========================================#

library(tidyverse)
library(stringr)
library(purrr)
#------------#



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

# SNPid for repSNPs file

#gsub("[a-zA-Z ]", "", vcf$MARKER_ID)
#str_extract(vcf$MARKER_ID, perl('(?<=[A-Z])\\d+'))
#vcf$SNPid <- noquote(str_split(vcf$MARKER_ID, "\\W", simplify=TRUE)[,6])
#repSNPs$SNPid <- as.character(noquote(str_split(repSNPs$MARKER_ID, "_", simplify=TRUE)[,1]))
#------------#

# Genotypes file in wide format
vcfmod <- 
  vcf %>%
  #mutate_at("MARKER_ID", str_replace, "chr", "")%>%
  mutate_at("MARKER_ID", 
            str_replace_all, 
            ":[A-T-C-G]+:[A-T-C-G]+", "")%>%
  rename(SNPid = MARKER_ID) %>% #head()
  pivot_wider(names_from = "SNPid",
              values_from= "Dosage")%>% 
  select(-c("chr11:78335892", "chr11:78392251"))
#-----------------------------------------------------#

# Combining dosages + kidney-thyroid traits + PCs
# This dataset is being used for mediation analysis.
# Later we will used this data to filter cases to be 
# used in interaction analyses with TSH and thyroid disease.

vcfReg <-
  chris[c("AID",
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



#-----------------------------------------------------#
#----------------- phenom-wide scan -------------------
#-----------------------------------------------------#

# Mediation analysis steps are as follows:

# step 1: ln(eGFR) ~ SNP + kinship_matrix via emmax model in EPACTS
# step 2: ln(eGFR) ~ SNP + trait    + PC1:10
# step 3:    trait ~ SNP + PC1:10
# step 3:    trait ~ SNP + ln(eGFR) + PC1:10


# Storing the variants name (SNPid)
targets <- vcfReg %>% select(starts_with("chr")) %>% colnames()



#-----------------------------------------------------#
#------------------- Testing model -------------------
#-----------------------------------------------------#

# Testing LR model
summary(lm(eGFRw.log.Res ~ eGFR + `chr1:10599281` + Operator + Municipality , data=vcfReg))#+ HbA1c
#------------#

# For Alexander to be replicated in SHIP Study
# Default seed for consistency of the results
set.seed(0)
#------------#

summary(lm(eGFRw.log.Res ~ `chr1:10599281` * TSH_cat + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = vcfReg_TSHmod))
summary(lm(eGFRw.log.Res ~ `chr1:10670853` + T3      + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = vcfReg_TSHmod))

#------------#
leadingSNPs <- c("chr1:10670853",  "chr2:15648568", 
                 "chr4:76492991",  "chr5:39385539",
                 "chr5:177386403", "chr7:77733187", 
                 "chr8:23894869",  "chr9:68819791",
                 "chr11:78339803", "chr15:98729803")
#------------#

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

# getCoefs2table <- 
#   function(mytrait, mytarget){
#     results1 <- lapply(mytrait,
#                        function(trait){
#                          map_df(mytarget,
#                                 function(SNP){
#                                   myformula <- as.formula(eGFRw.log.Res ~ SNP + trait + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
#                                   m <- lm(myformula, data = vcfReg)
#                                   s <- coef(m)[2]
#                                   p <- summary(m)$coefficients[2, c(1,2,4)] #Beta, se, Pvalue
#                                   #t <- tidy(m)[2, c(1,2,4)]#broom
#               })
#   })
#   results2 <- as.data.frame(do.call(cbind, results1)) #%>% clean_names() #library(janitor)
#   #names(results2) <- gsub(x = names(results2), pattern = "Pr\\([^\\(]*\\)", replacement = "_")
#   names(results2) <- str_replace_all(names(results2), c("Estimate"="Beta","Std. Error"="SE","Pr\\([^\\(]*\\)"="Pvalue"))
#   #colnames(results2) <- names(mytrait)
#   results3 <- cbind(SNPid = names(mytarget), results2)
#   return(results3)
#   }
#------------#

#Test the above function
getCoefs2table(vcfReg[quantVars[1:3]],
               vcfReg[targets[1:3]])

#------------#
getCoefs2table <-
  function(mytrait,
           mytarget,
           myformula,
           data){
    #    res1 <- 
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
                                #t <- tidy(m)[2, c(1,2,4)]#broom   OR   #%>% clean_names() #library(janitor)
                                #p$snpid <- mySNP  append(mySNP, p)
                                names(p) <- str_replace_all(names(p),
                                                            c("Estimate"   = "Estimate",
                                                              "Std. Error" = "SE",
                                                              "Pr\\([^\\(]*\\)" = "Pvalue"))
                                p$SNPid <- mySNP
                                return(p)
                              })
              #res1 <- res0
              res0$Trait <- mytrait
              res1 <- 
                repSNPs %>%
                select(SNPid, Locus) %>% 
                inner_join(res0,
                           by = c("SNPid")) %>%
                select(SNPid, Locus, Trait, everything())
              return(res1)
            })
    # nloc <- length(mytarget)
    # res3 <- cbind(SNPid = names(mytarget),
    #                   Locus = repSNPs$Locus[1:nloc],
    #                   res2)
    #return(res0)
  }

#Test the above function
getCoefs2table(quantVars[1:4],
               targets[1:3],
               paste("eGFRw.log.Res ~ SNP + trait + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
               vcfReg)

cat(paste0("PC",1:9," +"), "PC10")

#------------#
#preparing the results of step 2 (covariate) for merging with step3 (direct)

# Step 2: long format
results_Step2_long <- 
  getCoefs2table(quantVars,
                 targets,
                 paste("eGFRw.log.Res ~ SNP + trait + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
                 vcfReg) %>%
  group_by(SNPid) %>%
  mutate(outlier = ifelse(is_outlier(Estimate),
                          "Yes",
                          "No")) %>% 
  ungroup()

#------------#
# Step 2: wide format

#results_Step2_wide <- 
results_Step2_long %>% 
  pivot_wider(id_cols     = c(SNPid, Locus),
              names_from  = Trait,
              values_from = c(Estimate, SE , Pvalue, outlier),
              names_glue  = "{Trait}_{.value}") %>%
  select(matches("[_](SE|Pvalue)")) %>% View()
#relocate( ends_with(c("Estimate", "SE", "Pvalue"))) 
#write.csv(., "27-Dec-22_Step2_SNPs association with log(eGFR) adjusted for health traits.csv", quote = F , row.names = FALSE)

#------------#
# Step 2: wide to long format
resModelN %>%
  pivot_longer(cols = !SNPid,   
               names_to = c("trait", "value"),
               names_pattern = "(.+).(Beta|SE|Pvalue)$", #([A-Za-z0-9_]+)
               #names_sep = '.',
               values_to = c("score")) %>% 
  pivot_wider(names_from = "value",
              values_from = "score") %>% View

#------------#

#Outlier traits grouped by Locus and SNPid
results_Step2_long %>% 
  #mutate(outlierTrait = ifelse(outlier == "Yes", trait, "")) %>%
  filter(outlier == "Yes") %>%
  group_by(Locus) %>%
  count(Trait) %>% View()
#write.csv("08-Oct-2022_outlier traits in Step2 grouped by trait.csv", row.names = FALSE)
ggplot(aes(y = forcats::fct_rev(forcats::fct_infreq(Trait)), fill = Locus)) +
  geom_bar(width = .85,  position = "dodge")+
  scale_x_continuous(breaks = seq(1,10, 1), limits = c(0,10.2), expand = c(0,0)) +
  #scale_fill_brewer(palette = "Paired")+#"Accent"
  #scale_fill_viridis_d(option = "magma")+
  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.y  = element_text(size = 10,  face="bold"),
        legend.position = c(.88, .4),
        legend.key.size = unit(0.8, 'cm'),
        panel.grid.major.x = element_line(color = "lightgray", size = .25)) +
  labs(x = "Frequency of the SNPs with outlier trait in step 2", y = NULL)

ggsave("08-Oct-22_Frequency of the SNPs with outlier trait in step 2.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")



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

getCoefs3table <- 
  function(mytrait, 
           mytarget, 
           myformula,
           mydata){
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

#Step 3 results: Long format
results_Step3_long <-
  getCoefs3table(vcfReg[quantVars],
                 vcfReg[targets],
                 paste("trait ~ SNP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
                 vcfReg) %>%
  mutate(related = ifelse(Pvalue <= 0.05/770,
                          "Yes",
                          "No"))
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
results_Step3_wide <-
  step3_to_table(vcfReg[quantVars],
                 vcfReg[targets],
                 paste("trait ~ SNP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
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
results_Step4_long <-
  getCoefs3table(vcfReg[quantVars],
                 vcfReg[targets],
                 paste("trait ~ SNP + eGFRw.log.Res + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
                 vcfReg)



#-----------------------------------------------------#
#----------------- Mediation Analysis ----------------
#-----------------------------------------------------#

# Collecting the 3 Steps of the Mediation Analysis

#step 1: ln(eGFR) ~ SNP + PC1:10
#step 2: ln(eGFR) ~ SNP + trait    + PC1:10
#step 3:    trait ~ SNP + PC1:10
#step 3:    trait ~ SNP + ln(eGFR) + PC1:10


sum3steps_long <-
  repSNPs %>%
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


write.csv(sum3steps_long, "11-Jan-23_SNP-wise summary of mediation analysis steps.csv", row.names = FALSE, quote = FALSE)







#---------#
#Saving model result containing Beta, SE, Pvalue, and outlier Yes/No
results_Step2 <-
  results_Step2_long %>%
  pivot_wider(id_cols = c(SNPid, Locus),
              names_from = trait,
              values_from = c(Beta, SE , Pvalue, outlier),
              names_glue = "{trait}_{.value}") %>% 
  #select(contains(c("SNPid", "Locus", "T3", "T4"))) %>% 
  write.csv(., "25-Aug-2022_phase3_OutlierTraits_groupedBy_Locus_SNPid+Beta+SE+Pvalue.csv", row.names = FALSE, quote = FALSE)

#---------#
#Merging Step 2 and Step 3 result
results_Step2 %>% 
  inner_join(results_Step3,
             by = c("SNPid", "Locus")) %>%
  #write.csv(., "12-Aug-2022_results_step2-outlier_and_step3-related_Beta+SE+Pvalue.csv", row.names = FALSE)
  select(contains(c("SNPid",
                    "Locus",
                    "outlier",
                    "related"))) #%>%
#write.csv(., "12-Aug-2022_results_step2-outlier_and_step3-related.csv", row.names = FALSE)

#---------#
#sum3steps_outlier <-
select(-ends_with(c("Estimate", "Pvalue" , "SE")),#"Beta"
       -starts_with(c("Model",
                      "GWAS",
                      "Calcium_Corrected_mg",
                      "T_Iron_Binding",
                      "Transferrin_saturation",
                      "AntiTPO",
                      "Specific_weight"))) %>%
  #select_all(~gsub("\\_Pvalue$", "_Step3", .)) %>% #"\\_Estimate$"
  select_all(~gsub("\\.Beta$", "_Step2", .)) %>% #"\\.Pvalue$"
  rename_with(~str_replace(.x,
                           "TSH.q",
                           "TSH"), matches("TSH")) %>%
  pivot_longer(cols      = ends_with("Step2"),
               names_to  = "feature",
               values_to = "Beta") %>%
  group_by(SNPid, Locus) %>%   #select(contains("T3")) %>% View()
  mutate(Outlier = if_else(is_outlier(Beta), "Mediator", "No")) %>%
  ungroup() %>%
  #filter(SNPid == "chr1:10670853" & Locus == "CASZ1" | Outlier == "Mediator")
  pivot_wider(id_cols     = c(SNPid, Locus),
              names_from  = feature,
              values_from = Outlier) #%>% View()
# mutate(SNPid, Locus,
#        across(ends_with("Step2"),
#               ~ if_else(.x >= 0.05,
#               #~ if_else(is_outlier(.x),
#                         "Mediator",
#                         "No"))) %>% #, .names = "{.col}.fn{.fn}")
# mutate(across(ends_with("Step3"),
#               ~ if_else(.x <= 0.05/620,
#                         "Mediator",
#                         "No"))) #%>% View()
# map(.x = ., .f = function(x) #if_else(x >= 0.05,
#   if_else(is_outlier(x), "Mediator", "No"))
#cbind(resAllModels[c("SNPid")], .) %>% View()

#---------#
# Finding the frequency of SNPs of which trait led to an outlier SNP effect size

#2st try
sum3steps_long %>%
  group_by(Trait, related) %>%
  #map("outlier", ~ count(.x))
  #summarise( n.out = n()) %>% View()
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(id_cols = Trait,
              names_from = related,
              values_from = count,
              names_prefix = "related_") %>%
  mutate_all(~replace_na(., 0)) %>%
  arrange(factor(Trait, levels = quantVars)) %>%
  #slice(order(factor(Trait, levels = quantVars)))%>%
  #left_join(data.frame("Trait" = quantVars), ., by = "Trait")
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



#-----------------------------------------------------#
#---------------- Figure 4: line range ---------------
#-----------------------------------------------------#

#New line range plot for the paper
#Changing the order and the label of the Locus top SNPs
Locus_factor <- function(x) {
  factor(x,
         levels = c("CASZ1","DDX1","SHROOM3","DAB2","SLC34A1",
                    "TMEM60","STC1","PIP5K1B","GAB2","IGF1R","PDILT"),
         labels = c("CASZ1\n1:10670853",
                    "DDX1\n2:15642347",
                    "SHROOM3\n4:76480299",
                    "DAB2\n5:39393631",
                    "SLC34A1\n5:177386403",
                    "TMEM60\n7:77714744",
                    "STC1\n8:23885208",
                    "PIP5K1B\n9:68817258",
                    "GAB2\n11:78324786",
                    "IGF1R\n15:98733292",
                    "PDILT\n16:20381010"))
}

# with rsID
Locus_factor <- function(x) {
  factor(x,
         levels = c("CASZ1","DDX1","SHROOM3","DAB2","SLC34A1",
                    "TMEM60","STC1","PIP5K1B","GAB2","IGF1R","PDILT"),
         labels = c("CASZ1\nrs74748843",
                    "DDX1\nrs807624",
                    "SHROOM3\nrs28817415",
                    "DAB2\nrs10062079",
                    "SLC34A1\nrs3812036",
                    "TMEM60\nrs57514204",
                    "STC1\nrs819196",
                    "PIP5K1B\nrs2039424",
                    "GAB2\nrs7113042",
                    "IGF1R\nrs59646751",
                    "PDILT\nrs77924615"))
}

#---------#
#Taking effect size of the top SNPs in GWAS by subsetting summary of Mediation-A
sum3steps_long_eGFR <-
  sum3steps_long %>%
  filter(SNPid %in% tagSNPs$SNPid,
         Trait == "SCr") %>% 
  mutate(Locus_reordered = Locus_factor(Locus)) #%>% View()

#---------#
# Reorder Trait based on frequency of outlier effect for the lead SNPs
results_Step2_long_freq <- 
  results_Step2_long %>%
  filter(SNPid %in% tagSNPs$SNPid) %>%
  group_by(Trait, outlier) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(id_cols = Trait,
              names_from = outlier,
              values_from = count,
              names_prefix = "outlier_") %>%
  mutate_all(~replace_na(., 0)) %>%
  mutate(Feature_freq = fct_reorder(factor(Trait), outlier_Yes, .desc = T)) %>%
  arrange(Feature_freq) %>% View

#---------#
#Step2: Line range plot of the CKDGen lead variants for the paper

results_Step2_long %>%
  filter(SNPid %in% tagSNPs$SNPid) %>%
  right_join(results_Step2_long_freq, by = "Trait") %>%
  mutate(Trait = recode(Trait,
                        "Body_Fat"   = "Body Fat",
                        "Visceral_Fat" = "Visceral Fat",
                        "Pulse_Rate" = "Pulse Rate",
                        "INR_PT"     = "INR PT",
                        "APTT"       = "aPTT",
                        "APTT_ratio" = "aPTT ratio",
                        "AST_GOT"    = "AST GPT",
                        "ALT_GPT"    = "ALT GPT",
                        "Urine_pH"   = "Urine pH"),
         Feature         = fct_reorder(Trait, Estimate,    .desc = F),
         Feature_freq    = fct_reorder(Trait, outlier_Yes, .desc = F),
         Locus_reordered = Locus_factor(Locus)) %>%
  ggplot(aes(Estimate,
             Feature_freq,
             xmin = Estimate - SE,
             xmax = Estimate + SE)) +
  # geom_vline(data = sum3steps_long_eGFR,
  #            aes(xintercept = Estimate_GWAS), 
  #            color = "steelblue3",
  #            lty = 1) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50")+
  geom_pointrange(aes(color = outlier), show.legend = F, fatten = 1.2) +
  scale_color_manual(values = c("grey50", "tomato3")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  coord_cartesian() +
  facet_wrap(~ Locus_reordered, nrow = 1, scales = "free_x", shrink = T)+
  labs(x = "Trait-adjusted effect of CKDGen lead SNP at each locus") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside",
        axis.text.x = element_text(size = 6,  face = "bold"),
        axis.text.y = element_text(size = 6,  face = "bold"),
        axis.title  = element_text(size = 10, face = "bold"),
        #axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray", size = .25),
        strip.text.x = element_text(size = 8, color = "Black", face = "bold.italic"))

#ggsave("Top_SNPs_LinerangePlot_CASZ1_10-Apr-22.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

ggsave("28-Mar-23_Linerange plot for the leading SNP in each locus.png",
       last_plot(), width = 12, height = 7, pointsize = 2, dpi = 300, units = "in")

#---------#

#Step2: Line range plot of the entire 163 variants
pdf('12-Jan-23_Linerange plot of the entire 163 replicated SNPs.pdf', width=8, height = 8)

map(1:length(targets), function(i){
  results_Step2_long %>%
    filter(SNPid %in% targets[i]) %>%
    mutate(Trait = recode(Trait,
                          "Body_Fat"   = "Body Fat",
                          "Visceral_Fat" = "Visceral Fat",
                          "Pulse_Rate" = "Pulse Rate",
                          "INR_PT"     = "INR PT",
                          "APTT_ratio" = "APTT ratio",
                          "AST_GOT"    = "AST GPT",
                          "ALT_GPT"    = "ALT GPT",
                          "Urine_pH"   = "Urine pH"),
           Locus_reordered = Locus_factor(Locus),
           Feature = fct_reorder(Trait, Estimate, .desc = F)) %>%
    ggplot(aes(Estimate, 
               Feature,
               xmin = Estimate - SE, 
               xmax = Estimate + SE)) +
    geom_vline(data = repSNPs[i,],
               aes(xintercept = Beta_CHRIS),
               color = "steelblue3",
               lty = 1) +
    geom_vline(xintercept = 0, lty = 2, color = "grey50")+
    geom_pointrange(aes(color = outlier), show.legend = F, fatten = 1.2) +
    scale_color_manual(values = c("grey50", "tomato3")) +
    #scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    coord_cartesian() +
    labs(x = paste0("Trait-adjusted effect of the ", targets[i], " SNP in *", repSNPs[i,"Locus"], "* locus on ln(eGFRcreat)")) + # in SHROOM3 Locus
    theme_classic() +
    theme(strip.background = element_blank(), strip.placement = "outside",
          axis.text.x = element_text(size = 6,  face = "bold"),
          axis.text.y = element_text(size = 6,  face = "bold"),
          axis.title  = element_text(size = 10, face = "bold"),
          axis.title.x = ggtext::element_markdown(),
          axis.title.y.left = element_blank(),
          panel.grid.major.y = element_line(color = "lightgray", size = .25),
          strip.text.x = element_text(size = 8, color = "Black", face = "bold.italic"))
}
)

dev.off()


#---------#
#-----------------------------------------------------------------#

#outlier detection function for SNPs effects of eGFRw.res.ln
is_outlier <- function(x) { 
  return(x < quantile(x, 0.10) - 1.5 * IQR(x) | x > quantile(x, 0.90) + 1.5 * IQR(x))#x > quantile(x, 0.75) + 1.5 * IQR(x)
}
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
#---------#
resAllModels2 %>%
  gather(key = "feature", value = "coef", -SNPid, -Locus, na.rm = FALSE) %>% 
  group_by(Locus) %>% 
  mutate(outlier = ifelse(is_outlier(coef), feature, NA)) %>% 
  ungroup() %>%
  ggplot(aes(Locus, coef))+
  geom_boxplot(alpha = 0.6 , outlier.color = "red") +
  geom_text(aes(label = outlier, color = outlier), 
            show.legend = FALSE, na.rm = TRUE, size = 2, vjust = -0.5,
            position = position_jitter(width = -0.1)) +
  theme_classic()

#ggsave("Top_SNPs_Boxplot.png", ., width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#-----------------------------------------------------------------#
#---------#

#Box plot for paper
library(ggtext)
library(ggrepel)

resAllModels3_162_grouped <-
  resAllModels2_Beta_162 %>% 
  inner_join(resAllModels2_SE_162, 
             by = c("SNPid", "Locus", "trait")) %>% 
  select(-var1, -var2) %>%
  mutate(trait = replace(trait, trait == "SerumCreatinine", "SCr"),
         trait = replace(trait, trait == "Uric_Acid",       "Urate"),
         trait = replace(trait, trait == "PTT_sec",         "PTT"),
         trait = replace(trait, trait == "UrinaryAlbumin",  "UAlb"),
         trait = replace(trait, trait == "Homocysteine",    "Homocyst"),
         trait = replace(trait, trait == "Magnesium_mg",    "Magnesium"),
         trait = replace(trait, trait == "Calcium_mg",      "Calcium"),
         trait = replace(trait, trait == "pH",              "Urine pH"),
         trait = replace(trait, trait == "BodyFat",         "Body Fat"),
         trait = replace(trait, trait == "PulseRate",       "Pulse Rate"),
         trait = replace(trait, trait == "BPsystolic",      "SBP"),
         trait = replace(trait, trait == "BPdiastolic",     "DBP")) %>% #head()
  group_by(Locus, SNPid) %>% View()
mutate(outlier = ifelse(is_outlier(Beta),
                        "Yes", #trait,
                        "No")) %>% #NA #"All Other Traits"
  ungroup() #%>% 
#filter(SNPid == "chr11:78319000") %>%
#arrange(desc(Beta)) %>%
#View()

#---------#
#resAllModels3_162_grouped %>% count(outlier)
# A tibble: 2 x 2
#outlier     n
#1 No      11,177
#2 Yes      1,135

#---------#
#---------#
#---------#



#-----------------------------------------------------#
#-------------------- Heatmap plots ------------------
#-----------------------------------------------------#

#Heatmap of step2 and step3
mediatorSNPs <- 
  results_Step3_long %>%
  inner_join(results_Step2_long,
             by = c("SNPid", "Locus", "trait"),
             suffix = c("_step3", "_step2")) %>% 
  mutate(Mediator = case_when(
    outlier == "Yes" & related == "Yes" ~ "Yes/Yes",
    outlier == "Yes" & related != "Yes" ~ "Yes/No",
    outlier != "Yes" & related == "Yes" ~ "No/Yes",
    outlier != "Yes" & related != "Yes" ~ "No/No"),
    Moderator = ifelse(outlier == "Yes" & related == "Yes", "Yes", "No"),
    SNP = factor(SNPid,
                 levels = str_sort(unique(SNPid),
                                   numeric = TRUE,
                                   decreasing = TRUE))) %>% 
  filter(#Mediator != "No/No" & Mediator != "No/Yes"
    #Mediator == "Yes/Yes"
    Moderator == "Yes",
    trait     != "SCr") %>% #View()
  #count(Mediator)#outlier, related)
  inner_join(repSNPs[c("SNPid",
                       "Locus",
                       "MARKER_ID",
                       "BETA",
                       "SEBETA",
                       "PVALUE")],
             by = c("SNPid", "Locus")) %>%
  rename(Beta_GWAS = BETA, SE_GWAS = SEBETA, Pvalue_GWAS = PVALUE) %>%
  select(SNPid, Locus, MARKER_ID, Beta_GWAS, SE_GWAS, Pvalue_GWAS, everything()) %>% View()
#write.csv(., "25-Feb-2022_Step2&3_12 SNPs successfully passed both steps_12SNPs.csv", row.names = FALSE)
ggplot(aes(x = trait, y = SNP, fill = Moderator)) +
  geom_tile() +
  theme_classic()+
  labs(x    = "Health Trait", #changing the legend title
       y    = "SNPid")+#, fill = "Mediator\n(Outlier/Related)") +
  #scale_fill_grey(start = 0.8, end = 0.2) + #'#56B4E9', "#E7B800"
  scale_fill_manual(values=c('#999999', "grey20" , "grey50", "#E69F00", "#FF6666"))+
  theme(legend.title= element_text(size=7, face="bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text  = element_text(size = 8),
        axis.text.x = element_text(size=6, face="bold", angle=35, vjust=1.05, hjust=1),
        axis.text.y = element_text(size=4, face="bold"),
        axis.title  = element_text(size=12,face="bold"))

ggsave("25-Aug-22_Heatmap of Step 2&3_Yes-No Mediator(Outlier_Related) traits_fliped3.png", 
       last_plot(), width = 10, height = 9, pointsize = 4, dpi = 300, units = "in")
#-----------------------------------------------------#

#Summary of outliers for Cristian deduction
pdf('10-Aug-22_Top_SNPs_LinerangePlot_GAB2_3SNPs.pdf', width=8, height = 8)

lapply(c("chr11:78315470", "chr11:78319000", "chr11:78319103"),
       function(pos) resAllModels3_162_grouped %>% 
         filter(SNPid == pos) %>% 
         arrange(desc(Beta)) %>% #View()
         mutate(Feature = fct_reorder(trait, Beta)) %>% 
         ggplot(aes(Beta, Feature, xmin = Beta - SE, xmax = Beta + SE))+
         geom_point(aes(color = outlier))+
         geom_linerange()+
         geom_vline(aes(xintercept = quantile(Beta, probs = c(.10))), lty = 2, size = .4, color = "indianred3") +
         geom_vline(aes(xintercept = quantile(Beta, probs = c(.90))), lty = 2, size = .4, color = "indianred3") +
         labs(x = paste("Adjusted effect size of the",
                        pos,
                        "SNP in GAB2 locus"), 
              y = "Health Trait")+
         theme_classic() +
         theme(axis.text.y = element_text(size=4, face="bold"),
               axis.title = element_text(size=10, face="bold"),
               legend.position = c(0.8, 0.3)))

dev.off()

#---------#
#Heatmap for the Yes/No outlier adj effect sizes

library(wesanderson)
resAllModels3_162_grouped %>%
  mutate(id = fct_reorder(SNPid, trait)) %>% #View()
  ggplot(aes(x = trait, y = SNPid, fill = outlier)) +
  geom_tile() +
  theme_classic() +
  labs(x = "Health Trait") +
  #scale_fill_grey(start = 0.8, end = 0.2) +
  scale_fill_manual(values=c('#56B4E9', '#999999'))+
  theme(axis.text.x = element_text(size=4, face="bold", angle=40),
        axis.text.y = element_text(size=4, face="bold"),
        axis.title  = element_text(size=10, face="bold"))

ggsave("10-Aug-22_Heatmap of the phase3_Yes-No Outlier Traits grouped by Locus and SNPid.png", 
       last_plot(), width = 12, height = 10, pointsize = 4, dpi = 300, units = "in")

#--------------------------------------------#
#subseting GWAS results
resAllModels3_eGFR <- 
  resAllModels3 %>% 
  filter(SNPid %in% tagSNPs$SNPid, trait == "eGFR")
resAllModels3_eGFR %>% View()
#---------#

dodge_with = 0.6

ggplot(data = resAllModels3,
       aes(x = Locus,
           y = Beta)) +
  geom_hline(aes(yintercept = 0),
             lty = 2,
             size = 1,
             color = "grey50") +
  geom_point(data = resAllModels3_eGFR,
             aes(x = Locus,
                 y = Beta),
             color = "black",
             size = 1.5,
             alpha = 0.9) +
  geom_linerange(data = resAllModels3_eGFR,
                 aes(ymin = Beta - 1.96*SE,
                     ymax = Beta + 1.96*SE),
                 color = "black",
                 size = 0.8,
                 lty = 1,
                 alpha = 0.8) +
  geom_segment(data = resAllModels3_eGFR,
               aes(x = 1:10 - dodge_with/2,
                   xend = 1:10 + dodge_with/2, y = Beta,
                   yend = Beta),
               size = 2,
               alpha = 0.8,
               lineend = "round",
               inherit.aes = FALSE) +
  geom_pointrange(data = resAllModels3,
                  aes(x = Locus,
                      y = Beta,
                      ymax = Beta + 1.96*SE,
                      ymin = Beta - 1.96*SE,
                      color = outlier), #, fill=trait
                  size = 0.25,
                  position = position_dodge(dodge_with),
                  show.legend = TRUE,
                  inherit.aes = FALSE) + 
  scale_x_discrete(limits = resAllModels3_eGFR$Locus,
                   breaks = resAllModels3_eGFR$Locus,
                   labels = paste0(resAllModels3_eGFR$Locus,
                                   "\n",
                                   resAllModels3_eGFR$SNPid)) +
  scale_color_manual(values=
                       c("grey50",      "gold3",       "darkorchid1",
                         "skyblue2",    "tomato",      "orange2",
                         "pink2",       "darkorchid3", "maroon1",
                         "brown2",      "orangered2",  "springgreen1",
                         "darkviolet",  "royalblue2",  "royalblue4",
                         "springgreen3","chartreuse3", "olivedrab3", "gold4"))+ #turquoise2
  coord_cartesian(ylim = c(-0.05, 0.05)) +
  labs(y = "Effect on log(eGFR)",
       fill = "Trait") +
  #scale_fill_discrete(name = "Title", labels = resAllModels3$outlier) +
  theme_classic() + 
  theme(legend.title = element_blank(),#legend.position = "top",
        axis.text.y  = element_text(size = 8,  face = "bold"),
        axis.text.x  = element_text(size = 7,  face = "bold"),
        axis.title   = element_text(size = 11, face = "bold"),
        legend.key.size = unit(0.6, 'cm'),
        legend.text  = element_text(size = 12))

#geom_boxplot(color = "Black", alpha = 0.6, outlier.shape = NA, position = "dodge2")+
#geom_point(position = position_dodge(width = 0.6), size = 1.2, alpha=0.7) +
#geom_linerange(aes(ymin = Beta - SE, ymax = Beta + SE), 
#               position = position_dodge(width = 0.6), size = 0.9, alpha=0.7, linetype = 1)+
#geom_text_repel(aes(color = trait, label = trait), position = position_jitter(width=c(0.4, 1.2)), size = 3, fontface=2)
##+ coord_flip()
#geom_text(aes(color=outlier, label=outlier), fill=NA, label.color=NA, position = position_dodge(c(0.6,1.2)), size=1.5, fontface=2)+
#+#outlier.color = "red",
#geom_jitter(aes(color = outlier), position = position_jitter(0.2), show.legend = FALSE, hjust="outward")+
#geom_text_repel(aes(color = outlier, label = outlier), position = position_jitter(0.1), size = 3, fontface=2, hjust="outward")+
#geom_text(aes(label = outlier, color = outlier), show.legend = FALSE, na.rm = TRUE, size = 2,position = position_jitter(0.15))+

ggsave("Top_SNPs_Boxplot_20-Apr-22_poster.png", last_plot(), width = 11, height = 5.5, pointsize = 0.5, dpi = 600, units = "in")

#---------#
resAllModels2 %>%
  select(TSH.Beta, T3.Beta, Locus) %>% 
  filter(Locus == "CASZ1") %>% 
  select(-Locus) %>%
  pivot_longer(cols = everything(),
               names_to = c("BETA", "set"),
               names_pattern = "(.).(.+)") #%>% View()
#pivot_longer(cols = contains("Beta"), names_pattern = ".", names_to = c("var1", "var2"), values_to = "coef") %>% View()
#gather(key = "feature" , value = "coef", -SNPid, -Locus, na.rm = FALSE) 

#Watch the below videos which are helpful to proceed the analysis
#CodeClub: CC0134 & CC0037
#Julia Silge: Episode 43 - Predicting housing prices in Austing TX 
#-----------------------------------------------------#

pdf('resAllModels_SNPsCoefs_Comparison_29Jan.pdf', width=8, height = 6)

lapply(c(2,4,6,8), function(i)
{
  ggplot(resAllModels, aes(BETA_GWAS, resAllModels[,i])) + 
    geom_abline(slope = 1, intercept = 0, color = "gray") + 
    geom_point(aes(shape = Locus, color = Locus), size = 2) +
    scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) + 
    ylab(names(resAllModels[i])) + theme_classic()})

#pdf('resAllModels_SNPsCoefs_Heatmaps_29Jan.pdf', width=18, height = 18)
pheatmap(resAllModels[-c(1:3)] / resAllModels$BETA,
         cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = T, labels_row = resAllModels$SNPid, 
         border_color = 'Black', fontsize_row = 6, fontsize_col = 6, angle_col = "45")

dev.off()











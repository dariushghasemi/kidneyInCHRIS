#=============================================#
# Interrogating replicated SNPs in databases
#=============================================#


library(MendelianRandomization)
library(gwasrapidd)

#-----------------------------------------------------#
#---------------- Mendelian Randomization ------------
#-----------------------------------------------------#

# To query the variants, we are using "MendelianRandomization" package.
# results of PhenoScanner database for genotype-phenotype associations
phenoscanner(snpquery = "rs807624")$results %>% View()

#
pheno_input(snps      = c("rs807624"),
            exposure  = "Serum creatinine",
            pmidE     = "20383146", 
            ancestryE = "European")
#exposure = "Low density lipoprotein", pmidE = "24097068", 
#ancestryE = "European", outcome = "eGFR", pmidO = "26343387", ancestryO = "Mixed")

# bmi_snps <- c("rs1558902", "rs2867125", "rs571312", "rs10938397", "rs10767664", "rs2815752", "rs7359397", "rs9816226",
#               "rs3817334", "rs29941", "rs543874", "rs987237", "rs7138803", "rs10150332", "rs713586", "rs12444979", "rs2241423",
#               "rs2287019", "rs1514175", "rs13107325", "rs2112347", "rs10968576", "rs3810291", "rs887912", "rs13078807",
#               "rs11847697", "rs2890652", "rs1555543", "rs4771122", "rs4836133", "rs4929949", "rs206936")

ps_result <- pheno_input(snps     = ps_object$rsid[1:100],
                         exposure = "Body mass index",
                         pmidE    = "UKBB", ancestryE = "European", 
                         outcome  = "Type II diabetes",
                         pmidO    = "28566273", ancestryO = "European")#%>% View()

BMISNPs <- read.table("C:\\Users\\dghasemisemeskandeh\\Desktop\\rsid.txt", header = TRUE)

mr_ivw(ps_result)
mr_plot(ps_result, orientate=TRUE)
mr_allmethods(ps_result, method="main")



#-----------------------------------------------------#
#-------------------- Phenoscanner -------------------
#-----------------------------------------------------#

# Looking up the 162 SNPs in the databases for finding 
# an evidence of the association with other traits
# To query the variants, we are using "MendelianRandomization" package.
#----------#

# looking up 100 + 62 + PDILT replicated SNPs + Mediator SNPs
ps_object0 <- phenoscanner(snpquery = repSNPs$SNPid[1:100],   build = 38, pvalue = 5E-8)
ps_object1 <- phenoscanner(snpquery = repSNPs$SNPid[101:162], build = 38, pvalue = 5E-8)
ps_object3 <- phenoscanner(snpquery = repSNPs$SNPid[163],     build = 38, pvalue = 5E-8)
ps_object2 <- phenoscanner(snpquery = mediatorSNPs$SNPid,     build = 38, pvalue = 5E-8)

#----------#
ps_object0$results %>%
  #filter(grepl("magnesium|Magnesium", trait)) 
  filter(
    str_detect(trait,
               regex("magnesium|Magnesium|PTT|aPTT|thromboplastin|Thromboplastin|coagulation|Coagulation"))) %>%
  write.csv("25-Aug-2022_GenomeWideSig SNPs associated with Coagulation and Magnesium queried using phenoscanner.csv",
            row.names = FALSE, quote=FALSE)

#----------#
# Frequency of SNPs with outlier trait-adjusted effect on kidney
results_Step2_long %>%
  filter(outlier == "Yes") %>%
  count(Trait) %>%
  write.csv("01-Feb-23_Summary of Step 2_Number of SNPs with outlier trait-adjusted effect.csv",
            row.names = FALSE, quote = FALSE)

#----------#
# Gathering phenoscanner look-up
results_step3_long_pscanner <-
  ps_object0$results %>%
  bind_rows(ps_object1$results) %>%
  bind_rows(ps_object3$results) %>%
  merge(repSNPs[c("Locus","SNPid")], ., by.x = "SNPid", by.y = "snp", sort = F) %>%
  mutate(p = as.numeric(p)) %>%
  rename(EA = a1, OA = a2) %>%
  select(-SNPid) %>%
  #select(snp, study, pmid,trait, beta, se, p) %>% 
  filter(p < 5e-08) %>% 
  #count(trait) %>% View
  #write.csv("31-Jan-23_163 SNPs queried using phenoscanner at genome-wide sig associated with trait.csv", row.names = FALSE, quote=FALSE)
  #add_count(trait, name = "N_of_gw_snp", sort = T) %>% View
  #arrange(desc(n)) %>% View()
  count(trait, name = "N_of_gw_snp", sort = T) %>% View
  #filter(n>10) %>%
  write.csv("31-Jan-23_163 SNPs queried using phenoscanner at genome-wide sig associated with trait_Frequency.csv", 
            row.names = FALSE, quote=FALSE)
  
  ggplot(aes(x = trait, y = n)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           mapping = aes(x = reorder(trait, n), y = n),
           nshow.legend = FALSE,
           width = 0.7,
           fill = "steelblue2",
           color = "grey50") +
  geom_text(aes(label = n),
            hjust = -.3,
            color = "Black", 
            size = 2) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face="bold"),
        axis.text.y  = element_text(size = 5,  face="bold")) +
  labs(x = NULL,
       y = "Frequency of the traits found in Phenoscanner")+
  #y = "Number of SNPs whom traits led to a nonsignificant effect size for Dosage",
  coord_flip()

ggsave("23-Jun-22_Frequency of the traits found in Phenoscanner at GWS level.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")



#-----------------------------------------------------#
#--------------------- OpenTargets -------------------
#-----------------------------------------------------#

# Reading the traits found to be related 
# with the 162 SNPs in OpenTargets by David

opentargets_David <- 
  read.table("Dariush_SNPs_OpenTargets_4_REPORT.txt", header = T, sep = "\t")

opentargets <- 
  repSNPs[,c("SNPid", "Locus")] %>% 
  right_join(opentargets_David, by = c("SNPid")) %>%
  filter(pval <= 5E-8) %>%
  mutate(trait = factor(traitReported)) #%>% #View()
#mutate(trait = fct_reorder(traitReported, n, .desc = TRUE))
#----------#

#opentargets_count <- 
opentargets %>%
  #group_by(Locus) %>%
  count(Locus, traitReported, sort = TRUE) %>% 
  #ungroup(Locus) %>% 
  #arrange(Locus, desc(n)) %>% #View()
  ggplot(aes(n, traitReported, fill = Locus)) +
  geom_bar(aes(n, fct_reorder(traitReported, n, .fun = sum, .desc = FALSE)),
           stat = "identity") +
  geom_text(aes(label = n),
            position = position_stack(vjust = .97),
            colour = "Black", size = 2)+
  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.y  = element_text(size = 5,  face="bold"),
        legend.position = c(.85, .4))+
  labs(x = "Frequency of the traits in OpenTargets", y = NULL)
#----------#

ggplot(opentargets, aes(y = forcats::fct_rev(forcats::fct_infreq(trait)), fill = Locus)) +
  geom_bar(width = .85) +
  # geom_text(opentargets_count, aes(label = n), 
  #           position = position_stack(vjust = 0.5),
  #           #hjust = -.3, 
  #           color = "Black",
  #           size  = 2) +
  scale_x_continuous(limits = c(0,80), expand = c(0,0)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face="bold"),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.y  = element_text(size = 5,  face="bold"),
        legend.position = c(.85, .4),
        panel.grid.major.x = element_line(color = "lightgray", size = .25)) +
  labs(x = "Frequency of the traits in OpenTargets", y = NULL)

ggplot(opentargets)+ #, fill = Locus
  geom_col(aes(reorder(trait, n), n, fill = Locus), position = "dodge")

ggsave("13-Aug-22_Frequency of the traits found in OpenTargets at GWS level grouped by Locus4.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")



#-----------------------------------------------------#
#--------------------- Thyroidomix -------------------
#-----------------------------------------------------#

#Reading the summary results of the Thyroidomix Meta-GWAS study
thyroidomics_TSH      <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\Thyroidomix\\formatted_invnormTSH_overall_150611_invvar1.txt-QCfiltered_GC.txt.gz", header = T, sep = ",")
thyroidomics_FT4      <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\Thyroidomix\\formatted_invnormFT4_overall_150611_invvar1.txt-QCfiltered_GC.txt.gz", header = T, sep = ",")
thyroidomics_HyperTyd <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\Thyroidomix\\formatted_decTSH_overall_150602_invvar1.txt-QCfiltered_GC.txt.gz", header = T, sep = ",")
thyroidomics_HypoTyd  <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\Thyroidomix\\formatted_incTSH_overall_150602_invvar1.txt-QCfiltered_GC.txt.gz", header = T, sep = ",")

thyroidomics_HypoTyd %>% head(50) %>% View()

thyroidomics_TSH      <- thyroidomics_TSH %>% mutate_at("MarkerName", str_replace_all, ":SNP", "")
thyroidomics_FT4      <- thyroidomics_FT4 %>% mutate_at("MarkerName", str_replace_all, ":SNP", "")
thyroidomics_HyperTyd <- thyroidomics_HyperTyd %>% mutate_at("MarkerName", str_replace_all, ":SNP", "")
thyroidomics_HypoTyd  <- thyroidomics_HypoTyd  %>% mutate_at("MarkerName", str_replace_all, ":SNP", "")
#----------#

repSNPs %>%
  mutate(SNPid_37 = str_c(Chr37, ':', Pos_b37)) %>% 
  inner_join(thyroidomics_FT4,
             by = c("SNPid_37" = "MarkerName"), 
             suffix = c(".Kidney", ".Thyroid")) %>% #View()
  filter(P.value.Thyroid < 0.05) %>% #View()
  select(SNPid, Locus, RSID, 
         Allele1.log.Std, CHRIS.Allele, Effect.allel, Allele1.Thyroid,
         Allele2.log.Std, Allele2.Thyroid, 
         BETA, Effect.CHRIS, Effect.ckdgen, Effect.Kidney, Effect.Thyroid,
         SEBETA, StdErr.Kidney, StdErr.Thyroid, 
         PVALUE, P.value.Kidney, P.value.Thyroid) %>% View()
  #write.csv(., "25-Aug-2022_Significant FT4-associated SNPs found by Thyroidomics at alpha05.csv", row.names = FALSE)



#-----------------------------------------------------#
#---------------- step 3 + phenoscanner --------------
#-----------------------------------------------------#

# phenoscanner ready to join with step 3
pscanner <-
  results_step3_long_pscanner %>%
  select(Locus,  hg38_coordinates, trait, beta, se, p, pmid) %>%
  rename(SNPid = hg38_coordinates) %>%
  mutate(
    trait = str_replace(trait, "Creatinine levels$|Serum creatinine$|Creatinine$", "SCr"),
    trait = str_replace(trait, "Serum urate",                                   "Urate"),
    trait = str_replace(trait, "Magnesium levels|Serum magnesium|Magnesium",    "Magnesium"),
    trait = str_replace(trait, "Activated partial thromboplastin time",         "APTT"),
    trait = str_replace(trait, "Diastolic blood pressure",                      "DBP")
  ) %>%
  filter(
    str_detect(trait,"SCr|Urate|Magnesium|APTT|DBP"))#,
    #!if_all(c(beta, se), is.na) %>% #View
  # Joining the results of OpenTrgets look-up
  #full_join(opentargets, ., by = "SNPid", suffix = c("_OT", "_PS")) %>%
  # Joining the results of step 3
  full_join(results_Step3_long, ., by = c("SNPid" = "SNPid",
                                          "Locus" = "Locus",
                                          "pheno" = "trait")) %>%
  filter(!is.na(as.numeric(p)),
         #Locus == "SLC34A1",
         #!if_any(c(p, pval), is.na)
         ) %>% View



#-----------------------------------------------------#
#------------- GWAScatalog + phenoscanner ------------
#-----------------------------------------------------#

# To query the variants in GWAS catalog database,
# we are using the "gwasrapidd" package.

# Fetch the studies
my_studies  <- get_studies(efo_trait = 'autoimmune disease')

# You could have also used get_associations(efo_trait = 'autoimmune disease')
my_associations <- get_associations(study_id = my_studies@studies$study_id)

# Get association ids for which pvalue is less than 1e-6.
# Extract column association_id
association_ids <-
  gwasCatalog@associations %>% 
  filter(pvalue < 5e-8) %>% # Filter by p-value
  drop_na(pvalue) %>%
  pull(association_id)

# Extract associations by association id
my_associations2 <- my_associations[association_ids]

gwasrapidd::n(my_associations2)


gwasCatalog <- get_studies(variant_id = repSNPs[repSNPs$Locus == "SLC34A1", "RSID"])

# Filtering the studies by traits
gCat <- 
  gwasCatalog@studies %>% 
  filter(reported_trait %in% c("Creatinine levels", "Serum creatinine levels"))

gCat_Results <- get_associations(study_id = gCat$study_id) #%>% View   gwasCatalog@studies$study_id



# Subset the output
# Risk allele info
# Beta, SE, Pvalue
# Store association_id
gCat_id <-
  gCat_Results@genes %>%
  select(association_id, gene_name) %>%
  full_join(gCat_Results@risk_alleles %>%
              select(association_id, 
                     variant_id, 
                     risk_allele, 
                     genome_wide), by = "association_id" ) %>% 
  full_join(gCat_Results@associations %>%
              select(association_id,
                     beta_number, 
                     standard_error, 
                     pvalue),      by = "association_id") %>% 
  filter(variant_id == repSNPs[repSNPs$Locus == "SLC34A1", "RSID"]) %>% 
  pull(association_id)
  #Map an association accession identifier to an EFO trait id.

association_to_trait(gCat_id, verbose = FALSE, warnings = TRUE)
association_to_variant(gCat_id, verbose = FALSE, warnings = TRUE)

#----------#


results_step3_long_pscanner %>% colnames()

opentargets %>% 
  filter(str_detect(traitReported, "serum"))
#----------#


library(httr)
library(jsonlite)

# Define your list of variants of interest
variants <- c("rs3812036", "rs41274482")
#variants <- c("rs123456", "rs789012", "rs345678")

# Define the base URL for the GWAS catalog API
base_url <- "https://www.ebi.ac.uk/gwas/summary-statistics/api"

# Loop through each variant and query the GWAS catalog API
for (i in 1:length(variants)) {
  
  # Define the API endpoint for the variant
  endpoint <- paste0(base_url, "/single-variant/", variants[i])
  
  # Query the API and parse the JSON response
  response <- GET(endpoint)
  if (status_code(response) == 200) {
    result <- fromJSON(content(response, "text"), simplifyDataFrame = TRUE)
    if (result$total > 0) {
      print(paste0("Variant: ", variants[i]))
      print(paste0("Traits: ", paste(result$data$trait, collapse = ", ")))
    } else {
      print(paste0("No associated traits found for variant ", variants[i]))
    }
  } else {
    print(paste0("Error: Unable to retrieve data for variant ", variants[i]))
  }
  
  # Pause for a few seconds to avoid overloading the API
  Sys.sleep(3)
}
#----------#

#----------#

library(gwasrapidd)

# Loop through each variant and query the GWAS catalog
for (i in 1:length(variants)) {
  results <- gwasrapidd::get_associations(variant_id = variants[i])
  
  # Check if any associations were found
  if (nrow(results$data) > 0) {
    print(paste0("Variant: ", variants[i]))
    print(paste0("Traits: ", paste(results$data$trait, collapse = ", ")))
  } else {
    print(paste0("No associated traits found for variant ", variants[i]))
  }
  
  # Pause for a few seconds to avoid overloading the API
  Sys.sleep(3)
}
#----------#

library(httr)

# Loop through each variant and query the GWAS catalog
for (i in 1:length(variants)) {
  url <- paste0("https://www.ebi.ac.uk/gwas/summary-statistics-api/variant/", variants[i], "/associations")
  response <- httr::GET(url)
  
  # Check if the request was successful
  if (httr::status_code(response) == 200) {
    data <- httr::content(response, "text", encoding = "UTF-8")
    data <- jsonlite::fromJSON(data)$data
    
    # Check if any associations were found
    if (length(data) > 0) {
      print(paste0("Variant: ", variants[i]))
      print(paste0("Traits: ", paste(data$trait, collapse = ", ")))
    } else {
      print(paste0("No associated traits found for variant ", variants[i]))
    }
  } else {
    print(paste0("Error: Unable to retrieve data for variant ", variants[i]))
  }
  
  # Pause for a few seconds to avoid overloading the API
  Sys.sleep(3)
  
  return(data)
}
#----------#




#-----------------------------------------------------#
#------------- GWAScatalog + phenoscanner ------------
#-----------------------------------------------------#

gwas_cat <- read.csv("D:\\Dariush\\PhD\\Analysis\\Outputs_tables\\GWAScatalog\\18-Feb-23_summary_results_all.csv", sep = "\t")

gwas_cat %>%
  filter(Pvalue <= 5e-8) %>%
  rename(trait = Trait,
         beta = Beta,
         se = Stderr) %>% 
  mutate(
    trait = str_replace(trait, "Creatinine levels$|Serum creatinine$|Creatinine$", "SCr"),
    trait = str_replace(trait, "Serum urate",                                   "Urate"),
    trait = str_replace(trait, "Magnesium levels|Serum magnesium|Magnesium",    "Magnesium"),
    trait = str_replace(trait, "Activated partial thromboplastin time",         "APTT"),
    trait = str_replace(trait, "Diastolic blood pressure",                      "DBP")
  ) %>%
  filter(str_detect(trait, "SCr|Urate|Magnesium|APTT|DBP")) %>%
  right_join(repSNPs %>% select(SNPid, Locus, RSID),
             ., by = c("RSID" = "variant"),
            suffix = c("_CHRIS", "_GC")) %>%
  rename(Locus = Locus_CHRIS) %>%
  right_join(pscanner, 
             by = c("SNPid", "Locus", "trait"),
             suffix = c("_GC", "_PS")) %>%
  mutate(across(c(beta_PS, se_PS), as.numeric)) %>%
  # Remove duplicates on selected columns
  #distinct(SNPid, beta_PS, se_PS, .keep_all = T) %>%
  mutate(missing = case_when(
    Locus == "IGF1R"  & if_all(c(beta_PS, se_PS), is.na) ~ 1,
    Locus == "SHROOM3" & if_all(c(beta_PS, se_PS), is.na) & pmid == "20700443" ~ 1,
                             TRUE ~ 0)) %>%
  filter(missing != 1) %>%
  select(-missing) %>%
  write.csv(., "18-Feb-23_Summary of interrogation in GWAS catalog merged with Phenoscanner.csv", row.names = FALSE)
  
  # full_join(results_Step3_long, ., by = c("SNPid" = "SNPid",
  #                                         "Locus" = "Locus",
  #                                         "pheno" = "trait"))
  

 
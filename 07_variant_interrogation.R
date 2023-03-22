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
           show.legend = FALSE,
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
#----------#

# Preparing Phenoscanner summary results for joining to step 3
pscanner <-
  results_step3_long_pscanner %>% 
  #filter(pmid == "20383146") %>%  
  #count(trait)
  #write.csv("23-Feb-23_pubmed_id 20383146.csv", row.names = FALSE, quote=FALSE)
  select(Locus,  hg38_coordinates, rsid, EA, trait, beta, se, p, ancestry, pmid, efo) %>%
  rename(SNPid = hg38_coordinates) %>%
  mutate(
    across(c(beta, se), as.numeric),
    trait = str_replace(trait, "Creatinine levels$|Serum creatinine$|Creatinine$", "SCr"),
    trait = str_replace(trait, "Serum urate",                                      "Urate"),
    trait = str_replace(trait, "Magnesium levels|Serum magnesium|Magnesium",       "Magnesium"),
    trait = str_replace(trait, "Activated partial thromboplastin time",            "APTT"),
    trait = str_replace(trait, "Diastolic blood pressure",                         "DBP"),
    beta  = case_when(rsid == "rs13146355" & efo == "NCIT_C61035" & pmid == "20700443" ~ 0.005,   TRUE ~ beta),
    se    = case_when(rsid == "rs13146355" & efo == "NCIT_C61035" & pmid == "20700443" ~ 0.001,   TRUE ~ se),
    beta  = case_when(rsid == "rs13146355" & pmid == "22797727" & efo  == "EFO_0004518" ~ 0.0047, TRUE ~ beta),
    se    = case_when(rsid == "rs13146355" & pmid == "22797727" & efo  == "EFO_0004518" ~ 0.0007, TRUE ~ se),
    beta  = case_when(Locus == "SLC34A1" & rsid == "rs3812036"  & pmid == "22703881" ~ -0.419985, TRUE ~ beta),
    se    = case_when(Locus == "SLC34A1" & rsid == "rs3812036"  & pmid == "22703881" ~ 0.0481015, TRUE ~ se),
    missing = case_when(Locus == "IGF1R"    & if_all(c(beta, se), is.na) ~ 1,
                        pmid  == "19430482" & if_all(c(beta, se), is.na) ~ 1,
                        pmid  == "20383146" & if_all(c(beta, se), is.na) ~ 1,
                        pmid  == "20700443" & p == 6.000e-13 ~ 1,
                        pmid  == "20700443" & p == 6.270e-13 ~ 1,
                        pmid  == "UKBB" ~ 1, # duplicate associations
                        TRUE ~ 0),
    dataset = "PS"
    ) %>%
  filter(str_detect(trait,"SCr|Urate|Magnesium|APTT|DBP"),
         missing != 1) %>%
  select(-missing, -efo)



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
#------------- GWAScatalog (gwasrapidd) --------------
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

#-----------------------------------------------------#
#---------------- GWAScatalog Rest API ---------------
#-----------------------------------------------------#

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

# Summary results of SNPs look-up via GWAS Catalog REST API in python
gwas_cat0 <- read.csv("D:\\Dariush\\PhD\\Analysis\\Outputs_tables\\GWAScatalog\\18-Feb-23_Summary results of interrogation in GWAS catalog.csv", sep = "\t")

# with ancestry info
gwas_cat <- read.csv("D:\\Dariush\\PhD\\Analysis\\Outputs_tables\\GWAScatalog\\24-Feb-23_GWAS Catalog interrogation summary results.csv")
#----------#

# Pooling GC interrogation with PS variants look-up
gc_pscanner <-
  gwas_cat %>%
  rename(beta  = Beta,
         se    = Stderr,
         p     = Pvalue,
         trait = Trait,
         EA    = Risk_A,
         rsid  = variant,
         ancestry = Ancestry,
         pmid  = Pubmed_ID) %>%
  filter(rsid != "rs5020545", #For this variant, the result shows the mGWAS association
         p <= 5e-8) %>%
  mutate(
    #across(c(beta, se), as.numeric),
    trait = str_replace(trait, "Serum creatinine levels$|Creatinine levels$|Serum creatinine$|Creatinine$", "SCr"),
    trait = str_replace(trait, "Urate levels$|Urea levels|Serum urate|Serum uric acid levels", "Urate"),
    trait = str_replace(trait, "Magnesium levels|Serum magnesium|Magnesium", "Magnesium"),
    trait = str_replace(trait, "Activated partial thromboplastin time", "APTT"),
    trait = str_replace(trait, "Diastolic blood pressure", "DBP"),
    trait = str_replace(trait, "Systolic blood pressure",  "SBP"),
    # For pmid 20700443, the original EA was G. For consistency, we didn't
    # change incorrect EA = A, but negate the effect size of rs13146355
    beta  = case_when(rsid == "rs13146355" & pmid == "20700443" & trait == "Magnesium" ~ 0.005, TRUE ~ beta),
    se    = case_when(rsid == "rs13146355" & pmid == "20700443" & trait == "Magnesium" ~ 0.001, TRUE ~ se),
    p     = case_when(rsid == "rs13146355" & pmid == "20700443" & trait == "Magnesium" ~ 6.270E-13, TRUE ~ p),
    p     = case_when(rsid == "rs13146355" & pmid == "34899825" & trait == "Urate"     ~ 9.113E-142,TRUE ~ p),
    p     = case_when(rsid == "rs28817415" & pmid == "36329257" & trait == "SCr"       ~ 1.161E-23, TRUE ~ p),
    p     = case_when(rsid == "rs7042786"  & pmid == "36329257" & trait == "SCr"       ~ 2.417E-14, TRUE ~ p),
    p     = case_when(rsid == "rs55940751" & pmid == "27841878" & trait == "SBP"       ~ 3.600E-8,  TRUE ~ p),
    p     = case_when(rsid == "rs56019566" & pmid == "35213538" & trait == "SCr"       ~ 3.7E-09,   TRUE ~ p),
    p     = case_when(rsid == "rs77924615" & pmid == "31578528" & trait == "Urate"     ~ 1.27E-11,  TRUE ~ p),
    p     = case_when(rsid == "rs28394165" & pmid == "34594039" & trait == "SCr"       ~ 2.26E-189, TRUE ~ p),
    p     = case_when(rsid == "rs77924615" & pmid == "34594039" & trait == "SCr"       ~ 4.36E-204, TRUE ~ p), 
    # Aligning the effect allele and change the effect direction
    EA    = case_when(rsid == "rs13146355" & pmid == "20700443" ~ "A", TRUE ~ EA),
    EA    = case_when(rsid == "rs77924615" & pmid == "34594039" ~ "A", TRUE ~ EA),
    EA    = case_when(rsid == "rs56019566" & pmid == "35213538" ~ "T", TRUE ~ EA),
    EA    = case_when(rsid == "rs10008637" & pmid == "35213538" ~ "C", TRUE ~ EA),
    EA    = case_when(rsid == "rs4744712"  & pmid == "33356394" ~ "C", TRUE ~ EA),
    EA    = case_when(rsid == "rs28394165" & pmid == "34594039" ~ "C", TRUE ~ EA),
    beta  = case_when(rsid == "rs77924615" & pmid == "34594039" ~ -beta, TRUE ~ beta),
    beta  = case_when(rsid == "rs56019566" & pmid == "35213538" ~ -0.0218485, TRUE ~ beta),
    beta  = case_when(rsid == "rs10008637" & pmid == "35213538" ~ 0.0520574,  TRUE ~ beta),
    beta  = case_when(rsid == "rs77924615" & pmid == "31578528" ~ -0.0273 ,   TRUE ~ beta),
    beta  = case_when(rsid == "rs807624"   & pmid == "33356394" ~ -beta ,   TRUE ~ beta),
    beta  = case_when(rsid == "rs819196"   & pmid == "31015462" & trait == "SCr" ~ -0.04046, TRUE ~ beta),
    beta  = case_when(rsid == "rs4237268"  & pmid == "31015462" & trait == "SCr" ~ -0.04639, TRUE ~ beta),
    dataset = "GC") %>% 
  filter(str_detect(trait, "SCr|Urate|Magnesium|APTT|DBP|SBP")) %>%
  select(-Locus) %>%
  left_join(.,
            repSNPs %>% select(Locus, RSID, SNPid),
            by = c("rsid" = "RSID")) %>%
  # Remove duplicates association for each trait by 
  # taking the one with smallest p-value
  #distinct(SNPid, beta_PS, se_PS, .keep_all = T) %>%
  add_count(SNPid, trait) %>% 
  mutate(missing = case_when(n > 1 & if_any(c(beta, se), is.na) ~ 1,
                             # Remove duplicated variant in SHROOM3 
                             # available also in pscanner -> not now: removed it later
                             #Locus == "SHROOM3"   & p == 6.27E-13 ~ 1,
                             rsid  == "rs3812036" ~ 1,
                             TRUE ~ 0)
         ) %>%
  #count(SNPid, trait)
  filter(missing != 1) %>%
  group_by(SNPid, trait) %>%
  slice_min(p, n = 1) %>%
  #filter(Pvalue == min(Pvalue, na.rm = TRUE)) %>%
  ungroup() %>%
  # Checking if exclusion criteria correctly applied to associations -> 24 remained
  #select(SNPid, n, Pvalue, trait) %>% count(SNPid, trait) 
  select(Locus, SNPid, rsid, EA, trait, beta, se, p, ancestry, pmid, dataset) %>%
  # merging with phenoscanner interrogation results
  rbind(pscanner) %>%
  filter(pmid != "22797727") %>%
  rename(Estimate_step3b = beta,
         SE_step3b     = se,
         Pvalue_step3b = p,
         EA_step3b     = EA) #%>%
  #write.csv(., "26-Feb-23_GWAS Catalog interrogation summary results merged with Phenoscanner.csv", row.names = FALSE)
#----------#

# PheWeb-like scatter plot
gc_pscanner %>%
  ggplot(aes(Locus, y = -log10(p), color = trait, shape = trait)) +
  geom_jitter(size = 2.5) + #show.legend = FALSE
  ggrepel::geom_text_repel(aes(label = trait), size = 5, 
                           color = "grey10", max.overlaps = 5) +
  scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  labs(x = "") +
  #guides(color = guide_legend(), shape = "none") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
        panel.grid.major.x = element_blank(),
        axis.text  = element_text(size = 8,  face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position =  c(.92, .7), #"none",
        legend.key.size  = unit(0.99, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))

ggsave("26-Feb-23_GWAS Catalog merged with pscanner 1.png", 
       last_plot(), width = 9, height = 5.5, pointsize = 5, dpi = 300, units = "in")


  
#-----------------------------------------------------#
#---------------- step 3 + phenoscanner --------------
#-----------------------------------------------------#
  
# GWAS Catalog + Phenoscanner + OpenTargets
gc_pscanner %>%
  # Joining the results of OpenTrgets look-up
  full_join(opentargets, ., by = "SNPid", suffix = c("_OT", "_PS")) %>%
  # Joining the results of step 3
  full_join(results_Step3_long, ., by = c("SNPid" = "SNPid",
                                          "Locus" = "Locus",
                                          "pheno" = "trait")) %>%
  filter(!is.na(as.numeric(p))) %>% View
#----------#

# Joining the results of interrogation in Phenoscanner and 
# GWAS Catalog with mediation analysis steps

# Supplementary table 4: Mediation analysis results

sum3steps_long <-
  repSNPs %>%
  # Step1 or GWAS
  select(SNPid, Locus, RSID, EA_CHRIS_disc, RA_CHRIS_disc,                   
         Beta_CHRIS, SE_CHRIS, Pvalue_CHRIS) %>%
  # Joining with step2
  inner_join(results_Step2_long, by = c("SNPid", "Locus")) %>%
  # Joining with step3
  inner_join(results_Step3_long, by = c("SNPid" = "SNPid",
                                        "Locus" = "Locus",
                                        "Trait" = "pheno"),
             suffix = c("_step2", "_step3")) %>%
  # Joining with interrogation results
  full_join(gc_pscanner, by = c("SNPid" = "SNPid",
                                "Locus" = "Locus",
                                "Trait" = "trait",
                                "RSID"  = "rsid")) %>%
  # Joining with OpenTrgets look-up
  #full_join(opentargets, ., by = "SNPid", suffix = c("_OT", "_PS")) %>%
  # Taking the common associated variants in step3 and interrogation 
  #filter(!is.na(Pvalue_step3b)) %>%
  mutate(Pvalue_CHRIS_1s = Pvalue_CHRIS/2,
         # Check concordance of effect allele in CHRIS with interrogated study
         EA_step3ab              = ifelse(EA_CHRIS_disc == EA_step3b, "concordant", "discor"),
         Estimate_step3b_aligned = ifelse(EA_CHRIS_disc == EA_step3b, Estimate_step3b, -Estimate_step3b),
         Assoc_dir_step3ab       = ifelse(Estimate_step3 * Estimate_step3b_aligned < 0,
                                      "inconsis", "consistent"),
         # compute the missing standard error based on beta ann p-value of the association
         SE_estimated = abs(Estimate_step3b/ qnorm(Pvalue_step3b/2)), # tail = 2
         # replace the missing standard error with the computed one: only for two SNPs
         SE_step3b    = ifelse(is.na(SE_step3b), SE_estimated, SE_step3b),
         # Check if the replacement is done correctly
         SE_check     = ifelse(SE_step3b == SE_estimated, "Yes", " ")) %>%
  # drop the interrogated association results for the variants without outlier effect in step 2
  mutate_at(vars(SE_step3b, Pvalue_step3b, Estimate_step3b_aligned, pmid), ~replace(., outlier == "No", NA)) %>%
  # rename related column
  rename(associated_step3a = related) %>%
  # rename the association effect column and replace it with the corrected one
  mutate(Estimate_step3b = Estimate_step3b_aligned,
         # correct the binary indicator for significant association in step 3a or step 3b
         associated_step3b = ifelse(Pvalue_step3b < 5e-08, "Yes", NA)) %>%
  select(-c(EA_step3b, EA_step3ab, ancestry, dataset, Estimate_step3b_aligned, Assoc_dir_step3ab, SE_estimated, SE_check)) %>%
  #write.csv(., "06-Mar-23_Step 3 mediation analysis merged with interrogation results.csv", row.names = FALSE)
  #----------#
  # ggplot(aes(x = Estimate_step3b_aligned, y = Estimate_step3, color = Locus, shape = Locus)) +
  # geom_abline(slope = 1, lty = 2) +
  # geom_vline(xintercept = 0) +
  # geom_hline(yintercept = 0) +
  # geom_point(alpha = 0.9, size = 3.5) +
  # scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  # scale_color_manual(values = c(
  #            "maroon1",
  #            "darkorchid2",
  #            "orange2", 
  #            "green4", 
  #            "steelblue2",
  #            "darkturquoise",
  #            "tomato",
  #            "springgreen2",
  #            "royalblue2",
  #            "gold",
  #            "grey50")) +
  # labs(x = "Effect size of variants association with traits found in databases",
  #      y = "Effect size of variants association with traits found in CHRIS") +
  # theme(panel.background = element_rect(fill = "white"),
  #       panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
  #       panel.grid.major.x = element_line(linetype = 'solid', color = "grey80", size = .15),
  #       axis.text  = element_text(size = 8,  face = "bold"),
  #       axis.title = element_text(size = 12, face = "bold"),
  #       #legend.position = "none",
  #       legend.key.size  = unit(0.99, 'cm'),
  #       legend.key.width = unit(0.7, 'cm'),
  #       legend.text  = element_text(size = 12),
  #       legend.title = element_text(size = 14, face = "bold"))
  # ggsave("07-Mar-23_Association effect in Step 3 vs in interrogation only with outlier in step 2.png", 
  #        last_plot(), width = 9, height = 5.5, pointsize = 5, dpi = 300, units = "in")
  #----------#
  # Joining with step4
  inner_join(results_Step4_long %>% 
               rename(Estimate_step4 = Estimate, SE_step4 = SE, Pvalue_step4 = Pvalue),
             by = c("SNPid" = "SNPid",
                    "Locus" = "Locus",
                    "Trait" = "pheno")) %>%
  mutate(Mediator = case_when(outlier == "Yes" & associated_step3a == "Yes" ~ "Yes",
                              outlier == "Yes" & associated_step3b == "Yes" ~ "Yes",
                              TRUE ~ "No"),
         Change_P = round((Beta_CHRIS - Estimate_step2) / Beta_CHRIS * -100, 2),
         EA_OA = paste0(EA_CHRIS_disc, "/", RA_CHRIS_disc)) %>% 
  rename(Estimate_GWAS  = Beta_CHRIS,
         SE_GWAS        = SE_CHRIS,
         Pvalue_GWAS    = Pvalue_CHRIS) %>%
  select(Locus, RSID, SNPid, EA_OA, ends_with("GWAS"), Pvalue_CHRIS_1s, everything()) %>% 
  select(-c(EA_CHRIS_disc, RA_CHRIS_disc)) #%>%
  #write.csv(., "08-Mar-23_SNP-wise summary of mediation analysis steps.csv", row.names = FALSE, quote = FALSE)
#----------#

# Table 3 of the paper - mediator variants

sum3steps_long %>%
  filter(Mediator == "Yes",
         Trait    != "SCr") %>%
  as_tibble() #%>%
  #write.csv("07-Mar-23_Table 3_Step1&2&3&4_17 mediator SNPs either in CHRIS or in consortia.csv", row.names = F, quote = F)


# Finally the first part of the project finished on Tuesday at 16:30, 07-Mar-2023. 
# The last modification was to add the 1-sided p-values of the GWAS to the mediation analysis results on 16:40, 08-Mar-2023. 







 
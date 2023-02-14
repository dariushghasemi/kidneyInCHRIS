#=============================================#
# Interrogating replicated SNPs in databases
#=============================================#

#-----------------------------------------------------#
#---------------- Mendelian Randomization ------------
#-----------------------------------------------------#
library(MendelianRandomization)

#queries the PhenoScanner database of genotype-phenotype associations
phenoscanner(snpquery="rs807624")$results %>% View()

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
library(MendelianRandomization)
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
  write.csv("25-Aug-2022_GenomeWideSig SNPs associated with Coagulation and Magnesium queried using phenoscanner.csv", row.names = FALSE, quote=FALSE)

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
  write.csv("31-Jan-23_163 SNPs queried using phenoscanner at genome-wide sig associated with trait_Frequency.csv", row.names = FALSE, quote=FALSE)
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

ggsave("23-Jun-22_Frequency of the traits found in Phenoscanner at GWS level.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#-----------------------------------------------------#
#---------------- step 3 + phenoscanner --------------
#-----------------------------------------------------#

results_step3_long_pscanner %>% View
  select(Locus,  hg38_coordinates, trait, beta, se, p, pmid) %>%
  rename(SNPid = hg38_coordinates) %>%
  #filter(trait == "Creatinine levels" | trait == "Creatinine" | trait == "Serum creatinine")
  #filter(trait == "Serum urate")
  mutate(
    trait = str_replace(trait, "Creatinine levels$|Serum creatinine$|Creatinine$", "SCr"),
    trait = str_replace(trait, "Serum urate",                                   "Urate"),
    trait = str_replace(trait, "Magnesium levels|Serum magnesium|Magnesium",    "Magnesium"),
    trait = str_replace(trait, "Activated partial thromboplastin time",         "APTT"),
    trait = str_replace(trait, "Diastolic blood pressure",                      "DBP")
    ) %>%
  filter(
    str_detect(trait,"SCr|Urate|Magnesium|APTT|DBP"),
    !if_all(c(beta, se), is.na)
         ) %>%
         #!is.na(as.numeric(p)))
  full_join(results_Step3_long, ., by = c("SNPid" = "SNPid",
                                          "Locus" = "Locus",
                                          "pheno" = "trait")) %>%
  filter(Locus == "SLC34A1") %>% 
  View




   #-----------------------------------------------------#
#--------------------- OpenTargets -------------------
#-----------------------------------------------------#

#Reading the traits found to be related 
#with the 162 SNPs in OpenTargets by David
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

ggsave("13-Aug-22_Frequency of the traits found in OpenTargets at GWS level_grouped_by_Locus4.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")


#-----------------------------------------------------#
#--------------------- Thyroidomix --------------------
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

 
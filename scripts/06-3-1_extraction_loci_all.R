#!/usr/bin/Rscript

#=========================================#
#  Grepping all 147 CKDGen Loci in CHRIS
#=========================================#

# Initiated on December 22, 2023.
# Last modifications on 07-Dec-23

#---------------#
library(tidyverse)

# take the date and time
format(Sys.time(), "%a %b %e %H:%M:%S %Y")

#---------------#
# input files
base.dir = "/home/dghasemisemeskandeh/projects/gwas"
GWAS = paste0(base.dir, "/Output/ReformedScheme/eGFRw.log.Res.txt.gz")
SNPs = paste0(base.dir, "/04_pairways_LD/output/test_convert_prepared.txt")
BED38 = paste0(base.dir, "/02_liftover_ckdgen/output/MetaBED38.txt")
Meta_EA = paste0(base.dir, "/01_subtract_ckdgen/output/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz")
Meta_TA = paste0(base.dir, "/01_subtract_ckdgen/output/meta.results_corrected.with.MetaSubtract.txt.gz")
ckdgen_loci = paste0(base.dir, "/Data/Supptable.txt")
snpSigGFR = paste0(base.dir, "/06_replication_analysis/output/replicatedSNPs_eGFR.Std.txt")
out.dir = paste0(base.dir, "/06_replication_analysis/output")
OUT1 = paste0(out.dir, "/31-Dec-22_Replication_of_147_EA_Loci_in_CHRIS_NOTsorted.txt")
OUT2 = paste0(out.dir, "/07-Dec-23_Replication_of_147_Trans-A_loci_in_CHRIS.txt")

no.rows=100000
no.rows=-1

#--------------------------------------------------------------#
#                      Read input files                        #
#--------------------------------------------------------------#

names.gwas <- c("CHR","POS38","END","MARKER_ID_CHRIS","NS_CHRIS", "AC_CHRIS","CALLRATE_CHRIS","GENOCNT_CHRIS","MAF_CHRIS","STAT_CHRIS","Pvalue_CHRIS","Beta_CHRIS","SE_CHRIS","R2_CHRIS")
names.meta <- c("Chr37","Pos_b37","RSID","EA_CKDGen","RA_CKDGen", "EAF_CKDGen","Beta_CKDGen","SE_CKDGen","Pvalue_CKDGen","n_total_sum_CKDGen")

#---------------#
cat("\n\nImporting the GWAS results ...\n\n")
gwas    <- read.table(GWAS, header = F, col.names = names.gwas, sep = '\t', stringsAsFactors = F) #, nrows=no.rows

cat("\n\nImporting the MetaSubtract results ...\n\n")
metaB37 <- read.table(Meta_TA, header = F, col.names = names.meta, sep = '\t', stringsAsFactors = F) #, nrows=no.rows

cat("\n\nImporting MetaSubtract lifted positions to GRCh38 ...\n\n")
POS 	<- read.table(BED38, header = T, sep = '\t', stringsAsFactors = F, fill = T) #, nrows=no.rows

cat("\n\nImporting the 147 CKDGen Loci for Trans-A or EA ...\n\n")
loci_147 <- read.delim(ckdgen_loci, header = T, sep= '\t', row.names = NULL, stringsAsFactors = F)

#---------------#
# drop unnecessary columns
gwas$END <- NULL

head(gwas, 5)
head(metaB37)
head(POS)

#--------------------------------------------------------------#
#                     Replication Analysis                     #
#--------------------------------------------------------------#

cat("\n\nAppending the GRCh38 cooardiantes to MetaSubtract results...\n\n")
cat("\n\nExtracting 147 loci from subtracted CKDGen Meta-GWAS...\n\n")

#---------------#
# 147 CKDGen Loci for Trans-A and later for EA
meta_loci <- metaB37 %>%
   # Adding POS38 to subtracted Meta-GWAS
   merge(., POS, by.x = c("Chr37", "Pos_b37"), all = FALSE) %>%
   # drop unnecessary columns
   select(- Chr37) %>%
   # Making column names consistant for merge with GWAS in CHRIS
   rename(CHR = Chr38, POS37 = Pos_b37, POS38 = Pos_b38) %>%
   # Adding locus name to meta-gwas
   inner_join(loci_147 %>% select("RSID", "Closest.Gene"), by = "RSID") %>%
   # change class from character to integer for merge
   mutate(CHR = as.integer(CHR)) %>%
   # Extracting 147 CKDGen Loci from CHRIS GWAS
   left_join(gwas, by = c("CHR", "POS38")) %>%
   # indicate if a locus is replicated in CHRIS
   mutate(status_CHRIS = if_else(Pvalue_CHRIS  < 6.8e-4, "Sig", "Non-Sig")) %>%
   # relocate columns
   select(CHR, POS37, POS38, RSID, Closest.Gene, everything()) %>%
   # sort results based on chrom and positions
   arrange(CHR, POS38)


#--------------------------------------------------------------#
#                    Save analysis output                      #
#--------------------------------------------------------------#

# save results of the replication analysis in CHRIS 
# on all CKDGen 147 loci Multi-A or later EA
write.table(meta_loci, file = OUT2, sep = "\t", row.names = FALSE, quote = FALSE)

head(meta_loci)

#---------------#
# Number of replicated (significant) loci in CHRIS
cat("\n\nNumber of replicated loci in CHRIS out of 147 Trnas-A or EA CKDGen Loci is ...\n\n")
table(meta_loci$status_CHRIS)

# take the date and time
format(Sys.time(), "%a %b %e %H:%M:%S %Y")

#---------------#
#sbatch --wrap 'Rscript 06-2_extraction_loci_all.R' --mem=8GB -J "06-2.R"

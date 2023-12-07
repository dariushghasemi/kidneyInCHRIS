#!/usr/bin/Rscript

#=========================================#
#   Replication of CKDGen Loci in CHRIS
#=========================================#

#---------------#
library(tidyverse)

# take the date and time
format(Sys.time(), "%a %b %e %H:%M:%S %Y")

#---------------#
# input files
base.dir = "/home/dghasemisemeskandeh/projects/gwas"
GWAS = paste0(base.dir, "/Output/ReformedScheme/eGFRw.log.Res.txt.gz")
SNPs = paste0(base.dir, "/04_pairways_LD/output/test_convert_prepared.txt")
Meta = paste0(base.dir, "/01_subtract_ckdgen/output/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz")
BED38 = paste0(base.dir, "/02_liftover_ckdgen/output/MetaBED38.txt")
snpSigGFR = paste0(base.dir, "/06_replication_analysis/output/replicatedSNPs_eGFR.Std.txt")
out.dir = paste0(base.dir, "/06_replication_analysis/output")
OUT1 = paste0(out.dir, "/10-Jan-23_Meta_replicatedSNPs_eGFRw.log.Res.all.txt")
OUT2 = paste0(out.dir, "/10-Jan-23_replicatedSNPs_eGFRw.log.Res.txt")

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
metaB37 <- read.table(Meta, header = F, col.names = names.meta, sep = '\t', stringsAsFactors = F) #, nrows=no.rows

cat("\n\nImporting MetaSubtract lifted positions to GRCh38 ...\n\n")
POS 	<- read.table(BED38, header = T, sep = '\t', stringsAsFactors = F, fill = T) #, nrows=no.rows

cat("\n\nImporting the SNPs with strong LD ...\n\n")
snps    <- read.delim(SNPs,  header = FALSE, col.names = c("CHR", "POS38"), sep= '\t', row.names = NULL, stringsAsFactors = F)

cat("\n\nImporting the 147 CKDGen Loci for Trans-A or EA ...\n\n")
loci_147 <- read.delim(ckdgen_loci, header = T, sep= '\t', row.names = NULL, stringsAsFactors = F)

#---------------#
# Making column names consistant
meta$Chr37 <- NULL
names(meta) [names(meta) == "Chr38"]   <- "CHR"
names(meta) [names(meta) == "Pos_b37"] <- "POS37"
names(meta) [names(meta) == "Pos_b38"] <- "POS38"
head(meta)
   
#--------------------------------------------------------------#
#                      Taking Proxy SNPs                       #
#               Extracting SNPs from GWAS results              #
#--------------------------------------------------------------#

# N. rows of proxy snps
cat("\n\n N. proxy snps is:", nrow(snps), "\n\n")

# Find the row name of the last proxy SNP RPL3L, the locus before PDILT locus
RPL3L <- which(snps$POS38 == 1972803)

cat("\n\nExtracting SNPs having strong LD from GWAS results...\n\n")
cat("\n\nThis is the merged file consisting of significants SNPs with strong LD...\n\n")

# Merge GWAS and 6337 proxy variants with CKDGen hits
snpGFRlog <- snps %>%
   # Add missing SNP in PDILT locus after RPL3L locus
   add_row(CHR = 16, POS38 = 20381010, .after = RPL3L) %>%
   # Taking in-strong LD variants from GWAS in CHRIS
   inner_join(gwas, by = c("CHR", "POS38")) %>%
   # Filtering for significant variants
   filter(Pvalue_CHRIS < 6.8e-4)


head(snpGFRlog, 5)

#--------------------------------------------------------------#
#                     Replication in CHRIS                     #
#--------------------------------------------------------------#

meta <- metaB37 %>%
   # Adding POS38 to subtracted Meta-GWAS
   merge(., POS, by.x = c("Chr37", "Pos_b37"), all = FALSE) %>%
   select(- Chr37) %>%
   # Making column names consistant for merge with GWAS in CHRIS
   rename(CHR = Chr38, POS37 = Pos_b37, POS38 = Pos_b38) %>%
   # taking all GWAS signif. SNPs in CHRIS from subtracted GWAS in CKDGen EA
   right_join(snpGFRlog, by = c("CHR", "POS38")) %>%
   # indicate significant SNPs
   mutate(status_CHRIS  = if_else(Pvalue_CHRIS  < 6.8e-4, "Sig", "Non-Sig"),
          status_CKDGen = if_else(Pvalue_CKDGen < 5e-08,  "Sig", "Non-Sig")) %>%
   # Sort results based on CHR, POS
   arrange(CHR, POS38)


#--------------------------------------------------------------#
#                       Saving output                          #
#--------------------------------------------------------------#

# The old name of the last of the above files was:
# `Meta_replicatedSNPs_eGFRw.log.Res.all.txt`

# saving the CKDGen and CHRIS merged data, teh replication result 
write.table(result, file = OUT1, sep = "\t", row.names = FALSE, quote = FALSE)

# save the sorted result
write.table(result, file = OUT2, sep = "\t", row.names = FALSE, quote = FALSE)

#---------------#
# sort in the server if needed
#cat 10-Jan-23_Meta_replicatedSNPs_eGFRw.log.Res.all.txt | sort -g -k 1 > 10-Jan-23_replicatedSNPs_eGFRw.log.Res.txt

# take the date and time
format(Sys.time(), "%a %b %e %H:%M:%S %Y")

#---------------#
# run on the server
#sbatch --wrap 'Rscript extraction.R' --mem=65536 -J "repInCHRIS"
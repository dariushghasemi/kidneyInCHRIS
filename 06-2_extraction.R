#!/usr/bin/Rscript


#---------------#
# input files
base.dir = "/home/dghasemisemeskandeh/projects/gwas"
GWAS = paste0(base.dir, "/Output/ReformedScheme/eGFRw.log.Res.txt.gz")
SNPs = paste0(base.dir, "/04_pairways_LD/output/test_convert_prepared.txt")
Meta = paste0(base.dir, "/01_subtract_ckdgen/output/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz")
BED38 = paste0(base.dir, "/02_liftover_ckdgen/outputMetaBED38.txt")
snpSigGFR = paste0(base.dir, "/06_replication_analysis/output/replicatedSNPs_eGFR.Std.txt")
OUT1 = paste0(base.dir, "/06_replication_analysis/output/10-Jan-23_Meta_replicatedSNPs_eGFRw.log.Res.all.txt")
OUT2 = paste0(base.dir, "/06_replication_analysis/output/10-Jan-23_replicatedSNPs_eGFRw.log.Res.txt")

max.no.rows=100000

#---------------#
library(tidyverse)

#--------------------------------------------------------------#
#                        Step 1: GWAS                          #
#--------------------------------------------------------------#
cat("\n\nImporting the GWAS results ...\n\n")

names   <- c("CHR","POS38","END","MARKER_ID_CHRIS","NS_CHRIS", "AC_CHRIS","CALLRATE_CHRIS","GENOCNT_CHRIS","MAF_CHRIS","STAT_CHRIS","Pvalue_CHRIS","Beta_CHRIS","SE_CHRIS","R2_CHRIS")
gwas    <- read.table(GWAS, header = FALSE, col.names = names, sep = '\t', stringsAsFactors = F)#, nrows=max.no.rows

# drop unnecessary columns
gwas$END <- NULL
head(gwas, 5)

#--------------------------------------------------------------#
#                Step 2,3: Merging Meta and BED38              #
#--------------------------------------------------------------#

cat("\n\nImporting the MetaSubtract results to add the GRCh38 lifted positions ...\n\n")
nameMeta <- c("Chr37", "Pos_b37", "RSID", "EA_CKDGen", "RA_CKDGen", "EAF_CKDGen", "Beta_CKDGen", "SE_CKDGen", "Pvalue_CKDGen", "n_total_sum_CKDGen")
metaB37  <- read.table(Meta,  header = FALSE, col.names = nameMeta, sep = '\t', stringsAsFactors = F)#, nrows=max.no.rows

head(metaB37)

#---------------#
cat("\n\nImporting CKDGen GRCh38 lifted positions ...\n\n")
POS 	<- read.table(BED38, header = TRUE, sep = '\t', stringsAsFactors = F, fill = TRUE)
head(POS)

#---------------#
cat("\n\nAppending the GRCh38 cooardiantes to MetaSubtract results ...\n\n")
meta <- merge(metaB37, POS, by.x = c("Chr37", "Pos_b37"), all = FALSE)

# Making column names consistant
meta$Chr37 <- NULL
names(meta) [names(meta) == "Chr38"]   <- "CHR"
names(meta) [names(meta) == "Pos_b37"] <- "POS37"
names(meta) [names(meta) == "Pos_b38"] <- "POS38"
head(meta)

#--------------------------------------------------------------#
#                      Step 4: Proxy SNPs                      #
#--------------------------------------------------------------#

cat("\n\nImporting the SNPs with strong LD ...\n\n")
snps <- read.delim(SNPs,  header = FALSE, col.names = c("CHR", "POS38"), sep= '\t', row.names = NULL, stringsAsFactors = F)

head(snps, 5)
nrow(snps)

# Find the row name of the last proxy SNP RPL3L, the locus before PDILT locus
RPL3L <- which(snps$POS38 == 1972803)

# Add missing SNP in PDILT locus after RPL3L locus
snps <- snps %>%
	add_row(CHR    = 16,
		POS38  = 20381010,
		.after = RPL3L)
		
# N. rows of proxy snps
cat("\n\n N. proxy snps is:", nrow(snps), "\n\n")

#print("...importing eGFR.Std.log Significant SNPs...")
#snpGFR <- read.table(snpSigGFR,  header=TRUE, sep= '\t', stringsAsFactors=F)
#head(snpGFR, 5)

#--------------------------------------------------------------#
#          Step 5: Extracting SNPs from GWAS results           #
#--------------------------------------------------------------#

cat("\n\nExtracting SNPs having strong LD from GWAS results ...\n\n")
snpGFRlog <- merge(gwas, snps, by = c("CHR", "POS38"), all = FALSE)

# Significant SNPs 
snpGFRlog <- snpGFRlog[snpGFRlog$Pvalue_CHRIS < 6.8e-4, ]

#sort by mpg (ascending) and cyl (descending)
#newdata <- mtcars[order(mpg, -cyl),] ##result[order(result[,"CHR"], result[,"BEG"]), ]

snpGFRlog <- snpGFRlog[order(snpGFRlog[,"CHR"]), ]

cat("\n\nThis is the merged file consisting of significants SNPs with strong LD ...\n\n")
head(snpGFRlog, 5)

#--------------------------------------------------------------#
#         Step 6: Extracting significant GWAS SNPs             #
#              from MetaSubtract GWAS results                  #
#--------------------------------------------------------------#

cat("\n\nNow merging significant CHRIS GWAS SNPs with MetaSubtract GWAS results ...\n\n")
result <- merge(meta, snpGFRlog, by = c("CHR", "POS38"), all.y = TRUE)

head(result, 5)

#--------------------------------------------------------------#
#                      Significant SNPs                        #
#--------------------------------------------------------------#

result$status_CHRIS  <- ifelse(result$Pvalue_CHRIS  < 6.8e-4, "Sig", "Non-Sig")
result$status_CKDGen <- ifelse(result$Pvalue_CKDGen < 6.8e-4, "Sig", "Non-Sig")
result <- result[order(result[,"CHR"]), ]
table(result$status_CKDGen)

#--------------------------------------------------------------#
#                       Saving output                          #
#--------------------------------------------------------------#

#write.table(snpGFRlog,			    file = "SigSNPs_eGFRw.log.Res.txt",	        sep = "\t", row.names = FALSE, quote=FALSE)
#write.table(result[result$status=="Sig",], file = "replicatedSNPs_eGFRw.log.Res.txt", 	sep = "\t", row.names = FALSE, quote=FALSE)

# saving the CKDGen and CHRIS merged data, teh replication result 
write.table(result, file = OUT1, sep = "\t", row.names = FALSE, quote = FALSE)

# The old name of the last of the above files was Meta_replicatedSNPs_eGFRw.log.Res.all.txt

# Sort the replication results based on CHR, POS
result_sorted <- result %>% arrange(CHR, POS)

# save the sorted result
write.table(result, file = OUT2, sep = "\t", row.names = FALSE, quote = FALSE)

#---------------#
# sort in the server if needed
#cat 10-Jan-23_Meta_replicatedSNPs_eGFRw.log.Res.all.txt | sort -g -k 1 > 10-Jan-23_replicatedSNPs_eGFRw.log.Res.txt

#---------------#
# run on the server
#sbatch --wrap 'Rscript extraction.R' --mem=65536 -J "repInCHRIS"
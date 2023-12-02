#!/usr/bin/Rscript
#---------------

#filtered_eGFR.Std.gz/ filtered_eGFR.log.Std.gz/ filtered_SerumCreatinine.log.Cleaned.gz/filtered_UACR.log.Std.gz

GWAS="/home/dghasemisemeskandeh/projects/gwas/Output/ReformedScheme/eGFRw.log.Res.txt.gz"
Meta="/home/dghasemisemeskandeh/projects/gwas/replicationAnalysis/scripts/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz"
BED38="/home/dghasemisemeskandeh/projects/gwas/replicationAnalysis/scripts/MetaBED38.txt"
snpSigGFR="/home/dghasemisemeskandeh/projects/gwas/pairways_LD/LD_range/replicatedSNPs_eGFR.Std.txt"
supptable_147_Loci="/home/dghasemisemeskandeh/projects/gwas/Data/Supptable.txt"
max.no.rows=100000

#---------------
library(tidyverse)

#--------------- Step1: GWAS --------------

cat("\n\nImporting the GWAS results ...\n\n")

names   <- c("CHR","POS38","END","MARKER_ID_CHRIS","NS_CHRIS", "AC_CHRIS","CALLRATE_CHRIS","GENOCNT_CHRIS","MAF_CHRIS","STAT_CHRIS","Pvalue_CHRIS","Beta_CHRIS","SE_CHRIS","R2_CHRIS")
gwas    <- read.table(GWAS, header = FALSE, col.names = names, sep = '\t', stringsAsFactors = F) #, nrows=max.no.rows

#drop unnecessary columns
gwas$END <- NULL
head(gwas, 5) ##names(gwas) [names(gwas)=="#CHROM"] <- "CHR"

#--------------- Step2: Subtracted Meta-GWAS --------------

cat("\n\nImporting the MetaSubtract results ...\n\n")
nameMeta <- c("Chr37", "Pos_b37", "RSID", "EA_CKDGen", "RA_CKDGen", "EAF_CKDGen", "Beta_CKDGen", "SE_CKDGen", "Pvalue_CKDGen", "n_total_sum_CKDGen")
metaB37  <- read.table(Meta,  header = FALSE, col.names = nameMeta, sep = '\t', stringsAsFactors = F) #, nrows=max.no.rows

head(metaB37)
#--------------- Step3: POS38 to be appended --------------

cat("\n\nImporting MetaSubtract lifted positions to GRCh38 ...\n\n")
POS 	<- read.table(BED38, header = TRUE, sep = '\t', stringsAsFactors = F, fill = TRUE) #, nrows=max.no.rows
head(POS)

#--------------- Step4: Adding POS38 to Subtracted Meta-GWAS --------------

cat("\n\nAppending the GRCh38 cooardiantes to MetaSubtract results ...\n\n")
meta <- merge(metaB37, POS, by.x = c("Chr37", "Pos_b37"), all = FALSE)

# Making column names consistant
meta$Chr37 <- NULL
names(meta) [names(meta) == "Chr38"]   <- "CHR"
names(meta) [names(meta) == "Pos_b37"] <- "POS37"
names(meta) [names(meta) == "Pos_b38"] <- "POS38"

head(meta)

#--------------- Step5: 147 CKDGen Loci for EA --------------

cat("\n\nImporting the 147 CKDGen Loci for EA ...\n\n")
EA_Loci <- read.delim(supptable_147_Loci,  header = TRUE, sep= '\t', row.names = NULL, stringsAsFactors = F)

head(EA_Loci, 5)

#selecting needed columns
EA_Loci <- EA_Loci[c("RSID", "Closest.Gene")]

#--------------- Step6: Extracting 147 Loci from Subtracted CKDGen EA Meta-GWAS results -------------

cat("\n\nExtracting 147 Loci from subtracted CKDGen EA Meta-GWAS results ...\n\n")
meta_Loci <- merge(EA_Loci, meta, by = c("RSID"), all = FALSE)

tail(meta_Loci, 5)

#--------------- Step7: Finding 147 Loci in CHRIS GWAS -------------

cat("\n\nNow merging 147 CKDGen Loci with CHRIS GWAS results ...\n\n")
result <- merge(meta_Loci, gwas, by = c("CHR", "POS38"), all.x = TRUE)

result <- result %>%
		mutate(status_CHRIS = ifelse(result$Pvalue_CHRIS  < 6.8e-4,
					     "Sig",
					     "Non-Sig")) %>%
		select(CHR, POS37, POS38, RSID, everything())

result <- result[order(result[,"CHR"]), ]

head(result, 5)
tail(result, 5)


#--------------- Significant Loci --------------

#Number of replicated loci in CHRIS
cat("\n\nNumber of replicated loci in CHRIS out of 147 EA CKDGen Loci is ...\n\n")

table(result$status_CHRIS)

#--------------- saving output --------------

#saving the CKDGen and CHRIS merged data, the replication result for 147 Loci
write.table(result, file = "31-Dec-22_Replication_of_147_EA_Loci_in_CHRIS_NOTsorted.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------

#sbatch --wrap 'Rscript extraction_147_Loci.R' --mem=65536 -J "147replicationInCHRIS"
#cat 31-Dec-22_Replication_of_147_EA_Loci_in_CHRIS_NOTsorted.txt | sort -g -k 1,2 > 31-Dec-22_Replication_of_147_EA_Loci_in_CHRIS.txt
#!/usr/bin/Rscript

#------------#
# Where to find the inputs and output

output.dir <- "/home/dghasemisemeskandeh/projects/gwas/01_subtract_ckdgen/output/"
CHRIS  <- "/home/dghasemisemeskandeh/projects/gwas/Data/CHRIScohort/CHRIS_EA_eGFR_overall_HRC_20170216.gwas.gz"
CKDGen <- "/home/dghasemisemeskandeh/projects/gwas/Data/CHRIScohort/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz"
#------------#

cat("\nReading CKDGen EA ...\n")
ckdgen <- read.table(CKDGen, header = TRUE, sep = ' ', stringsAsFactors = F, fill = TRUE)

cat("\nReading CHRIS 5K ...\n")
chris  <- read.table(CHRIS,  header = TRUE, sep = ' ', stringsAsFactors = F)

cat("\nHead of CKDGen EA ...\n\n")
tail(ckdgen)

cat("\nColumns of CKDGen EA ...\n\n")
str(ckdgen)

cat("\nHead of CHRIS 5K ...\n\n")
head(chris)

cat("\nColumns of CHRIS 5K ...\n\n")
str(chris)

# ----------------------------------------#
# Genomic-Control (GC) Creteria Estimation
# ----------------------------------------#

# Function to estimate the GC factor lambda 
# on the CKDGen meta-analysis results

cat("\nCalculating CKDGen EA Lambda ...\n\n")

#------------#
gc_lambda1 <- function(p)
  {
   # p = p-value vector
   z <- qnorm(0.5 * p,
	      mean = 0.0,
	      sd   = 1.0,
	      lower.tail = F)
   median(z * z, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = F)
  }

#------------#
# Lambda estimation
lambda_CKDGen <- gc_lambda1(ckdgen$P.value)
lambda_CHRIS  <- gc_lambda1(chris$pval)

#GIF <- 
cat("\nCKDGen EA Lambda is ", round(lambda_CKDGen, 3))
cat("\nCHRIS 5K Lambda is ",  round(lambda_CHRIS,  3))

# -----------------------------------------#
# Sebtracting CHRIS5K from CKDGen Meta-GWAS
# on the Overall eGFR of European Ancestry
# -----------------------------------------#

# install.packages("MetaSubtract")
library(MetaSubtract)

cat("\n\nSebtracting CHRIS 5K from CKDGen EA Meta-GWAS ...\n\n")

ckdgen_subtracted <- meta.subtract(
		metafile = CKDGen, 
		cohortfiles = CHRIS,
		lambda.meta = 1.05,
		lambdas.cohort = 0.99, 
		metamethod = "FIV",
		#metamethod="FSZ",
		gc_meta = FALSE, 
		calculate_lambda.meta = FALSE,
		calculate_lambdas.cohort = FALSE, 
		logfile = paste0(output.dir, "19-Dec-22_MetaSubtract_CKDGen_EA_CHRIS5K.log"),
		save.as.data.frame = TRUE,
		savefile = paste0(output.dir, "19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz"),
		dir = output.dir)

cat("\n\nHead of subtracted CKDGen EA ...\n\n")

head(ckdgen_subtracted)

#------------#
# save the subtracted meta-analysis summary stats
write.table(ckdgen_subtracted,
            paste0(output.dir, "19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt"),
	    sep = "\t",
	    row.names = FALSE,
	    quote = FALSE)

cat("\n\nSubtracted CKDGen EA was saved in ", output.dir)

#sbatch --wrap 'Rscript metasubtract.R' --mem=65536 -J "CkDGen_sub"
#!/usr/bin/Rscript


library(tidyverse)
#------------#
output.dir <- "/home/dghasemisemeskandeh/projects/gwas/01_subtract_ckdgen/output/"
Original   <- "/home/dghasemisemeskandeh/projects/gwas/Data/CHRIScohort/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz"
Subtract   <- paste0(output.dir, "19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz")

# Reading the files
ckdgen_original <- read.table(Original, header = TRUE, sep = ' ', stringsAsFactors = F, fill = TRUE)
ckdgen_subtract <- read.table(Subtract, header = TRUE, sep = '\t',stringsAsFactors = F)

head(ckdgen_original)
head(ckdgen_subtract)

# ----------------------------------------#
# Comparing CKDGen p-values before and 
# after subtracting CHRIS 5K from CKDGen
# ----------------------------------------#

# Merging CKDGen EA Meta-GWAS results with subtracted version
merged_CKDGen <- ckdgen_original %>%
		 inner_join(ckdgen_subtract,
			    by = c("Chr", "Pos_b37", "RSID", "n_total_sum"),
			    suffix = c("_Original", "_Subtracted"))

cat("\nHead of merged CKDGen EA before and after subtraction...\n")
head(merged_CKDGen)

#------------#
# qq-plot
ggplot(merged_CKDGen, aes(P.value_Original, P.value_Subtracted)) + 
       geom_abline(intercept = 0, slope = +1, color = "gray", size = 1) +
       geom_point(color = "skyblue", size = 0.5, alpha = 0.6) +
       geom_vline(xintercept = 0) + #, linetype="dashed", color = "darkgray", size=1.5)+
       geom_hline(yintercept = 0) + #, linetype="dashed", color = "darkgray", size=1.5)+
       xlab("Association coefficient P-Values for log(eGFRcrea) in CKDGen") +
       ylab("Association coefficient P-Values for log(eGFRcrea) in CKDGen subtracted from CHRIS 5K") +
       #coord_cartesian(xlim = c(0,1)) +
       theme_classic() +
       theme(axis.line = element_blank(),
	     panel.background = element_rect(fill = "white"),
             axis.text = element_text(size = 8),
             axis.title = element_text(size = 8, face = "bold"))

#------------#
# save the qq plot  
ggsave(paste0(output.dir, "19-Dec-22_Comparing effect size of the CKDGen variants before and after subtraction.png"),
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#------------#
# Pivoting the merged file
pivoted_CKDGen <- merged_CKDGen %>%
		  select(Chr, Pos_b37, P.value_Original, P.value_Subtracted) %>%
		  pivot_longer(cols = -c("Chr", "Pos_b37"),
			       names_to = c("trait", "dataset"),
			       names_pattern = "(.+)_(Original|Subtracted)$",
			       values_to = c("pval"))

head(pivoted_CKDGen)

#------------#
# density plot of the p-values
ggplot(pivoted_CKDGen, aes(pval, fill = dataset)) + 
       geom_density(alpha = 0.7) +
       scale_fill_manual(values = c("#E7B800", "#00AFBB")) + 
       scale_color_manual(values = c("#E7B800", "#00AFBB")) +
       labs(fill = "CKDGen",
	    x = "Association coefficient P-Values for log(eGFRcrea)") + 
       theme_classic() +
       theme(axis.line = element_blank(),
	     panel.background = element_rect(fill = "white"),
	     axis.text = element_text(size =10),
	     axis.title = element_text(size = 8, face = "bold"),
	     legend.position = c(0.88, 0.88))

#------------#
# save the density plot of the p-values
ggsave(paste0(output.dir, "19-Dec-22_Comparing distribution of the CKDGen variants before and after subtraction.png"),
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#sbatch --wrap 'Rscript metasubtractQC.R' --mem=65536 -J "CkDGen_sub_comp"
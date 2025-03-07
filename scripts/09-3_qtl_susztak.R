#==================================#
#      eQTLs in kidney tissue
#          Feb 2, 2025
#==================================#


library(tidyverse)
library(data.table)
library(R.utils) # to be able to read zip file
source("00-0_configuration.R")


#------------#
# The script aims to check if the variants for which we found significant mediation effect are eQTLs.

# This table provides 1,179,179 significant SNP~gene pairs identified after 
# meta-analysis of four eQTL studies (Sheng et al., Ko et al., GTEx, NephQTL).
# significant cis-eQTLs derived from meta-analysis (N = 686).
# GRCh37
# https://susztaklab.com/Kidney_eQTL/download.php


#------------#
# inputs: results of meta-analysis GWAS
path_meta  <- "~/projects/HaploReg/data/Kidney_eQTL_Meta_S686_Significant.q0.01.txt.gz"


# outputs
out_tbl_kidney <- "~/projects/kidneyInCHRIS/outputs/03-Mar-25_eqtls_susztak.csv"
out_plt_kidney <- "~/projects/kidneyInCHRIS/outputs/03-Mar-25_eqtls_susztak.png"


#------------#
# read meta-analysis results
meta <- data.table::fread(path_meta)

# specify variants in Table 3 of mediation paper
snps_in_table3 = c("rs28394165", "rs10025351", "rs28817415", "rs4859682", "rs13146355", "rs3812036")


#------------#
# take variants in Table 3 from meta-analysis
eqtls_kidney <- meta %>%
  dplyr::filter(RSID %in% snps_in_table3) %>%
  dplyr::rename(
    chr_pos37 = SNP_POS,
    rsid = RSID,
    assoc_gene = GeneSymbol,
    beta  = Effect2ALT,  # do NOT change the effect direction as it refers to EA.
    pvalue = PVAL,
    ) %>%
  dplyr::mutate(
    EA_OA = paste0(ALT, "/", REF) # flip alleles to have effect allele first
    )

#------------#
# save eqtls in kidney as cCREs in paper
write.csv(eqtls_kidney, file = out_tbl_kidney, row.names = F, quote = F)


#----------------------------------------#
#-----          eQTLs Plot          -----
#----------------------------------------#

# plot eQTLs in kidney human atlas
plt_kidney <- eqtls_kidney %>%
  ggplot(aes(x = gene_coder(rsid), y = assoc_gene, colour = beta)) +
  geom_point(size = 20) +
  scale_y_discrete(position = "right") +
  scale_color_gradient2(
    low="blue", mid="#AF7AC5", high="#E74C3C", guide="colorbar", midpoint=0, na.value=NA
  )+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size=10, face=2, angle=90, hjust=0, vjust=0.5),
    axis.text.y = element_text(size=10, face=2),
    legend.position  = c(1.05, -.07)
  )

#----------#
# save the plot
ggsave(filename = out_plt_kidney, plt_kidney, width = 8, height = 5.5, dpi = 150, units = "in", bg = "white")

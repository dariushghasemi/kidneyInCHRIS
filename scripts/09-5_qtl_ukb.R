

# libraries
library(tidyverse)
library(data.table)
library(readxl)
source("00-0_configuration.R")


#----------#
# inputs
path_literature <- "~/projects/kidneyInCHRIS/inputs/literature_file.xlsx"

# outputs
out_tbl_ukb <- "~/projects/kidneyInCHRIS/outputs/03-Mar-25_pqtls_ukb.csv" 
out_plt_ukb <- "~/projects/kidneyInCHRIS/outputs/03-Mar-25_pqtls_ukb.png"


# specify variants in Table 3 of mediation paper
snps_in_table3 = c("rs28394165", "rs10025351", "rs28817415", "rs4859682", "rs13146355", "rs3812036")
snps_in_tab3   = c("4:76472865_T_C", "4:76472942_C_T", "4_76480299_C_T", "4_76489165_C_A", "4_76490987_G_A", "5_177386403_C_T")


#----------#
# function to read excel file with multiple sheets
pqtl_studies <- lapply(
  readxl::excel_sheets(path_literature), 
  function(x) read_excel(path_literature, sheet = x)
)


#----------#
# mapping file including gene and protein names to merge
mapping_olink   <- pqtl_studies[9] %>% as.data.frame()
mapping_protein <- pqtl_studies[8] %>% as.data.frame()

# check if multiple rows with the same UniProt 
# do NOT differ in key columns, so we can use unique rows for merge
#mapping_protein %>% dplyr::filter(UniProt %in% pqtls_ukb$UniProt)


#----------#
# proteins QTLs in UKB as cCREs in paper
pqtls_ukb <- pqtl_studies[4] %>% 
  as.data.frame() %>%
  dplyr::filter(rsID %in% snps_in_table3 | rsID %in% snps_in_tab3) %>% # till here for checking UniProt ids
  left_join(
    mapping_protein %>% select(UniProt, Target_longname, Target_Name) %>% distinct(),
    join_by(UniProt)
  ) %>%
  left_join(
    mapping_olink %>% select(UniProt, Target_Name, Entrez_Gene_Name), # better to use protein names in Olink dataset
    join_by(UniProt),
    suffix = c("_protein", "_olink")
  ) %>%
  dplyr::mutate(
    gene_prot = str_c(Entrez_Gene_Name, " (", Target_Name_protein, ")"),
    gene_prot = if_else(is.na(gene_prot), str_c(Entrez_Gene_Name, " (", Entrez_Gene_Name, ")"), gene_prot),
    rsid_locus = gene_coder(rsID),
    gene_prot_ord = fct_reorder(UniProt, pos38)
  )


# save pqtls in UKB
write.csv(pqtls_ukb, file = out_tbl_ukb, row.names = F, quote = F)


#----------------------------------------#
#-----          pQTLs Plot          -----
#----------------------------------------#

# draw pQTLs plot
plt_ukb <- pqtls_ukb %>%
  ggplot(aes(x = rsid_locus, y = gene_prot_ord, colour = BETA)) +
  geom_point(size = 8) +
  scale_y_discrete(position = "right") +
  scale_color_gradient2(
    low="blue", mid="#AF7AC5", high="#E74C3C", guide="colorbar", midpoint=0, na.value=NA
  )+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size=10, face=2, angle=90, hjust=0, vjust=0.5),
    axis.text.y = element_text(size=10, face=2),
    legend.position  = c(1.45, -.07)
  )


#----------#
# save the plot
ggsave(filename = out_plt_ukb, plt_ukb, width = 8, height = 5.5, dpi = 150, units = "in", bg = "white")


# libraries
library(tidyverse)

#----------#
# inputs
path_query <- "/scratch/dariush.ghasemi/projects/decode/results"

# outputs
out_combined <- paste0(path_query, "/combined_quey_results.csv")
out_tbl_decode <- "~/projects/kidneyInCHRIS/outputs/03-Mar-25_pqtls_decode.csv" 
out_plt_decode <- "~/projects/kidneyInCHRIS/outputs/03-Mar-25_pqtls_decode.png"

#----------#
# find list of query results separated by task number
files_query <- list.files(
  path_query,
  pattern = "results_([0-9]+).tsv", 
  recursive = FALSE, 
  full.names = TRUE
  )

#----------#
# read list of the results
query_results <- lapply(
  files_query,
  data.table::fread, 
  header = T, 
  sep = "\t"
) #%>% bind_rows()


# merge them
query_combined <- do.call(rbind, query_results)

# save pqtls in decode as cCREs in paper
write.csv(query_combined, file = out_combined, row.names = F, quote = F)


#----------#
# cCREs in paper
query_signif <- query_combined %>%
  # take significant associations
  dplyr::filter(Pval < 5e-08) %>%
  # prepare columns for plot
  dplyr::mutate(
    protein_name = str_remove_all(
      filename, "Proteomics_SMP_PC0_(\\d+)_(\\d+)_([A-Za-z0-9]+)_" # fix failed matching for protein name
      ) %>% 
      str_remove("_10032022.txt.gz"),
    gene_prot = str_c(gene_name, " (", protein_name, ")"),
    rsid_locus = gene_coder(rsids),
    gene_prot_ord = fct_reorder(gene_prot, Chrom)
  )

# save pqtls in decode as cCREs in paper
write.csv(query_signif, file = out_tbl_decode, row.names = F, quote = F)

#----------#
# draw QTLs plot
plt_decode <- query_signif %>%
  ggplot(aes(x = rsid_locus, y = gene_prot_ord, colour = Beta)) +
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
    legend.position  = c(1.5, -.07)
  )


#----------#
# save the plot
ggsave(filename = out_plt_decode, plt_decode, width = 8, height = 5.5, dpi = 150, units = "in", bg = "white")


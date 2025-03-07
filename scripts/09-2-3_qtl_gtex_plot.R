

# libraries
library(tidyverse)
library(data.table)
source("00-0_configuration.R")


#----------#
# inputs
path_query_getx <- "/scratch/dariush.ghasemi/projects/gtex_v10/02-Mar-25_summary.tsv"

# outputs
out_plt_gtex <- "/home/dariush.ghasemi/projects/kidneyInCHRIS/outputs/03-Mar-25_eqtls_gtex_v10.png"
out_tbl_gtex <- "/home/dariush.ghasemi/projects/kidneyInCHRIS/outputs/03-Mar-25_eqtls_gtex_v10.csv"



#----------#
# read query results
res_gtex <- data.table::fread(path_query_getx) %>%
  dplyr::filter(P_Value < 5e-08)



# save eqtls in kidney as cCREs in paper
write.csv(res_gtex, file = out_tbl_gtex, row.names = F, quote = F)


#----------------------------------------#
#-----          eQTLs Plot          -----
#----------------------------------------#

# draw eQTLs plot
plt_gtex_v10 <- res_gtex %>%
  ggplot(aes(x = fct_reorder(Gene, POS), y = fct_rev(Tissue), colour = NES)) +
  geom_point(size = 3) +
  scale_y_discrete(position = "right") +
  facet_grid(~ gene_coder(rs_id)) +
  scale_color_gradient2(
    low="blue", mid="#AF7AC5", high="#E74C3C", guide="colorbar", midpoint=0, na.value=NA
  )+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size=9, face="bold.italic", angle=90, hjust=0, vjust=0.5),
    axis.text.y = element_text(size=10, face=2),
    legend.position  = c(1.25, -0.1)
  )


#----------#
# save the plot
ggsave(filename = out_plt_gtex, plt_gtex_v10, width = 12, height = 8.5, dpi = 200, units = "in", bg = "white")


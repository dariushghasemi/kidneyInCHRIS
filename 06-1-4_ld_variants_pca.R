#!/usr/bin/Rscript

#=========================================#
# PCA on the in-LD variants with GWAS hits
#=========================================#

# This script was initiated on November 28, 2023.
# The data here correspond TOPMedR2 imputed CHRIS10K.
# The genomic build is GRCh38.

#-----------------------------------------------------#
#------     Variance explained in 147 loci      ------
#-----------------------------------------------------#

library(tidyverse)

# Read the dosage file of the variants 
# in strong LD (r2 > 0.8) generated in
# `06-1-3_ld_variants_extraction.sh`
inLD_variants <- read.delim("27-Nov-23_variants_in_ld_dosage.txt", header = T)

#------------#
# reshaping the dosage file to have SNPs in columns
inLD_wide <- inLD_variants %>%
  pivot_wider(id_cols = AID,
              names_from = SNP_ID, 
              values_from = Dosage)

# run PCA on dosage levels of the variants
inLD_PCA <- prcomp(inLD_wide %>% select(- AID), scale. = T)

#------------#
# representing the cumulative variances explained by PCs
#install.packages("factoextra")
library(factoextra)

plot_pca <- fviz_eig(inLD_PCA, addlabels = F, ncp = 40, ylim = c(0, 15)) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  theme_classic() +
  theme(title = element_blank())

ggsave("06-Dec-23_pcs_vs_total_varinace_explained.png", 
       plot_pca, width = 12, height = 5.5, dpi = 300, units = "in")

#------------#
# Calculate the cumulative variance explained by each component
cumulative_var <- cumsum(inLD_PCA$sdev^2) / sum(inLD_PCA$sdev^2)

# Create a data frame for ggplot
expl_variance <- data.frame(Components = 1:length(cumulative_var),
                            Cumulative_Variance = cumulative_var)

# save the variance (standard deviations^2) explained by PCs
write.table(expl_variance,
            "06-Dec-23_cumulative_variance_explained_by_each_pc_147_loci.txt",
            quote = F, row.names = F)

#------------#
# Plot the cumulative variance explained using ggplot
um_var_plot <- function(df){
  
  df %>%
    ggplot(aes(x = Components, y = Cumulative_Variance)) +
    geom_line() +
    geom_point(color = "steelblue", size = 2) +
    geom_point(aes(x = 147, y = .98), color = "red2", size = 2.5) +
    #geom_hline(yintercept = 0.99, linetype = "dashed", color = "red") +
    #geom_vline(xintercept = 11,   linetype = "dashed", color = "purple") +
    scale_x_continuous(breaks = c(0, 200, 500, 1000, 2000, 4000, 6000), limits = c(0, 6500)) +
    #scale_x_continuous(breaks = seq(0, 165, 15), limits = c(0, 165)) +
    scale_y_continuous(breaks = seq(0, 1, .1), limits = c(0, 1)) +
    scale_color_manual(values = c(NA, "red2")) +
    labs(x = "\nNumber of Principal Components",
         y = "Cumulative Variance Explained\n") +
    theme_classic() +
    theme(axis.title  = element_text(face = "bold", size = 13),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10))
}

#------------#
# cumulative variance plot
expl_variance %>% um_var_plot()

# save the plot
ggsave("06-Dec-23_No_of_pcs_explaining_total_varinace_147_loci.png", 
       last_plot(), width = 9, height = 7.5, dpi = 300, units = "in")

#------------#
# Find the number of components needed for 95% and 99% cumulative variance explained
components_95 <- which(cumulative_var >= 0.95)[1]
components_99 <- which(cumulative_var >= 0.99)[1]

# Print the results
cat("Number of components for 95% variance:", components_95, "\n")
cat("Number of components for 99% variance:", components_99, "\n")

quit()
#-----------------------------------------------------#
#------      Variance explained in 11 loci      ------
#-----------------------------------------------------#

# run PCA on dosage levels of the 163 
# replicated variants fo 10146 individuals
repSNPs_PCA <- prcomp(vcfReg[targets], scale. = T)

# Calculate the cumulative variance explained by each component
cum_var_repsnps <- cumsum(repSNPs_PCA$sdev^2) / sum(repSNPs_PCA$sdev^2)

# Create a data frame for ggplot
expl_var_repsnps <- data.frame(Components = 1:length(cum_var_repsnps),
                               Cumulative_Variance = cum_var_repsnps) 

# save the variance (standard deviations^2) explained by PCs
write.csv(expl_var_repsnps,
          "05-Dec-23_cumulative_variance_explained_by_each_pcs_163_repsnps.csv",
          quote = F, row.names = F)

#------------#
# PCs explaining variability across 163 replicated variants
fviz_eig(repSNPs_PCA, addlabels = F, ncp = 40) +
  scale_y_continuous(breaks = seq(0, 55, 5)) +
  theme_classic() +
  theme(title = element_blank(),
        axis.title = element_text(face = "bold", size = 12))


ggsave("05-Dec-23_pcs_vs_percentage_of_explained_varinace_163_repsnps.png", 
       last_plot(), width = 12, height = 5.5, dpi = 300, units = "in")

#------------#
# cumulative variance plot
expl_var_repsnps %>% um_var_plot()

# save the plot
ggsave("05-Dec-23_No_of_pcs_explaining_total_variance_163_resnps.png", 
       last_plot(), width = 9, height = 7.5, dpi = 300, units = "in")

#------------#
# sbatch --wrap 'Rscript 06-1-4_ld_variants_pca.R'  -c 4 --mem-per-cpu=8GB -J "06-1-4.R"



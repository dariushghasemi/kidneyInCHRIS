#!/usr/bin/Rscript

#=========================================#
# PCA on the in-LD variants with GWAS hits
#=========================================#

# This script was initiated on November 28, 2023.
# The data here correspond TOPMedR2 imputed CHRIS10K.
# The genomic build is GRCh38.

#------------#
library(tidyverse)

# read the dosage file of the variants in strong LD (r2 > 0.8)
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

plot_pca <- fviz_eig(inLD_PCA, addlabels = F, ncp = 6100, ylim = c(0, 15)) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  theme_classic() +
  theme(title = element_blank())

ggsave("30-Nov-23_pcs_vs_total_varinace_explained.png", 
       plot_pca, width = 12, height = 6.5, dpi = 300, units = "in")

#------------#
# Calculate the cumulative variance explained by each component
cumulative_var <- cumsum(inLD_PCA$sdev^2) / sum(inLD_PCA$sdev^2)

# Create a data frame for ggplot
expl_variance <- data.frame(Components = 1:length(cumulative_var),
                           Cumulative_Variance = cumulative_var)

# save the variance (standard deviations^2) explained by PCs
write.table(expl_variance,
            "30-Nov-23_cumulative_variance_explained_by_each_pc.txt",
            quote = F, row.names = F)

#------------#
# Plot the cumulative variance explained using ggplot
ggplot(expl_variance, aes(x = Components, y = Cumulative_Variance)) +
  geom_line() +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(0, 200, 500, 1000, 2000, 4000, 6000)) +
  scale_y_continuous(breaks = c(0, .25, .50, .75, .95, 0.99, 1)) +
  labs(x = "Number of Components", y = "Cumulative Variance Explained") +
  theme_classic() +
  theme(axis.title  = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 6))

# save the plot
ggsave("30-Nov-23_No_of_pcs_explaining_total_varinace.png", 
       last_plot(), width = 9, height = 6.5, dpi = 300, units = "in")

#------------#
# Find the number of components needed for 95% and 99% cumulative variance explained
components_95 <- which(cumulative_var >= 0.95)[1]
components_99 <- which(cumulative_var >= 0.99)[1]

# Print the results
cat("Number of components for 95% variance:", components_95, "\n")
cat("Number of components for 99% variance:", components_99, "\n")

#------------#
# sbatch --wrap 'Rscript 05-4_pca_in-LD_variants.R'  -c 4 --mem-per-cpu=8GB -J "05-4.R"



#=========================================#
# PCA on the in-LD variants with GWAS hits
#=========================================#

# This script was initiated on November 28, 2023.

inLD_variants <- read.delim("C:/Users/dghasemisemeskandeh/Downloads/dosage.test",
                             header = T)

# reshaping the dosage file to have SNPs in columns
inLD_wide <- inLD_variants %>%
  pivot_wider(id_cols = AID,
              names_from = SNP_ID, 
              values_from = Dosage)

# run PCA on dosage levels of the variants
inLD_PCA <- prcomp(inLD_wide %>% select("chr1:10599281:C:G", "chr22:42744714:T:G"), scale. = T)

# Calculate the cumulative variance explained by each component
cumulative_var <- cumsum(inLD_PCA$sdev^2) / sum(inLD_PCA$sdev^2)
sum(cumulative_var <= 0.95)

# Find the number of first components explaining >80% of the cumulative variance
message <- paste(sum(cumulative_var <= 0.95), 
                 "principal components are needed to explain %95 of the total variance of the",
                 length(inLD_PCA$sdev),
                 "in-LD variants.")

# show the number PCA explaining %95 of teh total variants
cat("\n", message, "\n")

# save the variance (standard deviations^2) explained by PCs
write.table(inLD_PCA$sdev^2, "28-Nov-23_variance_explained_by_pcs.txt", quote = F)
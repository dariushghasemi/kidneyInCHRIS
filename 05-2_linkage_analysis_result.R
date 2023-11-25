#=========================================#
# Categorizing the variants in strong LD with GWAS hits
#=========================================#

# Started this script on 24-Nov-23 
# after receiving the first revision of the paper


# read the result of linkage analysis: LD between variants in proxy of CKDGen GWAS hits
result_linkage <- read.delim("F:/Dariush/PhD/Analysis/1st Paper/Outputs_tables/linkageAnalysis/test_convert.txt", header = T)

# read the 147 CKDGen identified eGFRcrea-associated loci for EU Ancestry
result_ckdgen <- read.delim("F:/Dariush/PhD/Analysis/1st Paper/Outputs_tables/Replication_eGFRw.log.Res_MultiA/Supptable.txt", header = T)

# adding locus name to the variants in LD

result_ckdgen %>% 
  select(CHR, BEG, RS.number, Closest.Gene) %>%
  rename(rsid = RS.number, locus = Closest.Gene) %>%
  mutate(CHR = as.character(CHR)) %>%
  right_join(result_linkage %>% select(CHR, POS1, POS2, R), by = c("CHR", "BEG" = "POS1")) %>% head()


intersect(unique(result_linkage$POS1), result_ckdgen$BEG)

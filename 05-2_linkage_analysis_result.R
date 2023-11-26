#=========================================#
# Categorizing the variants in strong LD with GWAS hits
#=========================================#

# Started this script on 24-Nov-23 
# after receiving the first revision of the paper


# read the result of linkage analysis (built 38):
# LD between variants in proxy of CKDGen GWAS hits
result_linkage <- read.delim("F:/Dariush/PhD/Analysis/1st Paper/Outputs_tables/linkageAnalysis/test_convert.txt",
                             header = T, colClasses = c("character", "character", "character", "numeric", "numeric", "numeric", "numeric"),
                             stringsAsFactors = FALSE)

# read the 147 CKDGen identified leading signals
# eGFRcrea-associated loci for EU Ancestry (built 37)
result_ckdgen <- read.delim("F:/Dariush/PhD/Analysis/1st Paper/Outputs_tables/Replication_eGFRw.log.Res_MultiA/Supptable.txt",
                            header = T, colClasses = "character", stringsAsFactors = FALSE)

# read the coordinates of meta-analysis summary stats both built 37 & 38
bed38 <- read.delim("F:/Dariush/PhD/Analysis/Data/MetaBED38.txt",
                    header = T, colClasses = c("character", "character", "character", "character"), stringsAsFactors = FALSE)

# adding locus name to the variants in LD
b38_ckdgen <- result_ckdgen %>% 
  select(CHR, BEG, RS.number, Closest.Gene) %>%
  rename(rsid = RS.number, locus = Closest.Gene) %>%
  mutate(BEG = as.character(BEG)) %>%
  # adding positions in GRCh38 to CKDGen results
  right_join(bed38, by = c("CHR" = "Chr37", "BEG" = "Pos_b37")) %>%
  mutate(Chr38 = str_c("chr", Chr38))


ld_varinats <- result_linkage %>%
  as_tibble() %>%
  select(CHR, POS1, POS2, R) %>%
  # add rsID and locus name of POS2
  left_join(b38_ckdgen %>% select(Chr38, Pos_b38, rsid, locus),
            by = c("CHR" = "Chr38", "POS2" = "Pos_b38"),
            relationship = "many-to-many") %>%
  # add rsIDs belonging to POS2
  rename(rsid_POS2 = rsid, locus_POS2 = locus) %>%
  # add rsID and locus name of POS1
  left_join(b38_ckdgen %>% select(Chr38, Pos_b38, rsid, locus),
            by = c("CHR" = "Chr38", "POS1" = "Pos_b38"),
            relationship = "many-to-many") %>%
  # add rsIDs belonging to POS1
  rename(rsid_POS1 = rsid) %>%
  relocate(locus) %>%
  mutate(locus = if_else(is.na(locus), locus_POS2, locus)) %>%
  select(- locus_POS2)


# list of variants in strong LD with hits
ld_varinats %>%
  mutate(SNP1 = str_c(CHR, "_", POS1),
         SNP2 = str_c(CHR, "_", POS2)) %>% 
  select(SNP1, SNP2) %>%
  unlist() %>%
  unique() %>%
  sort()

# create a variants file for extraction from VCF
ld_varinats %>%
  select(CHR, POS1, POS2) %>%
  pivot_longer(cols = c(POS1, POS2), names_to = "Position") %>% 
  distinct(CHR, value) %>%
  pivot_wider(names_from = Position, values_from = value)


ld_varinats %>%
  #filter(POS1 == 10670853 | POS2 == 10670853)
  #filter(locus == "SHROOM3")
  count(locus)

repSNPs_tmp %>% as_tibble() %>%
  # 0.05/(6337/2) = 1.578e-05
  filter(Pvalue_CHRIS <= 1.578e-5) %>%
  # 0.05/(6337) = 8.03471e-06
  filter(Pvalue_CHRIS <= 7.890e-6) %>%
  count(Locus)

# A tibble: 2 × 2
# Locus       n
# DDX1        2
# PIP5K1B    12

# A tibble: 1 × 2
# Locus       n
# PIP5K1B    12


intersect(unique(result_linkage$POS2), result_ckdgen$BEG)
length(
  unique(
    intersect(unique(result_linkage$POS1),
              unique(b38_ckdgen$Pos_b38)
              )))

# heatmap
ggplot(result_linkage[1:100,], aes(POS1, POS2, color = R)) +
  geom_tile()


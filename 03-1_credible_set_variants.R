#=========================================#
# Haplotype Analysis: Credible set variant
#=========================================#


#-----------------------------------------------------#
#---------------- Haplotype Analysis -----------------
#-----------------------------------------------------#

library(haplo.stats)

# Export chris phenotype data for haplotype Regression analysis
vcfReg_Mag %>% 
  select(AID, Age, Sex, eGFRw.log.Res, starts_with("PC")) %>% 
  write.table("06-Nov-2022_chris4HaploReg.txt",
              sep = "\t", row.names = FALSE, quote=FALSE)
#----------#

# Achieving hot region based on recombination rate
# for rsid 4:76480299 in SHROOM3 in GRCh38 -> n=10,758
recombSHROOM3 <- read.table("D:\\Dariush\\PhD\\Analysis\\Data\\recomb-hg38\\recomb-hg38\\plink.GRCh38.map\\plink.chr4.GRCh38.map", header = FALSE, sep = " ")
#----------#

# Genotype data
vcfSHROOM3 <- 
  read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\chr4-SHROOM3_DS.txt",
             col.names = c("AID", "CHR", "POS", "MARKER_ID", "REF", "ALT", "Dosage"), 
             #nrows = 107579, 
             sep = "\t", stringsAsFactors = FALSE) #%>% count(POS)
#----------#

vcfSHROOM3_long <-
  vcfSHROOM3 %>%
  #filter(POS > 76251700 & POS < 76516500) %>% #View()
  #mutate_at("MARKER_ID", str_replace_all, ":[A-T-C-G]+:[A-T-C-G]+", "")%>%
  select(AID, POS, MARKER_ID, Dosage) %>%
  rename(SNPid = MARKER_ID) %>% 
  pivot_wider(id_cols = AID,
              names_from = SNPid, 
              values_from = Dosage) %>% 
  inner_join(.,
             chris[c("AID",
                     "Sex",
                     "Age",
                     "eGFRw.log.Res")],
             by="AID")
#slice_head(n = 1000)#filter(AID != "0010197321")#%>% head()
#----------#

recombSHROOM3 %>% #head(100) %>% View() 
  rename(CHR = V1, POS = V4, recombRate = V3) %>% 
  ##filter(POS > 76213540 & POS < 76783753) %>% #dim()
  #filter(POS > 76214040 & POS < 76850000) %>% #dim()
  #filter(POS > 76480299 - 500000 & POS < 76480299 + 500000) %>% 
  #beg: FAM47E gen,chr4:76,251,700; end: SOWAHB gene, chr4:76,894,152
  #filter(POS > 76251700 & POS < 76894151) %>%
  #1st recomb hotspot = 76516500, 2nd=
  filter(POS > 76251700 & POS < 76516500) %>% #dim()
  mutate(diff = recombRate - lag(recombRate, 
                                 default = first(recombRate))) %>%
  #top_n(n = 100, wt = recombRate)
  ggplot(aes(POS, diff)) +
  geom_line(color = "steelblue4") +
  # scale_x_continuous(minor_breaks = seq(76255000, 76900000, 50000), 
  #                    breaks       = seq(76255000, 76900000, 50000))+
  labs(x = "Position", y = "Recombination Rate")+
  theme_light() +
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text.x  = element_text(size = 5,  face="bold"),
        axis.text.y  = element_text(size = 5,  face="bold"))
#----------#

# Clinically relevant variants in SHROOM3 locus
read.csv("C:\\Users\\dghasemisemeskandeh\\Documents\\haplotype_regression\\gnomAD_v3.1.2_4-76214040-76435228_20230312.csv", na.strings = c("", "NA")) %>%
  rbind(read.csv("C:\\Users\\dghasemisemeskandeh\\Documents\\haplotype_regression\\gnomAD_v3.1.2_SHROOM3.csv", na.strings = c("", "NA"))) %>%
  janitor::clean_names() %>%
  count(vep_annotation) %>%
  mutate(vep_annotation = str_trim(vep_annotation, side = "both"),
         vep_annotation = str_replace(vep_annotation, "_variant", ""),
         p = n/sum(n),
         structure = reorder(vep_annotation, n)) %>%
  filter(!str_detect(vep_annotation, "intron|NA")) %>%
  ggplot(aes(x = structure, y = n, fill = vep_annotation)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           #mapping = aes(x = , y = n),
           show.legend = FALSE,
           width = 0.7,
           fill = "steelblue2",
           color = "grey50") +
  labs(x = "variant", y = "count") + #"percentage"
  coord_flip() +
  theme_classic()

ggsave("16-Mar-23_Frequency of structural variants in SHROOM3 locus_count.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")




#-----------------------------------------------------#
#                    Example run
#-----------------------------------------------------#

data(hla.demo) 
names(hla.demo)
#----------#

# dosage level to 2 columns: major(1) & minor(2) alleles
geno1 <- matrix(c(0,0,1,
                  1,0,2,
                  2,1,0), ncol=3, byrow=TRUE)
geno1to2(geno1, locus.label=c("A", "B", "C"))

## demonstrate how NA and 3 will be coded
geno1[1,3] <- NA
geno1[1,1] <- 3
geno1to2(geno1)

# set up data for haplo.glm, include geno.glm,
# covariates age and male, and responses resp and y.bin
geno <- hla.demo[,c(17,18,21:24)]
label <-c("DQB","DRB","B")

# converting dosage level from one column to 2 columns
geno1to2(geno, locus.label=c("A", "B", "C"))

# The hla.demo data already had alleles in two columns for each locus. For many SNP datasets, the data
# is in a one column format, giving the count of the minor allele. To assist in converting this format to two
# columns, a function named geno1to2 has been added to the package. See its help file for more details

geno.glm <- setupGeno(geno, miss.val = c(0,NA), locus.label = label)
attributes(geno.glm)

attach(hla.demo)
detach(hla.demo)

y.bin <- 1*(resp.cat == "low")
glm.data <- data.frame(geno.glm, 
                       age = hla.demo$age, 
                       male = hla.demo$male,
                       y = hla.demo$resp, 
                       y.bin = y.bin)
#----------#
# glm fit with haplotypes, additive gender covariate on gaussian response
fit.gaus <- haplo.glm(y ~ male + geno.glm,
                      family = gaussian,
                      data   = glm.data,
                      na.action = "na.geno.keep", 
                      locus.label = label, 
                      x = TRUE,
                      control = haplo.glm.control(haplo.freq.min = .02))

summary(fit.gaus)

#-----------------------------------------------------#
#                      SHROOM3
#-----------------------------------------------------#

# Genotypes
SHROOM3.geno <- 
  vcfSHROOM3_long %>% 
  #select(starts_with("chr")) %>%
  select(
    "chr4:76300094:C:T", #existed in CKDGen
    "chr4:76421860:A:G",
    "chr4:76427807:G:A",
    "chr4:76436625:G:A",
    "chr4:76438522:T:C",
    "chr4:76438524:C:T",
    "chr4:76438636:G:A",
    "chr4:76439053:C:T",
    "chr4:76439083:G:A",
    "chr4:76439833:G:A",
    "chr4:76440368:T:C",
    "chr4:76440454:G:A",
    "chr4:76442486:T:G",
    "chr4:76442973:A:G",
    "chr4:76443244:A:T",
    "chr4:76446201:G:T",
    "chr4:76446251:G:A",
    "chr4:76447040:C:T",
    "chr4:76447317:T:C",
    "chr4:76448885:C:T",
    "chr4:76449125:C:T",
    "chr4:76457149:T:C",
    "chr4:76458036:G:A",
    "chr4:76458504:C:T",
    "chr4:76462648:G:C",
    "chr4:76463426:A:G",
    "chr4:76464644:A:G",
    "chr4:76480299:C:T", #CKDGen
    "chr4:76489165:C:A", #Wuttke
    "chr4:76490987:G:A" #Wuttke
  ) %>%
  as.matrix()

# changing the format ina way to have two columns for each varint
SHROOM3 <- geno1to2(round(SHROOM3.geno,0),
                    locus.label = colnames(SHROOM3.geno))
#----------#
# Haplotype frequency
save.em <- haplo.em(geno = SHROOM3, 
                    locus.label = colnames(SHROOM3.geno), 
                    miss.val = c(0,NA))

print(save.em, nlines=10)
summary(save.em, nlines=7)

# show full haplotypes, instead of codes
summary(save.em, show.haplo=TRUE, nlines=7)
#----------#

# Fitting regression model
SHROOM3.glm <- setupGeno(SHROOM3,
                         miss.val = c(0,NA),
                         locus.label = colnames(SHROOM3.geno))
attributes(SHROOM3.glm)

# GLM data
SHROOM3.data <- data.frame(
  SHROOM3.glm,
  "Age"  = vcfSHROOM3_long$Age,
  "Male" = vcfSHROOM3_long$Sex,
  "eGFR" = vcfSHROOM3_long$eGFRw.log.Res)

# glm fit with haplotypes, additive gender covariate on gaussian response
fit.gaus <- haplo.glm(eGFR ~ Age + Male + SHROOM3.glm,
                      family = gaussian,
                      data   = SHROOM3.data,
                      na.action = "na.geno.keep", 
                      locus.label = colnames(SHROOM3.geno),
                      x = TRUE,
                      control = haplo.glm.control(haplo.freq.min = .02))#haplo.min.count=10

SHROOM3_haplo <- summary(fit.gaus)$haplotypes

SHROOM3_haplo %>%
  mutate(Haplotype = row.names(SHROOM3_haplo)) %>% 
  pivot_longer(cols = -c (hap.freq, Haplotype),
               names_to = "SNP",
               values_to = "Allele") %>% 
  ggplot(aes(SNP, Haplotype, color = Allele)) +
  geom_point(size = 9, alpha = .75) +
  scale_color_discrete(name = "Allele",
                       type = c("white", "grey60", "#0099CC"),
                       labels = c("*", "Major", "Minor")) +
  #scale_shape_manual(values = c(8, 21, 21)) +
  #guides(color = TRUE, shape = FALSE) +
  scale_x_discrete(position = 'top') +
  scale_y_discrete(limits = rev) +
  labs(x    = "", 
       y    = "") +
  theme_classic() +
  theme(legend.title= element_text(size=7, face="bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text  = element_text(size = 8),
        axis.text.x = element_text(size=8, face="italic", angle=35, vjust=1.05, hjust=0.05),
        axis.text.y = element_text(size=9, face="bold"),
        axis.title  = element_text(size=12,face="bold"))

# ggsave("06-Dec-22_Haplotype analysis results illustration3.png", 
#        last_plot(), width = 12, height = 4, pointsize = 4, dpi = 300, units = "in")
#----------#
#----------#
#----------#

#=========================================#
#          CHRIS main workflow
#=========================================#


#-----------------------------------------------------#
#----------------------- LocusZoom -------------------
#-----------------------------------------------------#

#locuszoom data preparation for Serum Creat log
mergedSerCretLog <- read.delim("mergedtable_hg19_SerumCreatLog.txt", header= TRUE)
#---------#

mergedSerCretLog %>%
  mutate(END    = BEG,
         Zscore = BETA/SEBETA,
         RH0    = ifelse(PVALUE < 0.05/256, "Yes", "No"))
#---------#

mergedeGFR <- read.delim("mergedtable_hg19.txt", header= TRUE)

mergedeGFR %>%
  mutate(Phenotype = rep("eGFR", nrow(mergedeGFR)),
         Zscore = BETA/SEBETA,
         RH0    = ifelse(PVALUE < 0.05/256, "Yes", "No"))


write.table(mergedSerCretLog, file = "mergedSerCretLog.txt",
            sep = "\t", row.names = FALSE, quote=FALSE)
#---------#

mergedTable <- merge(mergedSerCretLog[c("CHR", "BEG", "BETA", "SEBETA", "PVALUE", "Phenotype","Zscore", "RH0")],
                     mergedeGFR[c("CHR", "BEG", "BETA", "SEBETA", "PVALUE", "Phenotype","Zscore", "RH0")],
                     by = c("CHR", "BEG"), suffixes = c(".LogSerumCreat",".eGFR"))
#---------#

merged2Table <- data.frame("Z-LogSerumCreat" = mergedSerCretLog$Zscore,
                           "Z-eGFR"          = mergedeGFR$Zscore)
#---------#

#mergedTable$Zscore <- mergedTable$BETA/mergedTable$SEBETA
plot(mergedSerCretLog$Zscore, mergedeGFR$Zscore)
#---------#

ggplot(mergedTable[order(mergedTable$Zscore.eGFR),], 
       aes(Zscore.LogSerumCreat, Zscore.eGFR)) + 
  geom_point(color = "steelblue", fill = "gray") + 
  geom_abline(intercept = 0, slope = -1, color="hotpink", size=.67)
#---------#

mergedMarkerIDhg19 <- read.delim("mergedtable_hg19_SerumCreatLog_MARKER_ID37.txt", header= TRUE)
#---------#

write.table(mergedMarkerIDhg19[, c("CHR", "BEG", "END", "MARKER_ID", "RS.number")], 
            file = "mergedMarkerIDhg19.txt", sep = "\t", row.names = FALSE, quote=FALSE)

write.table(mergedMarkerIDhg19[,"RS.number"], file = "rsIDs_mergedtable_hg19_SerumCreatLog.txt", 
            sep = "\t", row.names = FALSE, quote=FALSE)



#-----------------------------------------------------#
#------------------ Table of Variables ---------------
#-----------------------------------------------------#

library(readxl)
varTable <- read_excel("D://Dariush//PhD//Analysis//Data//CHRIS Baseline_total//Query.xlsx")
#---------#

varTable %>%
  filter(source=="Labdata") %>%
  filter(module_en!="Blood collection") %>%
  filter(subquest_en=="value") %>%
  select("quest_en", "qid") %>% 
  mutate(varName =paste0(quest_en, " = CHRISbase$",qid)) %>%
  write_csv("vars.csv", "varName")
#---------#

table(varTable[varTable$source=="Labdata" & varTable$module_en!="Blood collection", "quest_en"]) %>% View()

View(varTable[(varTable$source=="Labdata" & varTable$subquest_en=="value" |
               varTable$source=="Labdata" & varTable$subquest_en=="quantitative"), 
               c("quest_en", "qid")])



#-----------------------------------------------------#
#---------------- Supplementary Table6 ---------------
#-----------------------------------------------------#

Supptable <- read.delim("E://Dariush//PhD//Analysis//Data//ST6.txt", header= TRUE, sep = "\t")
str(Supptable)
head(Supptable)

Supptable$X <- NULL
Supptable <- Supptable[Supptable$Support.for.kidney.function.relevance==1,]

#library(stringr)
#str_split_fixed(Supptable$Chr.Pos..b37., ":", 2)
#do.call(rbind, str_split(Supptable$Chr.Pos..b37., ':'))
#-----------------------------------------------------#
library(data.table) ## v 1.9.6+ 
setDT(Supptable)[, paste0("Chr.Pos..b37.", 1:2) := tstrsplit(Chr.Pos..b37., ":")]
names(Supptable) [names(Supptable)=="Chr.Pos..b37.1"] <- "CHR"
names(Supptable) [names(Supptable)=="Chr.Pos..b37.2"] <- "BEG"
Supptable$BEG <- as.numeric(Supptable$BEG)
Supptable$chr <- as.integer(Supptable$chr)

data.frame(ch1$CHR[1:20], Supptable$chr[1:20], 
           ch1$BEG[1:20], Supptable$BEG[1:20], 
           ch1$BEG[1:20]==Supptable$BEG[1:20])
#-----------------------------------------------------#
#way01
Amerg <- merge(ch1, Supptable, by.x = "BEG")#, all.x = TRUE, all.y = TRUE)
Amerg <- merge(ch1, Supptable, by = c("CHR", "BEG"), all.x = FALSE, all.y = TRUE)
#way02
library(data.table)
setDT(ch1)[as.character(BEG) %chin% as.character(Supptable$BEG)]
#way03
ch1$BEG<-as.character(ch1$BEG)
library(dplyr)
ch1 %>% filter(BEG %in% Supptable$BEG)
#way04
subset(ch1, BEG %in% Supptable$BEG)

write.table(Supptable, "Supptable.txt", sep = "\t", row.names = FALSE, quote=FALSE)
#-----------------------------------------------------#

ch2 <- ch1[1:20000,]
ch3 <- ch1[20001:50000,]
ch4 <- ch1[50001:70000,]
ch5 <- ch1[70001:90000,]
ch6 <- ch1[90001:123000,]

ch0 <- list(ch2$PVALUE, ch3$PVALUE, ch4$PVALUE, ch5$PVALUE, ch6$PVALUE)
ch00 <- list("ch7", "ch8", "ch9")
#-----------------------------------------------------#
MH_plot_PNG <- function(A, MH_file_name = "MH_plot.png", mytitle = "", color1 = "darkred", color2 = "blue", pch.type = 16)
{
  A <- A[order(A$chr, A$position),]
  
  A$log10p <- -log10(A$pval_gc)
  n.snps <- nrow(A)
  ymax <- max(A$log10p)+1
  A$abs_position <- 1:n.snps
  
  png(MH_file_name, width=2400, height=1400)
  
  par(mar=c(3,5,1,1))
  plot(A$abs_position, A$log10p, type="n", yaxt="n", xaxt="n", xlab="", ylab="", main=mytitle, xlim=c(0,n.snps), ylim=c(0,ymax))
  chr <-  1:22
  middle.points <- numeric(22)
  for (i in 1:22)
  {
    idx <- A$chr==i
    points(A$abs_position[idx], A$log10p[idx], col=ifelse(i%%2==0, color1, color2), pch=pch.type, cex=0.8)   
    middle.points[i] <- median(A$abs_position[idx])
  }
  axis(side=2, at=seq(from=0, to=ymax, by=1), labels=seq(from=0, to=ymax, by=1), tick=T, cex.axis=1, las=1)
  axis(side=1, at=middle.points, labels=as.character(c(1:22)), tick=F, cex.axis=1.5, las=1, pos=0)
  mtext(text="Chromosomes", side=1, cex=1.5, line=1.3)
  mtext(text=expression(-log[10]* (p-value)), side=2, line=3, cex=1.5)
  abline(h=-log10(5*10^-8),lty=2,col="black")
  
  dev.off()
}

lapply( list(ch2, ch3), function(ch) { ch$Avg <- rowMeans(ch[1:2]); ch })
#-----------------------------------------------------#
gc_lambda1 <- function(p)
{
  # p = p-value vector
  z <- qnorm(0.5 * p, mean=0.0, sd= 1.0, lower.tail=F)
  median(z * z, na.rm=TRUE) / qchisq(0.5, df=1, lower.tail=F)
}
#-----------------------------------------------------#
# Lambda estimation

#mylambda <- append(mylambda, gc_lambda1(U$PVALUE))
mylambda <- character(0)
mylambda <- append(mylambda, lapply(ch0, FUN = gc_lambda1))
mylambda

#-----------------------------------------------------#
setwd("Q:/ExamData/2018")
filelist <- list.files()
#-----------------------------------------------------#
library(data.table)
Exams2018 <- rbindlist(sapply(filelist, fread, simplify = FALSE, use.names = TRUE, idcol = "FileName"))

list.files(pattern = "\\.gz$")

for(i in 1:5) {
  #oname = 
  paste("file", i, sep="")
  #assign(oname, read.csv(paste(oname, ".txt", sep="")))
}


address <- "C://Users//Dariush_G3//Desktop//BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt//BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt"
read.delim(address, header= TRUE, sep = " ", nrows = 10)







#-----------------------------------------------------#
#---------------- Characterization of loci -----------
#-----------------------------------------------------#

#Less significant SNPs in each locus: top_n(n=1) 
repSNPs %>% 
  select(SNPid, Locus, BETA.log.Std.Corrected, PVALUE)%>%
  group_by(Locus) %>% 
  top_n(n = 1, wt = PVALUE)

#------------#
#Most significant SNPs
repSNPs %>% 
  select(CHR, SNPid, Locus, RSID, BETA, PVALUE)%>%
  group_by(Locus) %>% 
  slice_min(n = 1, order_by = PVALUE) %>% 
  arrange(SNPid) %>% 
  select(Locus) -> loci_names

#------------#
file_path <- "D:\\Dariush\\PhD\\Analysis\\Data\\replicatedLoci\\"
file_path %>% list.files(full.names =TRUE)

#------------#
#Extracting name of the csv files
#loci_file_names <- 
file_path %>% 
  list.files()%>% 
  .[str_detect(., ".csv")] #-> loci_file_names

#------------#
# Load everything into the Global Environment
loci_file_names %>%
  purrr::map(function(file_name){
    # iterate through each file name
    assign(x = str_remove(file_name, ".csv"),
      # Remove file extension ".csv"
           value = read_csv(paste0(file_path, file_name)),
           envir = .GlobalEnv)  })

#------------#
#1st approach
library(vroom)
file_path %>% list.files(full.names = TRUE) %>% vroom(id = "FileName")

#------------#
#2nd approach
library(data.table)
library(stringi) #stri_pad_left function

lociTrait <- 
  file_path %>% 
  list.files(full.names = TRUE) %>%
  #sprintf(., fmt = "%02s") %>% 
  #gsub(., (pattern = " ", replacement = "0")
  #str_pad(., width = nchar(.) + 1, pad = "0", side = "left") %>%
  #stri_pad_left(str = ., width = nchar(.) + 1, pad="0")
  #str_replace(pattern = "\\b([1-9]{1})_", replacement = "0\\1_") %>%
  #sort() %>% 
  set_names(., loci_names$Locus) %>%
  map_df(fread, .id = "Locus")

write.csv(lociTrait, "21-Feb-2022_Locus associations with traits.csv", row.names = FALSE)




#-----------------------------------------------------#
#----------------- word-cloud generator --------------
#-----------------------------------------------------#

library(RColorBrewer)
library(wordcloud)
library(wordcloud2)
#------------#

# Reading the saved above outputs
lociTrait <- read.csv("D:\\Dariush\\PhD\\Analysis\\Outputs_tables\\21-Feb-2022_Locus associations with traits.csv", header = TRUE)
#------------#

# Contingency table
lociTrait.freq <- lociTrait %>% 
  #group_by(Locus) %>% count(Trait)
  filter(Locus == "GAB2") %>% 
  select(Trait) %>% 
  count(Trait)

wordcloud(words = lociTrait.freq$Trait, 
          freq = lociTrait.freq$n,
          min.freq = 1,
          max.words=200, 
          random.order=FALSE, 
          rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))
#------------#

# Word cloud of CKDGen Loci
supTable <- read.table("D:\\Dariush\\PhD\\Analysis\\Outputs_tables\\Supptable.txt", header = T)
#------------#

png("08-Dec-2022_world cloud CKDGen EU-Loci.png",
    units="in", res = 500, width = 2, height = 2, bg = "transparent")

wordcloud(words = supTable$Closest.Gene,
          freq  = abs(supTable$eGFR.Effect),
          min.freq  = 1,
          max.words = Inf,      # Set top n words
          random.order = FALSE, # Words in decreasing freq
          #random.color = TRUE,
          rot.per = 0.0,        # % of vertical words
          scale = c(1,0.1),     # Set min and max scale
          use.r.layout = FALSE, # Use C++ collision detection
          fixed.asp = TRUE,
          colors = brewer.pal(8, "Dark2"))

dev.off()
#------------#

wordcloud2(data = supTable,
           #backgroundColor = "grey",
           fontWeight = supTable$CKD.OR*1000,
           size = 0.06,
           color = 'random-dark')
#------------#

library(ggwordcloud)

# 2nd way for making cloud words plot
lociTrait %>% 
  filter(Locus == "GAB2") %>%
  mutate(Rank = order(abs(Beta), decreasing = TRUE)) %>% 
  top_n(n = 20, wt = Rank) %>% #select(Trait, Beta, Rank)
  ggplot(aes(label = Trait, size = Rank, col = as.character(Rank))) + 
  geom_text_wordcloud(rm_outside = TRUE, max_steps = 1, grid_size = 1, eccentricity = .9)+
  scale_size_area(max_size = 6)+
  #scale_color_brewer(palette = "Paired", direction = -1)+
  theme_void()
#------------#

# Highly associated SNPs in the Genome-Wide significance level
lociTrait %>% 
  filter(P.value <= 5e-8) %>% 
  #select(Locus, Trait, P.value, Beta)%>% 
  count(Trait) %>% 
  write.csv(., "21-Apr-2022_Associated traits with the leading SNPs_Traits_frequency.csv", row.names = FALSE)











#-----------------------------------------------------#
#---------------- Haplotype Analysis -----------------
#-----------------------------------------------------#

library(haplo.stats)

# Export chris phenotype data for haplotype Regression analysis
vcfReg_Mag %>% 
  select(AID, Age, Sex, eGFRw.log.Res, starts_with("PC")) %>% 
  write.table("06-Nov-2022_chris4HaploReg.txt", sep = "\t", row.names = FALSE, quote=FALSE)


# Achieving hot region based on recombination rate
# for rsid 4:76480299 in SHROOM3 in GRCh38 -> n=10,758
recombSHROOM3 <- read.table("D:\\Dariush\\PhD\\Analysis\\Data\\recomb-hg38\\recomb-hg38\\plink.GRCh38.map\\plink.chr4.GRCh38.map", header = FALSE, sep = " ")

# Genotype data
vcfSHROOM3 <- 
  read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\chr4-SHROOM3_DS.txt",
             col.names = c("AID", "CHR", "POS", "MARKER_ID", "REF", "ALT", "Dosage"), 
             #nrows = 107579, 
             sep = "\t", stringsAsFactors = FALSE) #%>% count(POS)

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
#Example run
#----------#

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

#-----------#
#  SHROOM3  #
#-----------#

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

















#-----------------------------------------------------#
#---------------- Polygenic Risk Score ---------------
#-----------------------------------------------------#
library(tidyverse)
repSNPs %>% 
  select(SNPid,Locus,BETA.log.Std.Corrected)%>% 
  #group_by(Locus)%>% 
  filter(BETA.log.Std.Corrected < 0)%>% 
  head()

#-----------------------------------------------------#
#---------------------- UMOD -------------------------
#-----------------------------------------------------#
UMODvcf <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\chr16-20348509-rs13335818.txt", 
                      sep = "\t", stringsAsFactors = FALSE)

UMOD <-  merge( chris[c("AID","Sex","Age","eGFRw.log.Res","eGFR","UACR")],
                UMODvcf[c("AID", "DS")], 
                by="AID", all=FALSE)

UMOD <- UMOD %>% 
  mutate(DS_Level = cut(DS, breaks=c(-Inf, 0.500, 1.500, Inf), labels=c("0", "1", "2"))) %>% 
  mutate(Age_cat = cut(Age, breaks=c(18, 30, 40, 50, 60, 70, 94),
                       labels=c("18-30","30-40","40-50","50-60","60-70","70-94"))) #%>%
  #group_by(Age_cat) %>%
  #summarise(n(), M = mean(DS), SD = sd(DS))
ggplot(aes(Age_cat, DS)) +
geom_violin(aes(fill = Age_cat), trim = FALSE, alpha=0.9) + theme_bw()

tab <- as.matrix(prop.table(table(UMOD$DS_Level, UMOD$Age_cat), margin = 2))



#-----------------------------------------------------#
#------------------ CKDGen REGENIE --------------------
#-----------------------------------------------------#

library(RNOmni) #4 Inverse Normal Transformation
library(nephro)
#---------#
#for CKDGen REGENIE

#data4Regenie <- 
  chris %>%
  select(AID,
         Age,
         Sex,
         eGFR,
         SerumCreatinine,
         UrinaryAlbumin,
         UrinaryCreatinine) %>% 
  filter(!is.na(SerumCreatinine), 
         !is.na(UrinaryAlbumin),
         !is.na(UrinaryCreatinine)) %>% #dim()
  mutate(UACR     = (UrinaryAlbumin/10)/(UrinaryCreatinine/1000),
         #eGFRrf   = CKDEpi.creat.rf(SerumCreatinine, Sex, Age),
         eGFRrf   = CKDEpi.creat.rf(SerumCreatinine, as.numeric(as.character(Sex)), Age),
         eGFR.INT = RankNorm(eGFR),
         UACR.INT = RankNorm(UACR),
         eGFR.res.INT = RankNorm(lm(eGFR ~ Age + Sex)$residuals),
         UACR.res.INT = RankNorm(lm(UACR ~ Age + Sex)$residuals),
         FID = AID,
         IID = AID) %>%  
  #select(AID, Age, Sex, SerumCreatinine, eGFRrf) %>%
  #write.table("4Ryo.txt", sep = "\t", row.names = FALSE,quote = FALSE)
  ggplot(aes(x = eGFR)) + geom_histogram(color="Grey50", bins = 50) + theme_classic() #geom_point(alpha = 0.4)
  inner_join(PCs_raw,
             by = "AID") %>% 
  select(FID,
         IID,
         eGFR.res.INT,
         #UACR.INT,
         Age,
         Sex,
         starts_with("PC")) %>% 
    rename(#UACR = UACR.INT,
           eGFR.I_ITN = eGFR.res.INT) #%>% head()

ggsave("eGFR_9-May-22.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#---------#
#Export data for GWAS via REGENIE 
data4Regenie <- data.frame(
  FID = rep(0,nrow(vcfReg)),
  IID = vcfReg$AID,
  eGFRw.log.Res = vcfReg$eGFRw.log.Res,
  PC1 = vcfReg$PC1, PC2 = vcfReg$PC2,
  PC3 = vcfReg$PC3, PC4 = vcfReg$PC4,
  PC5 = vcfReg$PC5, PC6 = vcfReg$PC6,
  PC7 = vcfReg$PC7, PC8 = vcfReg$PC8,
  PC9 = vcfReg$PC9,PC10 = vcfReg$PC10)

write.table(data4Regenie, "pheno_I-INT.txt", sep = "\t", row.names = FALSE, quote=FALSE)
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#

write.table(head(data4Regenie),".Renviron")

x <- 1:1000
library(tidyverse)
data.frame("x" = 1:1000,
           "n" = 1:1000) %>% 
  mutate(y     = 1 - (1-0.05)^x,
         y.adj = 1 - (1-0.05 / n)^x,) %>%
  ggplot() +
  geom_point(aes(x, y.adj), alpha = 0.7, color = "springgreen3") +
  labs(x = "Number of tests",
       y = "Bonferoni-adjusted type I error")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 14, face="bold"),
        axis.title.y = element_text(size = 14, face="bold"),
        axis.text.x  = element_text(size = 10, face="bold"),
        axis.text.y  = element_text(size = 10, face="bold"))

ggsave("Type I error at 1000 times of test_Bonferoni-Adjusted.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300)  
  



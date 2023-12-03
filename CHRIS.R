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
#---------------- Polygenic Risk Score ---------------
#-----------------------------------------------------#
library(tidyverse)
repSNPs %>% 
  select(SNPid,Locus,BETA.log.Std.Corrected)%>% 
  #group_by(Locus)%>% 
  filter(BETA.log.Std.Corrected < 0)%>% 
  head()


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
  



#=========================================#
#      Visualizing biochemical traits
#=========================================#

#-----------------------------------------------------#
#-------------------- Histogram ----------------------
#-----------------------------------------------------#

library(ggplot2)

#ReformedScheme
pdf('Histogram_ReformedScheme.pdf', width=15)

lapply(c(#"SerumCreatinine.Stdw.Res",
         #"eGFRw.Res",
         #"eGFRw.log.Res",
  "SCr", "SAlb", "UAlb","UCr", "UACR"), 
       function(i)
         ggplot(chris, aes(chris[,i])) + 
         geom_histogram(aes(y=..density..), 
                        color = "steelblue", fill = "mediumturquoise", bins = 100) + 
         geom_density(alpha=0.5, fill = "steelblue")+ 
         xlab(colnames(chris[i]))+ 
         theme_minimal()#+ xlim(quantile(phenodf[,i], probs= 0.01, na.rm = TRUE), quantile(phenodf[,i], probs= 0.999, na.rm = TRUE))
)

dev.off()
#-----------------------------------------------------#

ggplot(chris, aes(HbA1c)) + stat_bin(bins=30)+
  geom_histogram(color = "black", fill = "green") 
ggplot(chris, aes(HbA1c)) + #stat_bin(bins=30)+
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=0.2, fill = "#FF8636")

ggplot(chris, aes(UrinCreatinine)) + stat_bin(bins=30)+
  geom_histogram(color = "black", fill = "gray")

ggplot(chris, aes(UrinCreatinine)) + #stat_bin(bins=30)+
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=0.2, fill = "#FF2666")
#-----------------------------------------------------#

#week
ggplot(chris, aes(week)) + 
  geom_histogram(color = "steelblue", fill = "steelblue3", bins = 100)+
  theme_minimal()

ggplot(chris, aes(week)) + 
  geom_histogram(aes(y=..density..), color = "steelblue", 
                 fill = "steelblue3", bins = 100) + 
  geom_density(alpha=0.5, fill = "steelblue")+ theme_minimal()

#eGFR
ggplot(chris, aes(eGFR)) + 
  geom_histogram(color = "steelblue", fill = "steelblue3", bins = 100)+
  theme_minimal()

ggplot(chris, aes(eGFR.Adj)) + 
  geom_density(color = "steelblue", fill = "steelblue3", alpha=0.7)+
  theme_minimal()+
  geom_vline(aes(xintercept=mean(30, na.rm=TRUE)),
             color="#FC4E07", linetype="dashed", size=1)

#log(eGFR)
ggplot(chris, aes(log(eGFR))) + 
  geom_histogram(color = "steelblue3", fill = "steelblue2", bins = 200)+
  theme_minimal()

ggplot(chris, aes(eGFR.log.Adj)) + 
  geom_density(color = "steelblue", fill = "steelblue3", alpha=0.7)+
  theme_minimal()

#Serum Creatinine
ggplot(chris, aes(SerumCreatinine)) + 
  geom_histogram(color = "orange", fill = "hotpink", bins = 200)+
  theme_minimal()

ggplot(chris, aes(SerumCreatinine)) + 
  geom_density(color = "black", fill = "orange3")+theme_minimal()
log(eGFR)
ggsave("eGFR.Adj_density.png",last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-----------------------------------------------------#

#Serum Albumin 
ggplot(chris, aes(UrinaryAlbumin)) + 
  geom_density(color = "black", fill = "skyblue")+
  theme_minimal()

ggplot(chris, aes(SerumAlbumin)) + 
  geom_density(color = "black", fill = "red")+
  theme_minimal()

#UACR
ggplot(chris, aes(UACR)) +
  geom_density(color = "#00AFBB", fill = "hotpink3", 
               bins = 150, alpha = 0.5)+xlim(0,130)

#log-log(UACR) Residuals
ggplot(chris, aes(llUACR)) + 
  geom_density(color = "black", fill = "steelblue3")+
  theme_minimal()

#Standardized
ggplot(chris, aes(UrinaryAlbumin.Std)) + 
  geom_density(color = "steelblue2", fill = "skyblue", bins = 300)+
  theme_minimal()+xlim(0,100)

ggplot(chris, aes(UrinaryCreatinine.Std)) + 
  geom_density(color = "steelblue1", fill = "skyblue", bins = 300)+
  theme_minimal()#+xlim(0,100)

ggplot(chris, aes(SerumCreatinine.Std)) + 
  ggeom_density(color = "steelblue1", fill = "skyblue", bins = 300)+
  theme_minimal()#+xlim(0,100)

ggplot(chris, aes(SerumCreatinine.Std)) +
  geom_histogram(aes(y=..density..), color="steelblue2", bins = 200) + 
  geom_density(alpha=0.2, fill = "steelblue3")+ theme_minimal()

ggsave("log(eGFR)_density.png",last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-----------------------------------------------------#

#chris$Serum2Urine.Creat
ggplot(chris, aes(Serum2Urine.Creat)) + 
  geom_histogram(color = "#E7B100", 
                 fill = "#E7B800", 
                 bins = 200, alpha = 0.5)+
  theme_minimal()

#chris$Serum2Urine.Album
ggplot(chris, aes(log(log(Serum2Urine.Album)))) + 
  geom_histogram(color = "#00AFBB", 
                 fill = "steelblue3", 
                 bins = 200, alpha = 0.5)+theme_minimal()

###geom_vline( aes(xintercept=mean(UrinCreatinine)), color="#FC4E07", linetype="dashed", size=1) + stat_bin(bins=50)

#lnUACR-GFR & 2xlnUACR-GFR
ggplot(chris, aes(UACR.GFR.ln.Std)) + 
  geom_histogram(color = "#00AFBB", 
                 fill = "steelblue3", 
                 bins = 200, alpha = 0.5)+
  theme_minimal()

ggplot(chris, aes(UACR.GFR.2xln.Std)) + 
  geom_histogram(color = "#E7B100",
                 fill = "steelblue3",
                 bins = 200, alpha = 0.5)+
  theme_minimal()

#UACR vs. UACR-GFR
UACR.join <- data.frame("id" = chris$AID,"UACR" = chris$UACR, "UACR.GFR.Std" = chris$UACR.GFR.Std)
UACR.join.melt <- melt(UACR.join, id.vars = "id")
ggplot(UACR.join.melt, aes(value)) + geom_histogram(aes(color = variable), fill="lightgray",  bins = 100, alpha = 0.7)+theme_minimal()+xlim(-20,100)

#-----------------------------------------------------#
#pdf of histogram of prepared phenotypes for GWAS
phenodf <- read.delim("phenodf.ped", header = TRUE, sep = "\t", colClasses = c("character", "character", "integer", "integer", "numeric",
                                                                               "numeric","numeric","numeric","numeric","numeric","numeric",
                                                                               "numeric","numeric","numeric","numeric","numeric",
                                                                               "numeric","numeric","numeric","numeric","numeric"), quote = "")
#-----------------------------------------------------#
pdf('Rplot45-PhenoHistBines100.pdf', width=15)
lapply(5:21, function(i)
  ggplot(chris, aes(phenodf[,i])) +
    geom_histogram(color = "steelblue2", 
                   fill = "steelblue3",
                   bins = 100, alpha = 0.7)+
    xlab(colnames(phenodf)[i])+ 
    theme_minimal()#+ xlim(quantile(phenodf[,i], probs= 0.01, na.rm = TRUE), quantile(phenodf[,i], probs= 0.999, na.rm = TRUE))
)
dev.off()

#-----------------------------------------------------#
#------------------- Scatter plot --------------------
#-----------------------------------------------------#

ggplot(chris, aes(Age, SerumCreatinine.Std))+ 
  geom_point(aes(color=Sex), alpha=0.3)+
  theme_minimal()

ggplot(chris, aes(Age, eGFR))+ 
  geom_point(aes(color=Sex), alpha=0.3)+
  theme_minimal()

ggplot(chris, aes(Age, log(eGFR)))+ 
  geom_point(aes(color=Sex), alpha=0.3)+
  theme_minimal()

#-----------------------------------------------------#
#------------------- Contour plot --------------------
#-----------------------------------------------------#

#UACR vs eGFR
library(LaplacesDemon)
joint.density.plot(chris$Age, chris$SerumCreatinine)#, na.rm==TRUE

smoothScatter(chris$UACR, chris$eGFR, bandwidth=4)

ggplot(chris, aes(UACR, eGFR)) + 
  geom_point(color="#FC2E01", alpha=0.1) + 
  theme_bw()+ xlim(0,100)

ggplot(chris, aes(UACR, eGFR)) +
  geom_point(color="steelblue2", alpha=0.2) + 
  theme_bw()+ xlim(0,100)
#ggplot(chris, aes(UACR, eGFR, z = "density")) + 
#  geom_raster(aes(fill = 'density')) + geom_contour(colour = "white")#aes(colour = after_stat(level))

ggplot(chris, aes(Age, Age, z = 'density',  scales = "free_x")) + 
  geom_contour()#aes( colour = ..level..)

ggplot(chris, aes(UACR, eGFR)) + 
  geom_bin2d(bins = 70) + 
  scale_fill_continuous(type = "viridis") +
  theme_bw() + 
  xlim(0,100) + 
  geom_vline(aes(xintercept=30), 
             color="#FC4E07", linetype="dashed", size=1) + 
  geom_hline(aes(yintercept=60), 
             color="#FC4E07", linetype="dashed", size=1)

ggplot(chris, aes(llUACR.Res, eGFR)) + 
  geom_bin2d(bins = 70) + 
  scale_fill_continuous(type = "viridis") + 
  theme_bw() + 
  geom_vline(aes(xintercept=log(log(30))), 
             color="#FC4E07", linetype="dashed", size=1) + 
  geom_hline(aes(yintercept=60), 
             color="#FC4E07", linetype="dashed", size=1)

ggsave("Rplot28.png", width = 8, height = 5.5, pointsize = 5, dpi = 1500, units = "in")

ggplot(chris, aes(llUACR.Res, eGFR)) + 
  geom_point(color = "#00AFBB", alpha=0.2) + 
  geom_density_2d(color = "#E7B800") + 
  theme_minimal()+
  geom_vline(aes(xintercept=log(log(30))), 
             color="#FC4E07", linetype="dashed", size=1) + 
  geom_hline(aes(yintercept=60), 
             color="#FC4E07", linetype="dashed", size=1)

ggsave("Rplot29.png", width = 8, height = 5.5, pointsize = 5, dpi = 1500, units = "in")

ggplot(chris, aes(llUACR.Res, eGFR)) + 
  geom_point(color = "#00AFBB", alpha=0.1) + 
  geom_density_2d(color = "#E7B800") + 
  theme_minimal()+
  geom_point() + 
  stat_density_2d( aes(fill = ..level..), geom="polygon")+
  geom_vline(aes(xintercept=log(log(30))), 
             color="#FC4E07", linetype="dashed", size=1) + 
  geom_hline(aes(yintercept=60), 
             color="#FC4E07", linetype="dashed", size=1)

ggsave("Rplot30.png", plot=Rplot30, width = 8, height = 5.5, pointsize = 5, dpi = 1500, units = "in")

#Contour plot UACR vs eGFR adjusted by AGE and SEX
#Adjusting UACR by AGE and SEX

ggplot(chris, aes(UACR.Std, eGFR)) + 
  geom_point(color = "#00AFBB", alpha=0.3) + 
  geom_density_2d(color = "#E7B800") + 
  theme_minimal()+xlim(-15,65)+
  geom_vline(aes(xintercept=30), 
             color="#FC4E07", 
             linetype="dashed", 
             size=1) + 
  geom_hline(aes(yintercept=60), 
             color="#FC4E07", 
             linetype="dashed", 
             size=1)

#UACR.2xlog.Std
ggplot(chris, aes(UACR.2xlog.Std, eGFR)) + 
  geom_point(color = "#FC2E07", alpha=0.2) + 
  geom_density_2d(color = "#E7B800") + 
  theme_minimal()


library("ggExtra")  #install.packages("ggExtra")
require("gridExtra")

# Scatter plot of x and y + Marginal density plot of x and y
scatterPlot <- ggplot(chris, aes(mu+llUACR.Res, eGFR)) + 
  geom_point(color = "#00AFBB", alpha=0.2) + 
  geom_density_2d(color = "#E7B800") +
  geom_vline(aes(xintercept = log(log(30))), 
             color="#FC4E07", 
             linetype="dashed", 
             size=1) + 
  geom_hline(aes(yintercept=60), 
             color="#FC4E07", 
             linetype="dashed", 
             size=1)+ 
  xlab("log(log(UACR.residuals))")+
  theme(legend.position=c(0,1), 
        legend.justification=c(0, 1))

xdensity <- ggplot(chris, aes(llUACR.Res)) + 
  geom_density(aes(fill=llUACR.Res), alpha=0.8) + 
  theme(legend.position = "none")

ydensity <- ggplot(chris, aes(eGFR)) + 
  geom_density(alpha=0.8) + 
  theme(legend.position = "none") + 
  coord_flip()

blankPlot <- ggplot() + 
  geom_blank(aes(1, 1)) + 
  theme_void()

grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

#mean+UACR.Residuals
mu <- mean(chris$llUACR, na.rm=TRUE) show.legend = T

ggplot(chris, aes(mu+llUACR.Res, eGFR)) + 
  geom_point(color = "#00AFBB", alpha=0.2) + 
  geom_density_2d(color = "#E7B800") +
  geom_vline(aes(xintercept=log(log(30))), 
             color="#FC4E07", 
             linetype="dashed", 
             size=1) + 
  geom_hline(aes(yintercept=60), 
             color="#FC4E07", 
             linetype="dashed", 
             size=1)+ 
  xlab("log(log(UACR.residuals))") +
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave("Rplot34.png", plot=last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")


#ggplot(chris, aes(UACR, eGFR)) + geom_hex(bins = 70) + 
#  scale_fill_continuous(type = "viridis") + theme_bw()

#Area + Contour
#ggplot(chris, aes(UACR, eGFR)) + stat_density_2d(aes(fill = ..level..), geom="raster")

#eGFR <- na.omit(chris$eGFR)
#remove(eGFR)

#-----------------------------------------------------#
#------------------- Violin plot ---------------------
#-----------------------------------------------------#

library(ggplot2)

colnames(chris) <- columnsN
colnames(chris) <- c(columnsNs,"UAlb.Sensor")

#making Instruments as Factor
chris$Sex<-as.factor(chris$Sex)
chris$HbA1c.I<-as.factor(chris$HbA1c.I)
chris$UCre.I<-as.factor(chris$UCre.I)
chris$UAlb.I<-as.factor(chris$UAlb.I)
chris$SCre.I<-as.factor(chris$SCre.I)
chris$SAlb.I<-as.factor(chris$SAlb.I)

#HbA1c
ggplot(chris, aes(HbA1c.I, HbA1c)) + geom_violin(aes(fill = HbA1c.I))
ggplot(chris, aes(HbA1c.I, HbA1c)) + geom_boxplot(aes(color = HbA1c.I))

#Urine Creatinine
ggplot(chris, aes(UCre.I, UCre)) + geom_violin(aes(fill = UCre.I))
ggplot(chris, aes(UCre.I, UCre)) + geom_boxplot(aes(color = UCre.I))

#Urine Albumin
ggplot(chris, aes(UAlb.I, UAlb)) + geom_violin(aes(fill = UAlb.I))
ggplot(chris, aes(UAlb.I, UAlb)) + geom_boxplot(aes(color = UAlb.I))

#Serum Creatinine
ggplot(chris, aes(SerumCreatinine.Ins, SerumCreatinine)) + geom_violin(aes(fill = SerumCreatinine.Ins)) + geom_boxplot(width = 0.2)+xlab("Instrument")+ylim(0,2.5)+theme_minimal()
ggplot(chris, aes(SCre.I, SCre)) + geom_boxplot(aes(color = SCre.I))

#Serum Albumin
ggplot(chris, aes(SAlb.I, SAlb)) + geom_violin(aes(fill = SAlb.I))
ggplot(chris, aes(SAlb.I, SAlb)) + geom_boxplot(aes(color = SAlb.I))

#UACR
ggplot(chris, aes(UACR)) + geom_dotplot(filmotl = "steelblue") + theme_bw() #theme_minimal() 

ggplot(chris, aes(UACR)) + geom_histogram(aes(y=..density..), colour="orange") + 
  geom_density(alpha=0.5,fill = "steelblue") + theme_bw() +xlim(0,350)

ggplot(chris, aes(Sex, UACR)) + geom_violin(aes(fill = Sex), color = "steelblue") + 
  geom_boxplot(width = 0.18) + ylim(0,18)+ theme_minimal() + scale_x_discrete(labels=c("Male", "Female"))

#Homocyst
ggplot(chrisMessy, aes(Homocyst.Ins, as.numeric(Homocyst))) + 
  geom_violin(aes(fill = Homocyst.Ins))

#-----------------------------------------------------#
#----------------------- Heatmap ---------------------
#-----------------------------------------------------#

#heatmap
HMexclud<-c("AID", "Sex", "HbA1c.Ins", "HbA1c.Unit",
            "SerumCreatinine.Ins", "SerumCreatinine.Unit",
            "SerumAlbumin.Ins", "SerumAlbumin.Unit",
            "UrinaryCreatinine.Ins", "UrinaryCreatinine.Unit",
            "UrinaryAlbumin.Ins", "UrinaryAlbumin.Unit",
            "UACR.Ins", "Instrument", "UACR.Ins2", "eGFR.cat", "UACR.cat","SCrInstrument",
            "Diabetes", "Diabetes.Type", "Diabetes.Treat", "ethnicity")

#correlation matrix
chris.correlations <- cor(chris[, !colnames(chris) %in% HMexclud], method = "pearson", use = "complete.obs")
#-----------------------------------------------------#
library(pheatmap)
png("Rplot39-6.png", units="in", res = 300, width=6, height=4)
pheatmap(chris.correlations, cluster_cols = FALSE, fontsize = 6, angle_col = "45")
dev.off()
#-----------------------------------------------------#
png("Rplot41-2.png", units="in", res = 300, width=6, height=4)
library(corrplot)
corrplot(chris.correlations, type = "upper", order = "hclust", tl.col = "black", tl.srt = 65, tl.cex = 0.5)
dev.off()

#-----------------------------------------------------#
#----------------- Leiden Presentation ----------------
#-----------------------------------------------------#

#Serum Creatinine - Violin plot
chris$SCrInstrument <- ifelse(chris$SerumCreatinine.Ins==0, paste("ROCHE"), paste("ABBOTT"))

ggplot(chris, aes(SCrInstrument, SerumCreatinine)) + 
  geom_violin(aes(fill = SCrInstrument)) + 
  geom_boxplot(width = 0.2) + 
  ylim(0,2.5) + 
  theme(panel.background = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black")) %>% 
  ggsave(last_plot(), "RplotSCr35.png", width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#instrument adjusted Serum Creatinine -> SCr.Ins.Std
ggplot(chris, aes(SCrInstrument, SerumCreatinine.Std)) + 
  geom_violin(aes(fill = SCrInstrument)) + 
  geom_boxplot(width = 0.2) + ylim(0,2.5) + 
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) %>% 
  ggsave("RplotSCrStd38.png", width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#instrument, age, sex, week adjusted Serum Creatinine -> SCr.Ins.Week.Res
ggplot(chris, aes(SCrInstrument, mean(SerumCreatinine, na.rm=TRUE) + SCr.Ins.Week.Res)) + 
  geom_violin(aes(fill = SCrInstrument)) + 
  geom_boxplot(width = 0.2) + ylim(0,2.5) + 
  ylab("SCr adj Ins&Age-Sex-Week")+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) %>%
  ggsave("RplotSCrInsWeekStd39.png", width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#APTT or eGFR
ggplot(chris, aes(APTT_ratio.Ins, APTT_ratio)) + 
  geom_violin(aes(fill = APTT_ratio.Ins)) + 
  geom_boxplot(width = 0.2) + 
  theme(panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) #%>% 
  #ggsave("RploteGFR36.png", width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#UACR
chris$UACR.Ins2  <- ifelse(chris$UACR.Ins==0, paste("Machine0-ROCHE"), (ifelse(chris$UACR.Ins==1, paste("SIEMENS-ROCHE"), paste("Machine0-ABBOTT"))))
chris$Instrument <- ifelse(chris$UACR.Ins==0, paste("Machine0-ROCHE"), (ifelse(chris$UACR.Ins==1, paste("SIEMENS-ROCHE"), paste("Machine0-ABBOTT"))))

ggplot(chris, aes(Instrument, UACR)) + 
  geom_violin(aes(fill = Instrument)) + 
  geom_boxplot(width = 0.2) + ylim(0,50) + 
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) %>%
  ggsave("RplotUACR37.png", width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#-----------------------------------------------------#
#---------------------- MH plot ----------------------
#-----------------------------------------------------#

ch1<-read.delim("F://[] PhD Courses//GWAS course//myoutput//UACR.chr1.1.180000001.190000000.epacts", header=TRUE)

str(ch1)
summary(ch1)

#install.packages("qqman")
library(qqman)

names(ch1) [names(ch1)=="X.CHROM"] <- "CHR"
names(ch1) [names(ch1)=="PVALUE"] <- "P"
ch1$BP <- seq(1:nrow(ch1))
ch1$CHR<- ifelse(ch1$CHR=="chr1", 1, 0)

manhattan(ch1, chr = "CHR", bp = "BP", p = "P", snp = "BEG", ylim = c(0, 80), cex = 0.6, 
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "P", "Q"))
qq(ch1$P)

#ggplot(ch1, aes(sample = P)) + stat_qq() + stat_qq_line()

#--- MH for merged output ---

A$BP <- seq(1:nrow(A))
names(A) [names(A)=="X.CHROM"] <- "CHR"
names(A) [names(A)=="PVALUE"] <- "P"
A$CHR <- as.numeric(str_replace(A$CHR, "chr", ""))
A <- A[order(A$CHR),]

manhattan(A, chr = "CHR", bp = "BP", p = "P", snp = "BEG", ylim = c(0, 80), cex = 0.6, 
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:20, "P", "Q"))

#-----------------------------------------------------#
#---------------- Cristian's MH plot -----------------
#-----------------------------------------------------#

A <- read.delim("C://Users//Dariush_G3//Desktop//outclean2.txt", sep=" ", header = TRUE)

names(A) [names(A)=="X.CHROM"] <- "chr"
names(A) [names(A)=="BEG"]     <- "position"
names(A) [names(A)=="PVALUE"]  <- "pval_gc"

#sub("$chr", "", A$chr)
#chartr("","chr", x = strings)
#-----------------------------------------------------#
library(stringr)
A$chr <- str_replace(A$chr, "chr", "")
#-----------------------------------------------------#
# sort by chr
A <- A[order(A$chr),]
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


MH_plot_PNG(A, "MH_plot_eGFRcrea_overall_MAFget001_2GC_b36_MAFget005_Nget50.png")
#---------------------------------#

#install.packages("data.table")
library(data.table)

setwd("/home/dghasemisemeskandeh/projects/gwas/Output/UACR.log-log.Res/out_epacts/")
filelist <- list.files(pattern = ".*.epacts")
chrisRes <- rbindlist(sapply(filelist, fread, simplify = FALSE, use.names = TRUE, idcol = "FileName"))
#---------------------------------#

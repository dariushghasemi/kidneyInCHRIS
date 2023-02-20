#=========================================#
#          Sensitivity Analysis
#=========================================#

# ---- Sensitivity to Municipality ----

withMunicipality <-
  step3_to_table(vcfReg[c("FT3", "FT4")],
                 vcfReg[targets],
                 paste("log(eGFR) ~ SNP + trait + Sex + Age + Municipality"),
                 vcfReg)

withoutMunicipality <-
  step3_to_table(vcfReg[c("FT3", "FT4")],
                 vcfReg[targets],
                 paste("log(eGFR) ~ SNP + trait + Sex + Age"),
                 vcfReg)

# Merging two models results

withMunicipality %>% 
  inner_join(withoutMunicipality,
           by = c("SNPid", "pheno"),
           suffix = c("_withMunic", "_withoutMunic")) %>%
  mutate(Beta_to_SE_withMunic    = Estimate_withMunic/SE_withMunic,
         Beta_to_SE_withoutMunic = Estimate_withoutMunic/SE_withoutMunic) %>% 
  inner_join(repSNPs[c("SNPid", "Locus")],
             by = c("SNPid")) %>%
  # pivot_longer(cols = -c(SNPid, pheno),
  #              names_to = c(".value", "covariate"),
  #              names_pattern = "(.+).(withMunic|withoutMunic)$", #"(Estimate|SE|Pvalue)(\\d+)",
  #              values_to = "value") %>%
  select(SNPid, Locus, pheno, everything()) %>% #View()
  #write.csv(., "14-Dec-2022_Beta_to_SE of the SNPs with and without adjusting for municipality.csv", row.names = F, quote = F)
  ggplot(aes(x = Beta_to_SE_withoutMunic,
             y = Beta_to_SE_withMunic,
             shape = Locus,
             color = Locus)) +
  geom_abline(slope = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 0.9, size = 2.5) +
  scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  scale_color_manual(values= c("maroon1", "darkorchid2", "orange2", "green4", "steelblue2", "darkturquoise",
                                        "tomato", "springgreen2", "royalblue2", "gold", "grey50")) +
  facet_wrap(~pheno, nrow = 1) +
  labs(x = "Beta to SE of the SNPs",
       y = "Beta to SE of the SNPs adjusted for municipality")+
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.placement = "outside",
        axis.text.x  = element_text(size = 8,  face = "bold"),
        axis.text.y  = element_text(size = 7,  face = "bold"),
        axis.title   = element_text(size = 14, face = "bold"),
        legend.key.size  = unit(0.99, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")) #, #axis.title.x.bottom = element_blank(),

ggsave("17-Jan-23_Beta to SE of the SNPs with and without adjusting for municipality.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")


#-----------------------------------------------------#
#----------- Sensitivity to missing FT3/4 ------------
#-----------------------------------------------------#
chris %>%  group_by(TSH.Ins) %>% median(TSH, na.rm = TRUE)

PWS_SAMadj <- 
  lapply(targets, 
         function(mytarget){
           map2_df("eGFRw.log.Res", v1, 
                   function(y,x){
                     vcfRegX <- vcfReg %>% filter(!is.na(vcfReg[x]))
                     trait   <- vcfRegX[,x]
                     outcome <- vcfRegX[,y]
                     SNP     <- vcfRegX[,mytarget]
                     m <- lm(outcome ~ trait + SNP + Operator + Municipality, data = vcfRegX) #trait + SNP + Operator + Municipality
                     return(coef(m)[3])
                   })}#%>% unlist()
  ) #%>% do.call(rbind,.) or bind_rows(.id = "SNP") #

#---------#
Messy2Tidy <- function(df, traits = v1)
{
  n   <- length(df)
  dft <- lapply(1:n, function(i) t(df[[i]]))
  dff <- as.data.frame(do.call(rbind, dft))
  colnames(dff) <- traits
  rownames(dff) <- 1:n
  out <- cbind(SNPid = targets[1:n], dff)
  return(out)
}
#---------#
SNPinMissing    <- Messy2Tidy(PWS_SA)
SNPinNonMissing <- Messy2Tidy(PWS_SAM)
SNPinNonMissingAdj <- Messy2Tidy(PWS_SAMadj)

SensAna <- cbind(repSNPs[, c("SNPid","BETA")],
                 SNPinNonMissing[-1],
                 SNPinNonMissingAdj[-1],
                 SNPinMissing[-1])

write.csv(cbind(repSNPs["Locus"], SensAna) , "SensitivityAnalysis_v1+CHIRIS_GWAS.csv")

#---------#
#png("Rplot39-6.png", units="in", res = 300, width=6, height=4)

png("Heatmap162.png", units="in", res = 300, width=10, height=12)
#pdf('pheatmap2.pdf', width=18, height = 18)
pheatmap(PWS3[c(1:70,72:124,126:164),-c(1)] / PWS3[c(1:70,72:124,126:164),2], #1:70, 72:124, 126:164
         cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = T, 
         labels_row = PWS3$SNPid, border_color = NA, fontsize = 6, angle_col = "45")

pdf('pheatmap_SensitivityAnalysis_M1-M3.pdf', width=18, height = 18)
pheatmap(SNPinNonMissing[,-1] - SNPinMissing[,-1], 
         cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = T, 
         labels_row = SNPinMissing$SNPid, border_color = NA, fontsize = 6, angle_col = "45")

pdf('pheatmap_SensitivityAnalysis_M2onM1.pdf', width=18, height = 18)
pheatmap(SNPinNonMissingAdj[,-1] / SNPinNonMissing[,-1], 
         cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = T, 
         labels_row = SNPinNonMissingAdj$SNPid, border_color = NA, fontsize = 6, angle_col = "45")

pdf('pheatmap_SensitivityAnalysis_M1-3.pdf', width=18, height = 18)
pheatmap(SensAna[,-1] / SensAna[,2], 
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = T, 
         labels_row = SensAna$SNPid, border_color = NA, fontsize = 6, angle_col = "45")

dev.off()


#-----------------------------------------------------#
#---------- Wilcoxon testing missing FT3/4 -----------
#-----------------------------------------------------#

#Wilcoxon test for Dosage difference in missed T3 or non-missed T3

wilcox2table <- function(SNP, FT)
  {
  test <-
    vcfReg %>%
    mutate(FT_indicator = if_else(is.na(!!sym(FT)),
                                  "missed",
                                  "nonmissed")) %>%
    wilcox.test(vcfReg[,SNP] ~ FT_indicator, data = .)
  #result <- rbind(test$statistic, test$p.value)  #%>% broom::tidy()
  result <- data.frame("Wstatistic" = test$statistic,
                       "Pvalue"     = test$p.value,
                       trait = FT)
  return(result)
  }

sapply(targets, wilcox2table, "FT4") %>% 
  as.data.frame %>% 
  rownames_to_column(var = "variable") %>%
  pivot_longer(cols      = -variable,
               names_to  = "SNPid",
               values_to = "stat") %>% 
  pivot_wider(id_cols     = SNPid,
              names_from  = variable,
              values_from = stat) %>%
  left_join(repSNPs[c("SNPid", "Locus")], ., by = "SNPid") %>% #View()
  data.table::fwrite("17-Jan-23_Wilcoxon test sensing dosage of SNPs to having missing in FT4.csv", row.names = F, quote = F)



#-----------------------------------------------------#
#------------- Violin plot missing FT3/4 -------------
#-----------------------------------------------------#

#Changing the SNPid lables in facet_grid
SNPid_factor <- function(x) {
  factor(x,
         levels = c("chr1:10670853","chr2:15642347",
                    "chr4:76480299","chr5:39393631",
                    "chr5:177386403","chr7:77714744",
                    "chr8:23885208","chr9:68817258",
                    "chr11:78324786","chr15:98733292", "chr16:20381010"),
         labels = c("CASZ1\nchr1:10670853",
                    "DDX1\nchr2:15642347",
                    "SHROOM3\nchr4:76480299",
                    "DAB2\nchr5:39393631",
                    "SLC34A1\nchr5:177386403",
                    "TMEM60\nchr7:77714744",
                    "STC1\nchr8:23885208",
                    "PIP5K1B\nchr9:68817258",
                    "GAB2\nchr11:78324786",
                    "IGF1R\nchr15:98733292",
                    "PDILT\nchr16:20381010"))
}

#-----------#
# Checking the sensitivity of Dosage to having missing in T3/4 with violin plot

missPlot <- function(FT)
{
  #FT <- enquo(FT)
  vcfReg %>%
    select(tagSNPs$SNPid, all_of(FT), eGFRw.log.Res) %>% 
    pivot_longer(cols = -c(eGFRw.log.Res, !!FT),
               names_to  = "SNPid",
               values_to = "Dosage") %>%
    mutate(fT = ifelse(!is.na(!!sym(FT)), "Non-missed", "Missed"),
           Dosage_Level = cut(Dosage,
                              breaks = c(-Inf, 0.500, 1.500, Inf),
                              labels = c("0", "1", "2"))) %>%
    ggplot(aes(Dosage_Level, eGFRw.log.Res)) +
    geom_violin(aes(fill = fT), trim = F, alpha=0.9, color = "white", position = position_dodge(.9)) +
    geom_boxplot(aes(fill = fT), color = "grey", width = 0.1, position = position_dodge(0.9)) +
    scale_fill_manual(values = c("steelblue2", "royalblue4")) +
    facet_wrap(~ SNPid_factor(SNPid) , nrow = 2, shrink = T) +
    #ylim(-.5, .5) +
    labs(fill = FT,
         x = "Dosage level of the variant",
         y = "ln(eGFR)")+
    theme_classic() +
    theme(strip.background = element_blank(), strip.placement = "outside",
          strip.text.x = element_text(size = 8, color = "Black", face = "bold.italic"),
          legend.position = c(.93, .265))
    
  # Save the plot
  ggsave(paste0("13-Jan-23_Association of SNPs with ln(eGFR) respect to having missing in ", FT, ".png"),
         last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
}

lapply(c("FT3", "FT4"), missPlot)

#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#

# Histogram
vcfReg %>% 
  #select(p) %>% 
  #drop_na(p) %>% 
  #ggplot(aes(BirthWeight1)) +
  #ggplot(aes(Pregnancy1)) +
  #ggplot(aes(PTT_sec)) +
  #ggplot(aes(T3)) +
  ggplot(aes(T4ng.dL)) + 
  geom_histogram(aes(y=..density..), colour="gray", fill="hotpink2", bins = 100) + 
  geom_density(alpha=0.5, colour="gray", fill = "#FF2666") + 
  theme_minimal()

#-----------#
vcfReg %>%
  mutate(Pregnant = if_else(is.na(Pregnancy1), "No", "Yes")) %>%
  ggplot(aes(`chr1:10599281`)) +
  geom_histogram(aes(fill = Pregnant), colour="black", bins = 100) +
  theme_minimal()
#geom_density(aes(fill = Pregnant), alpha=0.5, fill = "#FF2666") +
#geom_violin(aes(fill = Pregnant)) + geom_boxplot(width = 0.01) + # for violin plot
#geom_dotplot(aes(fill = Pregnant), binaxis='y', stackdir='center') + # for violin plot

#-----------#
library(tidyverse)
library(ggplot2)
#-----------#

pdf('density of potential signals 5 SNPs.pdf', width=15)
#-----------#
# Histogram
lapply(c("BirthWeight1", "Pregnancy1", "PTT_sec", "T3", "T4"), function(p)
  vcfReg %>% 
    select(p) %>% 
    #drop_na(p) %>% 
    ggplot(aes_string(p)) + 
    geom_histogram(aes(y=..density..), colour="gray", fill="hotpink2", bins = 100) + 
    geom_density(alpha=0.5, colour="gray", fill = "#FF2666") +
    theme_minimal())

#-----------#
# Density plot
lapply(c("BirthWeight1", "Pregnancy1", "PTT_sec", "T3", "T4"), function(p)
  vcfReg %>% 
    select(eGFRw.log.Res, p) %>% 
    mutate(Missed = if_else(!is.na(vcfReg[p]), "No", "Yes")) %>% 
    ggplot(aes(eGFRw.log.Res))+ geom_density(aes(fill = Missed), color = "steelblue", alpha=0.4)+ 
    labs(fill = paste("Missed", p)) + 
    theme(panel.background = element_blank(), 
          legend.position="top"))# + theme_minimal())

#-----------#
# scatter
lapply(c("BirthWeight1", "Pregnancy1", "PTT_sec", "T3", "T4"), function(p)
  vcfReg %>% 
    ggplot(aes_string(p, "eGFRw.log.Res")) + 
    geom_point(color = "steelblue", fill = "gold1", alpha=0.4, scale = "free_x")+ theme_minimal())

#-----------#
# Combine with violin plot
lapply(c("BirthWeight1", "Pregnancy1", "PTT_sec", "T3", "T4"), function(p)
  vcfReg %>%
    mutate(Missed = if_else(!is.na(vcfReg[p]), "No", "Yes")) %>% 
    ggplot(aes(`chr15:98708127`, eGFRw.log.Res)) +
    geom_point(aes(color = Missed), alpha=0.4, size=2) + theme_minimal()+
    labs(color = paste("Missed", p)) + theme(legend.position="top"))

dev.off()

#-----------------------------------------------------#
#--------- violin plot: Traits vs Dosage level --------
#-----------------------------------------------------#

v1 <- c("BirthWeight1", "Pregnancy1", "PTT_sec", "T3", "T4")#"TSH",
v2 <- c("chr1:10599281", "chr2:15642347", "chr4:76444011", "chr5:39385539", "chr5:177386403", 
        "chr7:77714744", "chr8:23885208", "chr9:68782831", "chr11:78315470", "chr15:98708127")
v3 <- c("CASZ1","DDX1", "SHROOM3","DAB2", "SLC34A1","RSBN1L", "STC1","PIP5K1B","GAB2", "IGF1R")
#-----------#
pdf('Traits vs Dosage level for 10 SNPs.pdf', width=15)
#-----------#
#comparing violin for categorical dosage level
map2(rep(c("BirthWeight1", "Pregnancy1", "PTT_sec", "TSH", "T3", "T4"),10),
     rep(v2, each=6),
     function(p, SNP) 
       vcfReg %>%
       mutate(Dosage_Level = cut(vcfReg[,SNP], 
                                 breaks=c(-Inf, 0.500, 1.500, Inf), 
                                 labels=c("0", "1", "2"))) %>% 
       select(p, SNP, quote(Dosage_Level)) %>% 
       ggplot(aes_string(quote(Dosage_Level), p)) +
       geom_violin(aes(fill = Dosage_Level), trim = FALSE, alpha=0.9) + 
       scale_fill_manual(values=c("paleturquoise2", "wheat2", "hotpink3"))+
       geom_boxplot( width = 0.1, position=position_dodge(0.9)) + 
       #geom_dotplot(aes(fill = Pregnant), binaxis='y', stackdir='center', tackratio = .1, dotsize = .6, position=position_dodge(0.8)) +
       theme_bw() + xlab(paste("Dosage of", SNP)) + theme(legend.position="top"))  

dev.off()
#-----------#

library(recipes)#recipe(~.) %>% 

vcfReg %>%
  mutate(Dosage_Level = cut(vcfReg[, "chr1:10599281"], 
                            breaks=c(-Inf, 0.500, 1.500, Inf), 
                            labels=c("0", "1", "2"))) %>%
  group_by(Dosage_Level) %>% 
  #miss_var_run(var = BirthWeight1) #naniar:::miss_var_table()
  summarise(n= n(), #m = across(c(BirthWeight1, Pregnancy1, PTT_sec, TSH, T3, T4ng.dL), ~ mean(.x, na.rm = TRUE)),
            missing = across(c(BirthWeight1, Pregnancy1, PTT_sec, TSH, T3, T4ng.dL), ~sum(is.na(.x))/n)) %>%
  #pivot_longer(BirthWeight1, names_to = "missing")
  ggplot(aes(missing$T3)) + geom_bar(aes(fill = Dosage_Level)) + coord_flip()
#-----------#
vcfReg %>% 
  mutate(Dosage_Level = cut(vcfReg[, "chr1:10599281"], 
                            breaks=c(-Inf, 0.500, 1.500, Inf), 
                            labels=c("0", "1", "2"))) %>%
  select(all_of(v1), Dosage_Level) %>% #head()
  gg_miss_var(show_pct = TRUE, facet = Dosage_Level)

#-----------#
map_dbl(.x = vcfReg[10:20], .f = mean)
#-----------#

vcfReg %>%  
  mutate(Dosage_Level = cut(`chr1:10599281`, 
                            breaks=c(-Inf, 0.500, 1.500, Inf), 
                            labels=c("0", "1", "2"))) %>% 
  select(Dosage_Level, `chr1:10599281`, Pregnancy1) %>% 
  group_by(Dosage_Level) %>% 
  count(Pregnancy1) %>% tail()


#ggsave("Pregnancy vs eGFRw-log-Res.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-----------#
library(corrr)
vcfReg %>% select(eGFRw.log.Res, TSH) %>% correlate()

vcfReg %>%
  mutate(Dosage_Level = cut(vcfReg[,"chr2:15642347"], 
                            breaks=c(-Inf, 0.500, 1.500, Inf), 
                            labels=c("0", "1", "2"))) %>% 
  select(BirthWeight1, `chr2:15642347`, Dosage_Level) %>% #head()
  ggplot(aes(Dosage_Level, BirthWeight1)) +
  geom_violin(aes(fill = Dosage_Level), trim = FALSE, alpha=0.9) + 
  scale_fill_manual(values=c("paleturquoise2", "wheat2", "hotpink3"))+
  geom_boxplot( width = 0.1, position=position_dodge(0.9)) + theme_minimal()



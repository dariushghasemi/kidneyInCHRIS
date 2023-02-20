#=========================================#
#       Interaction with Thyroid
#=========================================#

library(interactions)
library(tidyverse)
library(reshape2)
library(tidyr)
library(knitr) #for printing html-friendly tables.

#-----------------------------------------------------#
#-------------- Preparing data for interaction -------
#-----------------------------------------------------#

#Excluding participants with serious 
#Thyroid problems for interaction model

vcfReg_TSHmod <-
  vcfReg %>%
  mutate(Thyroid_DrugName = na_if(Thyroid_DrugName, "")) %>%
  filter(
    !if_all(c(Thyroid_DrugName, TSH), is.na),                        # Removed are:   4 with NAs
    #Cancer                 != 1 | is.na(Cancer),                    # Removed are: 356 with Cancer
    Thyroid_cancer          != 1 | is.na(Thyroid_cancer),            # Removed are:  16 with Thyroid_cancer (only 3 more removed after removing Cancer cases)
    kidneyCancer            != 1 | is.na(kidneyCancer),              # Removed are:   1 with Kidney cancer
    Goiter                  != 1 | is.na(Goiter),                    # Removed are: 277 with Goitre
    Operation_thyroid_gland != 1 | is.na(Operation_thyroid_gland),   # Removed are: 312 with operation on thyroid gland
    #Alteration_thyroid_pregnancy != 1 | is.na(Alteration_thyroid_pregnancy), # Removed are: 64 removed (but 62 should be eliminated)
    ) %>%
  mutate(TSH_cat = replace(TSH_cat, Thyroid_DrugName == "Iodine therapy",       "HyperT"),
         TSH_cat = replace(TSH_cat, Thyroid_DrugName == "Levothyroxine sodium", "HypoT"),
         TSH_cat = replace(TSH_cat, Thyroid_DrugName == "Propylthiouracil",     "HyperT"),
         TSH_cat = replace(TSH_cat, Thyroid_DrugName == "Thiamazole",           "HyperT"),
         #Changing the reference level of TSH_cat to TSH_cat = "2" or "NormT"
         TSH_cat = relevel(as.factor(TSH_cat), ref = 2)) 
#count(TSH_cat)
# group_by(TSH_cat, Thyroid_DrugName) %>%
# summarise(n = n()) %>%
# spread(TSH_cat, n)



#-----------------------------------------------------#
#-------------- Emerging baches effect  --------------
#-----------------------------------------------------#

# TSH vs. TSH Instruments
chris %>% 
  select(TSH, TSH.Ins) %>% 
  pivot_wider(names_from = TSH.Ins, values_from = TSH, names_prefix = "TSH_") %>% 
  head()

#-------------#
# Some randome ideas for interaction
#-------------#

anova(lm(eGFRw.log.Res ~ TSH_cat * `chr5:177386403` + Municipality, data = vcfReg),
      lm(eGFRw.log.Res ~ TSH_cat + `chr5:177386403` + Municipality, data = vcfReg))

table(vcfReg$TSH_cat, vcfReg$Thyroid_DrugName)

lapply("TSH", function(x) quantile(chris[x], probs=seq(0,1,0.1), na.rm=TRUE))

vcfReg[vcfReg$Cancer == 2 | is.na(vcfReg$Cancer), "TSH_cat"]
vcfReg[vcfReg$Cancer != 1,]

#-----------------------------------------------------#
#      Visualizing TSH for different instruments
#-----------------------------------------------------#

chris %>%
  #mutate(TSH.Std = minMaxNorm(chris$TSH)) %>% 
  #group_by(TSH.Ins) %>% 
  #summarise(m1 = mean(TSH.Std, na.rm = TRUE), sd1 = sd(TSH.Std, na.rm = TRUE),
  #          m2 = mean(TSH, na.rm = TRUE),     sd2 = sd(TSH, na.rm = TRUE), n = n())
  ggplot(aes(x= TSH.Ins, y= TSH.q)) + #geom_histogram(bins = 150) + xlim(0, 0.15)
  geom_violin(aes(fill = TSH.Ins), alpha = 0.3, trim = FALSE) + #ylim(0, 8) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  theme_classic() #%>%
#ggsave(., filename = "TSH.q vs Instrument.png", width = 8, height = 5.5, dpi = 300)#, pointsize = 5

#-------------#
chris %>% 
  ggplot(aes(x= T3)) + 
  geom_histogram(color = "steelblue1", fill = "Gold2", bins = 70) +
  #geom_point(color = "steelblue1") +
  #xlim(0, 12) + ylim(0, 8) + 
  theme_classic()

#-------------#
# I'd be tempted to do 
ggsave2 <- function(plot, filename, ...) {
  ggplot2::ggsave(filename, plot, ...)
  invisible(plot)
}

# So you could then do 
plot %>% ggsave2("x.pdf") %>% ggsave2("x.png") %>% ggsave2("x.eps")


#-----------------------------------------------------#
#-------------- Thyroid Questionnaire ----------------
#-----------------------------------------------------#

# Contingency table
crossing_tables <- function(x){
  CHRISbase %>% 
    #mutate(TSH = replace(x0lp35, x0lp35 == "-89", NA)) %>%
    #mutate(TSH_cat = cut(as.numeric(TSH), breaks=c(-Inf, 0.401, 3.799, Inf), labels=c("0", "1", "2"))) %>%
    group_by(x0dd24, x) %>% #TSH_cat
    summarise(n = n()) %>%
    #mutate(p = n / sum(n))
    spread(x0dd24, n) %>% #TSH_cat # x0th01, x0dd24
    knitr::kable()
}

map(c("x0th00", "x0th01"), crossing_tables)



#-----------------------------------------------------#
#------------------- Kidney Questionnaire ------------
#-----------------------------------------------------#

phenodf$AID <- phenodf$FAM_ID
Kidney$AID  <- as.character(Kidney$AID)
chrisKidQue <- merge(Kidney, phenodf, by = "AID")
#---------#

chrisKidQue[] <- lapply(chrisKidQue, function(x) if(is.factor(x)) as.character(x) else x)


#xlim(quantile(phenodf[,i], probs= 0.01, na.rm = TRUE), quantile(phenodf[,i], probs= 0.999, na.rm = TRUE)) + xlab(colnames(phenodf)[i])

combn(c(2:5), 2) #combn(letters[1:4], 2)
cchris<-as.matrix(expand.grid(2:3, 23:24))
cchris<-t(as.matrix(expand.grid(colnames(chrisKidQue[,2:22]),colnames(chrisKidQue[,23:39]))))

#-----------------------------------------------------#
#library(gtools)

pdf('Rplot46-PhenoKidneyViolin.pdf', width=15)
for (i in 2:22){
  for (j in 23:39){
    print(ggplot(chrisKidQue, aes(chrisKidQue[,i], chrisKidQue[,j])) +
            geom_violin(aes(fill = chrisKidQue[,i]), show.legend=FALSE, color = "steelblue") + 
            geom_boxplot(width = 0.08) + theme_minimal() + #+ ylim(0,18)# + scale_x_discrete(labels=c("Male", "Female"))
            xlab(colnames(chrisKidQue)[i]) +  ylab(colnames(chrisKidQue)[j]))
  }
}

dev.off()
#-----------------------------------------------------#
sapply(2:22, function(i)
  table(chrisKidQue[,i]))
#print(colnames(chrisKidQue)[i])



#-----------------------------------------------------#
#--------------- TSH-SNP interaction -----------------
#-----------------------------------------------------#

#Interaction TEST MODEL 
summary(lm(eGFRw.log.Res ~ `chr1:10599281` * TSH + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
           data = vcfReg_TSHmod))$coefficients[c(2:3,16), c(1,2,4)]


#------------------#
regresModel <- function(data, mySNP, formula, vecTraits)
{
  SNP <- data[, mySNP]
  #SNPid <- as.character(noquote(mySNP))
  m <- lm(as.formula(formula), data = data)
  p <- summary(m)$coefficients[vecTraits, c(1,2,4)] #c(2:3,16) applies to TSH; c(2:4,17:18)applies to TSH_cat
  r <- as.data.frame(p) %>%
    rename(Beta = Estimate, SE = `Std. Error`, Pvalue = `Pr(>|t|)`) %>%
    rownames_to_column(var = "Term") %>%
    mutate(associated = ifelse(Pvalue < 0.05/11, "Yes", "No")) %>%
    pivot_wider(#id_cols     = SNPid,
                names_from  = Term,
                values_from = c(Beta, SE , Pvalue, associated),
                names_glue  = "{Term}_{.value}") %>%
    #as.data.frame() %>%
    mutate("SNPid" = mySNP) %>%
    inner_join(repSNPs[c("SNPid", "Locus")],
               by = "SNPid") %>% 
    select(SNPid, Locus, starts_with("SNP_"), starts_with("TSH_"), everything())
  return(r)
}
#---------#

# Testing regression model function
#regresModel(vcfReg_TSHmod,
#            targets[1],  #vcfReg_TSHmod[targets[1:3]]
#            paste("eGFRw.log.Res ~ SNP * TSH + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
#---------#

# Forming principal components term of the model 
PCs <- paste0("PC", 1:10, collapse = " + ")

#---------#
# Iterating regression model function
# Linear interaction: TSH-SNP

results_Kidney_TSH <-
  map_dfr(targets, function(SNP) {
    regresModel(vcfReg_TSHmod,
                SNP, #vcfReg_TSHmod[,SNP],
                paste("eGFRw.log ~ SNP * TSH + Sex + Age +", PCs),
                c(2:3,16))
  }
  )
#---------#

# Save the results of linear interaction model
write.csv(results_Kidney_TSH,
          "13-Feb-23_Interaction of kidney with TSH (eGFRw.log-TSH) adjusted for sex age PCs.csv",
          row.names = FALSE, quote = FALSE)
#---------#

# Table XX of the paper 
results_Kidney_TSH %>%
  filter(SNPid %in% tagSNPs$SNPid) %>% View



#-----------------------------------------------------#
#-------------- SNP:TSH_cat interaction ---------------
#-----------------------------------------------------#


# Non-linear interaction: TSH_cat-SNP

results_Kidney_TSH_cat <-
    map_dfr(targets, function(SNP) {
      regresModel(vcfReg_TSHmod,
                  SNP,
                  paste("eGFRw.log ~ SNP * TSH_cat + Sex + Age+", PCs),
                  c(2:4,17:18))
    }
    ) %>% 
  rename_with(~ gsub("TSH_catHyperT", "Hyperthyrodism", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("TSH_catHypoT",  "Hypothyrodism",  .x, fixed = TRUE)) %>%
  select(SNPid, Locus, starts_with(c("SNP_",
                                     "Hyperthyrodism_",
                                     "Hypothyrodism_",
                                     "SNP:Hyperthyrodism_",
                                     "SNP:Hypothyrodism_")))
#---------#

# Save the results of non-linear interaction model
write.csv(results_Kidney_TSH_cat,
          "13-Feb-23_Interaction of kidney with thyroid disease (eGFRw.log-TSH_cat) adjusted for sex age PCs.csv",
            row.names = FALSE, quote = FALSE)
#---------#

# Table XX of the paper 
results_Kidney_TSH_cat %>% 
  filter(SNPid %in% tagSNPs$SNPid) %>% View






#-----------------------------------------------------#
#--------------- TSH_cen-SNP interaction -------------
#-----------------------------------------------------#

# Model the interaction in different levels of TSH for STC1 locus

# By default, scale() function will standardize the 
# data (mean zero, unit variance). To indicate that
# we just want to subtract the mean, we need to turn 
# off the argument scale = FALSE.

#results_Kidney_TSH_cen <- 
vcfReg_TSHmod %>%
  select(AID, eGFR, eGFRw.log, TSH, TSH_cat, Age, Sex, starts_with(c("chr8", "PC"))) %>% 
  # mutate(eGFR.log = log(eGFR),
  #        #TSH_cen  = scale(TSH, scale = F)
  #        ) %>% 
  pivot_longer(cols      = starts_with("chr"),
               names_to  = c("SNPid"),
               values_to = c("Dosage")) %>%
  group_by(SNPid) %>%
  nest() %>%
  mutate(
    model_tsh = data  %>% map(~lm(eGFRw.log ~ Dosage * TSH     + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)),
    model_cat = data  %>% map(~lm(eGFRw.log ~ Dosage * TSH_cat + Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = .)),
    tidy_tsh  = model_tsh %>% map(broom::tidy),
    tidy_cat  = model_cat %>% map(broom::tidy)) %>%
  unnest(tidy_tsh) %>%
  #unnest(tidy_cat) %>%
  filter(!str_detect(term, "(Intercept)|PC|Age|Sex")) %>%
  select(-data, -model_tsh, -statistic) %>%
  inner_join(repSNPs[c("Locus", "SNPid")],
            .,
            by = "SNPid") %>%
  mutate(associated = ifelse(p.value < 0.05/11, "Yes", "No")) %>% View
  #filter(associated == "Yes") %>% View
  filter(str_detect(term, "Dosage:")) %>% View
  #write.csv("01-Feb-23_Kideny interaction with centralized thyroid (SNP-TSH)_long format.csv", row.names = FALSE, quote = FALSE)
  pivot_wider(id_cols = c(Locus, SNPid),
              names_from = term,
              names_glue = "{term}_{.value}",
              values_from = c(estimate,
                              std.error,
                              #statistic,
                              p.value,
                              associated))
#---------#
write.csv(results_Kidney_TSH_,
          "01-Feb-23_Kideny interaction with thyroid adj for sex and age (SNP-TSH)_wide format.csv",
          row.names = FALSE, quote = FALSE)
#---------#

pdf('02-Feb-23_interaction comparisons.pdf', width=18, height = 18)

map2(c("Dosage_estimate_nonCen", 
       "TSH_cen_estimate", 
       "`Dosage:TSH_nonCen_estimate`", 
       "Dosage_p.value_nonCen.log", 
       "TSH_nonCen_p.value.log", 
       "Dosage_TSH_nonCen_p.value.log"),
     c("Dosage_estimate_Cen", 
       "TSH_cen_estimate",
       "`Dosage:TSH_cen_estimate`", 
       "Dosage_p.value_Cen.log", 
       "TSH_cen_p.value.log", 
       "Dosage_TSH_cen_p.value.log"),
     function(myvar1, myvar2) {
       results_Kidney_TSH_ %>%
         rename_with(~ gsub("TSH", "TSH_nonCen", .x, fixed = TRUE))%>%
         right_join(results_Kidney_TSH_cen,
                    by = c("Locus", "SNPid"),
                    suffix = c("_nonCen", "_Cen")) %>%
         write.csv("02-Feb-23_Kideny interaction with and without centralized thyroid (SNP-TSH)_wide format.csv", row.names = FALSE, quote = FALSE)
         mutate(Dosage_p.value_nonCen.log = -log10(Dosage_p.value_nonCen),
                TSH_nonCen_p.value.log    = -log10(TSH_nonCen_p.value),
                Dosage_TSH_nonCen_p.value.log = -log10(`Dosage:TSH_nonCen_p.value`),
                Dosage_p.value_Cen.log    = -log10(Dosage_p.value_Cen),
                TSH_cen_p.value.log       = -log10(TSH_cen_p.value),
                Dosage_TSH_cen_p.value.log= -log10(`Dosage:TSH_cen_p.value`)) %>%
       ggplot(aes_string(x = myvar1, y = myvar2, color = "Locus")) +
           geom_abline(lty = 2) +
           geom_vline(xintercept = 0) +
           geom_hline(yintercept = 0) +
           geom_point(size = 2.5) +
           labs(x = paste(myvar1, "in ln(eGFR) ~ Dosage * TSH + Sex + Age + PC1:10"),
                y = paste(myvar2, "in ln(eGFR) ~ Dosage * TSH_cen + Sex + Age + PC1:10")) +
           theme_classic()
     })

dev.off()

# ggsave("02-Feb-23_Effect of SNP in model with TSH vs centralized TSH.png",
#        width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

#---------#

#-----------------------------------------------------#
#--------------- Depicting interaction ---------------
#-----------------------------------------------------#


# Visualizing SNP-TSH interaction
lapply(c("chr8:23885208",
         "chr8:23894869"), function(SNP)
  {
   vcfReg_TSHmod %>%
             select(eGFRw.log, SNP, TSH, Age, Sex, starts_with("PC")) %>% #starts_with("chr8:")
             #filter(!is.na(TSH_cat)) %>%
             mutate(#TSH_cen = scale(TSH, scale = F),
                    TSH_third = case_when(TSH < quantile(TSH, 1/3, na.rm = T) ~ "Hyperthyrodism",
                                          TSH > quantile(TSH, 2/3, na.rm = T) ~ "Hypothyrodism",
                                          TRUE ~ "Normal TSH")) %>% 
            mutate(across(c(Age, starts_with("PC")), ~ scale(., scale = F))) %>%
            ggplot(aes(x = !!sym(SNP), y = eGFRw.log, color = TSH_third)) +
            geom_smooth(method = "lm", se = T) + 
             #facet_wrap(~TSH_cat, nrow = 3) +
             scale_color_manual(name = "Stratified TSH",
               values = c("steelblue3", "green4", "turquoise2"),
               breaks = c("Hyperthyrodism", "Normal TSH", "Hypothyrodism")) +
             labs(x = paste("Dosage of", SNP, "variant at STC1"),
                  y = "log(eGFRcreat)") +
             theme_classic()
           
           
   # Save the plot
   ggsave(paste0("15-Feb-23_SNP-TSH interaction effect on eGFRw.log for ", str_replace(SNP, ":", "_"), " at STC1.png"),
          last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")      
   }
  )


#geom_point()+
#geom_violin()+
#geom_boxplot(aes(fill = TSH_cat), color = "grey", width = .01, position = position_dodge(2)) +


#---------#
library(interactions)
library(jtools)
#---------#

# Emerging SNP-TSH interaction via package `chr8:23885208`
int_SNP2 <-
  interact_plot(
  lm(paste0("eGFRw.log ~ `chr8:23885208` * TSH + Sex + Age +", PCs), vcfReg_TSHmod),
  pred = `chr8:23885208`,
  modx = TSH,
  interval = T,
  linearity.check = F,
  plot.points = F,
  #colors = "seagreen",
  #legend.main = "Custom Legend Title",
  x.label = "Dosage for 8:23885208 variant at *STC1* locus",
  y.label = "log(eGFRcreat)") +
  theme_classic() +
  #jtools::theme_apa()
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.placement = "outside",
        axis.text.x  = element_text(size = 8,  face = "bold"),
        axis.text.y  = element_text(size = 7,  face = "bold"),
        axis.title   = element_text(size = 14, face = "bold"),
        axis.title.x = ggtext::element_markdown(),
        legend.position = c(.8, .22),
        legend.key.size  = unit(0.8, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))
  

# 1st way
library(gridExtra)
grid.arrange(int_SNP2, int_SNP1, ncol = 2, labels = c("A", "B"))

# 2nd way
ggpubr::ggarrange(int_SNP2, int_SNP1, #bp + rremove("x.text"), 
          labels = c("A", "B"), ncol = 2, nrow = 1)

# 3rd way
cowplot::plot_grid(int_SNP2, int_SNP1, labels = c("A", "B"), ncol = 2, nrow = 1)


ggsave("20-Feb-23_SNP-TSH interactions chr8_23885208 and chr8_23894869 at STC1.png",
       width = 11, height = 5, pointsize = 5, dpi = 300, units = "in")




#-----------------------------------------------------#
#------------- Interaction with Magnesium ------------
#-----------------------------------------------------#

quantile(vcfReg$Magnesium_mg, na.rm = T, probs = seq(0,1,.1))

#---------#
#scatter plot for Magnesium_mg and PTT_sec interaction
ggplot(data = vcfReg, aes(x = Magnesium_mg, y = Age))+
  geom_point(color = "steelblue", alpha = 0.6) + 
  theme_classic()

#Scatter plot for PTT vs aPTT and its family traits in CHRIS
pdf("13-Sep-2022_Scatter plot of PTT_sec vs. others")
lapply(c("Prothrombin_Time",
         "INR_PT_INR",
         "APTT_ratio"), function(i)
           ggplot(data = vcfReg, aes_string(x = "PTT_sec", y = i))+
         geom_point(color = "steelblue", alpha = 0.6) + 
         theme_classic()+
         labs(x="aPTT_sec"))
dev.off()

#---------#
#histogram
ggplot(vcfReg, aes(PTT_sec)) +  
  geom_histogram(color = "steelblue2", fill = "steelblue3", bins = 100, alpha = 0.7)+
  theme_classic()

#---------#
#testing model for magnesium interaction
vcfReg_Mag <- 
  vcfReg %>% 
  mutate(Magnesium_cat = cut(Magnesium_mg,
                             right  = T, 
                             breaks = c(-Inf, 1.459, 2.680, Inf),
                             #labels = c("Hypomagnesemia", "Normal", "Hypermagnesemia")
                             labels = c("0", "1", "2")),
         Magnesium_cat = relevel(Magnesium_cat, ref = 2),
         PTT_cat       = cut(PTT_sec,
                             right = T, 
                             breaks = c(-Inf, 20, 35, Inf),
                             labels = c("0", "1", "2"))) #%>% count(Magnesium_cat)

#---------#
map2(vcfReg_Mag[c("chr4:76489165", "chr4:76490987", "chr5:177386403")],
     vcfReg_Mag[c("Magnesium_cat", "Magnesium_cat", "PTT_cat")],
     function(SNP, TRAIT) summary(lm(eGFRw.log.Res ~ SNP * TRAIT,
                                     data = vcfReg_Mag,
                                     na.action = na.exclude))) #%>% do.call(rbind, .)

library(jtools)

summary(lm(eGFRw.log.Res ~ `chr4:76444730` + Magnesium_mg, data = vcfReg_Mag))
summary(lm(eGFRw.log.Res ~ `chr4:76444730` * Magnesium_mg, data = vcfReg_Mag))
summary(lm(eGFRw.log.Res ~ `chr4:76490987` * Magnesium_mg, data = vcfReg_Mag))
summary(lm(Magnesium_mg ~ `chr4:76490987`,                 data = vcfReg_Mag))
summary(lm(Magnesium_mg ~ `chr4:76444730` + Age + Sex,     data = vcfReg))
summary(lm(Magnesium_mg ~ `chr4:76447694` + eGFRw.log.Res + Age + Sex, data = vcfReg))
#---------#

getInteractCoefs2table <- 
  function(mytrait, mytarget, myformula, myparameter){
    results1 <- lapply(mytrait,
                       function(trait){
                         map_df(mytarget,
                                function(SNP){
                                  #myformula <- as.formula(eGFRw.log.Res ~ SNP + trait + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
                                  m <- lm(as.formula(myformula), data = vcfReg)
                                  #s <- coef(m)[2]
                                  p <- summary(m)$coefficients[myparameter, c(1,2,4)] #Beta, se, Pvalue
                                  #t <- tidy(m)[2, c(1,2,4)]#broom
                                })
                       })
    results2 <- as.data.frame(do.call(cbind, results1)) #%>% clean_names() #library(janitor)
    #names(results2) <- gsub(x = names(results2), pattern = "Pr\\([^\\(]*\\)", replacement = "_")
    names(results2) <- str_replace_all(names(results2), c("Estimate"="Beta","Std. Error"="SE","Pr\\([^\\(]*\\)"="Pvalue"))
    #colnames(results2) <- names(mytrait)
    results3 <- cbind(SNPid = names(mytarget), results2)
    return(results3)
  }


Mag_SNP <- 
  getInteractCoefs2table(vcfReg[c("Magnesium_mg")], 
                         vcfReg[targets],#, "PTT_sec"
                         paste("eGFRw.log.Res ~ SNP * trait "), 
                         2) %>% #+ Age + Sex
  rename(SNP.Beta   = Magnesium_mg.Beta, 
         SNP.SE     = Magnesium_mg.SE, 
         SNP.Pvalue = Magnesium_mg.Pvalue)

Mag_trait <- getInteractCoefs2table(vcfReg["Magnesium_mg"], 
                                    vcfReg[targets],
                                    paste("eGFRw.log.Res ~ SNP * trait"), 
                                    3) %>% 
  rename(Magnesium.Beta   = Magnesium_mg.Beta,
         Magnesium.SE     = Magnesium_mg.SE, 
         Magnesium.Pvalue = Magnesium_mg.Pvalue)

Mag_interaction <- 
  getInteractCoefs2table(vcfReg["Magnesium_mg"], 
                         vcfReg[targets],
                         paste("eGFRw.log.Res ~ SNP * trait "), 
                         4) %>% 
  rename(Interaction.Beta   = Magnesium_mg.Beta,
         Interaction.SE     = Magnesium_mg.SE, 
         Interaction.Pvalue = Magnesium_mg.Pvalue)

#---------#
Mag_with_covs <-
  repSNPs %>% 
  select(SNPid, Locus) %>% 
  inner_join(Mag_SNP, by = "SNPid") %>%
  inner_join(Mag_trait, by = "SNPid") %>% 
  inner_join(Mag_interaction, by = "SNPid") %>%
  #inner_join(Mag_Step4, by = "SNPid") %>% #View()
  #write.csv(., "01-Sep-2022_Interaction between SNPs and Serum Magnesium (quantitative).csv", row.names = FALSE)
  filter(Interaction.Pvalue < 0.01) %>% 
  as_tibble() %>% View()

#---------#
getInteractCoefs2table(vcfReg[c("Magnesium_mg")],
                       vcfReg[targets],
                       paste("eGFRw.log.Res ~ SNP + trait "), 2) %>% 
  rename(SNP.Beta   = Magnesium_mg.Beta,
         SNP.SE     = Magnesium_mg.SE,
         SNP.Pvalue = Magnesium_mg.Pvalue) %>% 
  inner_join(
    getInteractCoefs2table(vcfReg[c("Magnesium_mg")],
                           vcfReg[targets],
                           paste("eGFRw.log.Res ~ SNP + trait + Age + Sex"), 2) %>%
      rename(SNP_with_covs.Beta   = Magnesium_mg.Beta, 
             SNP_with_covs.SE     = Magnesium_mg.SE, 
             SNP_with_covs.Pvalue = Magnesium_mg.Pvalue)) %>% 
  ggplot(aes(SNP.Beta, SNP_with_covs.Beta)) + 
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = +1, color = "grey50", size = 0.5, lty = 1)+
  theme_classic()

#---------#
#Step2: testing Magnesium as a mediator for eGFR 
Mag_Step2 <- getInteractCoefs2table(vcfReg["Magnesium_mg"],
                                    vcfReg[targets],
                                    paste("eGFRw.log.Res ~ SNP + trait"),
                                    2) %>% 
  rename(MagnesiumAdj.Beta   = Magnesium_mg.Beta, 
         MagnesiumAdj.SE     = Magnesium_mg.SE, 
         MagnesiumAdj.Pvalue = Magnesium_mg.Pvalue)

#Step3: testing direct association of SNPs with Magnesium
Mag_Step3 <- getInteractCoefs2table(vcfReg["Magnesium_mg"], 
                                    vcfReg[targets],
                                    paste("trait ~ SNP + Age + Sex"),
                                    2) %>%
  rename(SNP.Beta   = Magnesium_mg.Beta, 
         SNP.SE     = Magnesium_mg.SE, 
         SNP.Pvalue = Magnesium_mg.Pvalue)

#Step4: testing eGFR as a mediator for Magnesium
Mag_Step4 <- getInteractCoefs2table(vcfReg["eGFRw.log.Res"], 
                                    vcfReg[targets],
                                    paste("Magnesium_mg ~ SNP + trait + Age + Sex"),
                                    2) %>%
  rename(eGFRadj.Beta   = eGFRw.log.Res.Beta, 
         eGFRadj.SE     = eGFRw.log.Res.SE, 
         eGFRadj.Pvalue = eGFRw.log.Res.Pvalue)
#---------#

# Merging Step3 and Step4 of Magnesium models
# for testing mediation effect of eGFR
Mag_Step2 %>% 
  inner_join(Mag_Step3,
             by = "SNPid") %>% 
  inner_join(Mag_Step4,
             by = "SNPid") %>% 
  inner_join(repSNPs[c("SNPid",
                       "Locus",
                       "MARKER_ID",
                       "BETA",
                       "SEBETA",
                       "PVALUE")],
             by = "SNPid") %>%
  rename(Beta_GWAS = BETA, SE_GWAS = SEBETA, Pvalue_GWAS = PVALUE) %>%
  select(SNPid, Locus, MARKER_ID, Beta_GWAS, SE_GWAS, Pvalue_GWAS, everything()) %>% 
  select(- contains("SE")) %>% 
  mutate(Effect_Change = (SNP.Beta - eGFRadj.Beta)/ SNP.Beta,
         MARKER_ID = noquote(str_extract(repSNPs$MARKER_ID, "chr[0-9]+:[0-9]+_[A-Z]+/[A-Z]+"))) %>% #View()
  #write.csv(., "13-Sep-22_Mediatory effect of eGFR on association of SNPs with Magnesium (Steps 1 to 4).csv", row.names = FALSE)
  ggplot(aes(SNP.Beta, eGFRadj.Beta)) +
  geom_point(aes(color = Locus), size = 3, alpha = .75) +
  geom_abline(intercept = 0, slope = +1, color = "grey50", size = 0.5) +
  theme_classic() +
  labs(x = "SNP effect in lm(Magnesium ~ SNP + Age + Sex)",
       y = "SNP effect in lm(Magnesium ~ SNP + eGFRw.log.Res + Age + Sex)")

ggsave("13-Sep-22_Mediatory effect of eGFR on association of SNPs with Magnesium.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")



#-----------------------------------------------------#
#---------------------- Other ------------------------
#-----------------------------------------------------#

# retrieve & save a printed but unsaved console output history by: .Last.value
View(.Last.value)  
history(max.show = Inf) # full history
#---------#

#Plotting and comparing Betas and Pvalues of interaction terms
resModel3 %>% 
  as.data.frame() %>% 
  mutate(coef = rownames(resModel3)) %>% 
  rename(Beta = Estimate, SE = `Std. Error`, Pvalue = `Pr(>|t|)`) %>% 
  mutate(SNPid = rep(repSNPs$SNPid, each = 5),
         Locus = rep(repSNPs$Locus, each = 5)) %>% 
  mutate(coef = replace(coef, coef == "TSH_catHyperT",     "HyperThyrodism"),
         coef = replace(coef, coef == "TSH_catHypoT",      "HypoThyrodism"),
         coef = replace(coef, coef == "SNP:TSH_catHyperT", "SNP_HyperThyrodism"),
         coef = replace(coef, coef == "SNP:TSH_catHypoT",  "SNP_HypoThyrodism")) %>%
  mutate(tagHyper = if_else(coef == "SNP_HyperThyrodism", "GAB2", ""),
         tagHypo  = if_else(coef == "SNP_HypoThyrodism",  "GAB2      IGF1R, PIP5K1B", "")) %>% 
  filter(coef %in% c("SNP_HyperThyrodism", "SNP_HypoThyrodism")) %>% 
  ggplot(aes(x = Pvalue)) +
  geom_density(data = data.frame("rnd" = runif(100000)),
               aes(x = rnd),
               fill = "grey70",
               color = "grey50")+
  geom_density(aes(fill = coef),
               alpha = 0.6,
               color = "grey50") + #, show.legend = NULL
  scale_fill_manual(values = c("violetred1", "turquoise2"),
                    #name   = "SNP vs",
                    labels = c("Hyperthyrodism",
                               "Hypothyrodism")) + 
  geom_vline(aes(xintercept = 0.05), lty = 2) +
  geom_text(aes(label = tagHyper),
            x = 0.14,
            y = 7,
            color = "violetred1",
            size = 4,
            fontface = 2) +
  geom_text(aes(label = tagHypo),
            x = 0.68,
            y = 3,
            color = "turquoise2",
            size = 4,
            fontface = 2) +
  labs(title = "Density plot of the p-values for interaction of SNPs vs:",
       x = "Pvalue for the interaction term") + 
  facet_grid(coef ~ .) + #xlim(0, 0.0025) + scale_fill_discrete() +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank(), #"",
        legend.position = "top", #"none", "bottom"
        axis.text.y  = element_text(size = 10, face = "bold"),
        axis.text.x  = element_text(size = 10, face = "bold"),
        axis.title   = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.6, 'cm'),
        legend.text  = element_text(size = 10, face = "bold"),
        strip.text = element_blank()
  )

ggsave("07-Jun-222_SNP-Hyperthyroid_vs_SNP-Hypothyroid_interactions_Pvalues_ESHG_poster.png", 
       width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-----------------------------------------------------#

#Plots
ggplot(resAllModels, aes(y = (SNP_Model0 - T3_beta)))+
  geom_boxplot(alpha = 0.6, fill = "Orange", outlier.color = "red")
#geom_text(aes(label = Locus, color = Locus), size = 3, vjust = -1)

#scatter
ggplot(resAllModels, aes(x = TSH_beta , y = TSH.q_beta))+
  geom_point(aes(shape = Locus, color = Locus), size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) +
  geom_abline(slope = 1, intercept = 0, color="gray", linetype="dashed", size = 1.1) +
  theme_classic()
#geom_text(aes(label = Locus, color = Locus), size = 2, vjust = -1)

ggplot(resAllModels, aes(x = (SNP_Model0 - T3_beta) , y = (Effect.ckdgen)))+
  geom_point(aes(shape = Locus, color = Locus), size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) + theme_classic()
#geom_abline(slope = 1, intercept = 0, color="gray", linetype="dashed", size = 1.1) + 

ggsave("SNP_Model0-T3_beta) vs. (Effect.ckdgen).png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#---------#

library(pheatmap)
pdf('pheatmap_MediationAnalysis_TSH_cat_25Jan.pdf', width=18, height = 18)
#Betas
pheatmap(medAna_Merged[c("BETA", "SNP.Beta2", "SNP.Beta5")] / medAna_Merged$BETA,
         cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = T, labels_row = medAna_Merged$SNPid, 
         border_color = 'Black', fontsize_row = 6, fontsize_col = 12, angle_col = "0")
#P-values
pheatmap(medAna_Merged[c("PVALUE", "SNP.Pvalue2", "SNP.Pvalue5")],
         cluster_cols = FALSE, cluster_rows = FALSE, 
         show_rownames = T, labels_row = medAna_Merged$SNPid, 
         border_color = 'Black', fontsize_row = 6, fontsize_col = 12, angle_col = "0")
#---------#

#Scatter plots of Betas
p11 <- medAna_Merged %>% 
  ggplot(aes(BETA, SNP.Beta2)) + 
  geom_point(aes(shape = Locus, color = Locus), size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) + 
  geom_abline(slope = 1, intercept = 0, color="gray") + theme_classic()

p12 <- medAna_Merged %>% 
  ggplot(aes(BETA, SNP.Beta5)) + 
  geom_point(aes(shape = Locus, color = Locus), size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) + 
  geom_abline(slope = 1, intercept = 0, color="gray") + theme_classic()

library(gridExtra)
p13 <- grid.arrange(p11, p12, ncol = 1)

dev.off()
#-----------------------------------------------------#

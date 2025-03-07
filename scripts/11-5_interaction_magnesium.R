



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



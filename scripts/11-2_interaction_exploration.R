
library(tidyverse)
library(knitr) #for printing html-friendly tables.

# inputs

# outputs

#-------------#
# Find frequency of cases with thyroid disorder
# group_by(TSH_cat, Thyroid_DrugName) %>%
# summarise(n = n()) %>%
# spread(TSH_cat, n)



#-----------------------------------------------------#
#-------         Emerging baches effect        -------
#-----------------------------------------------------#

# TSH vs. TSH Instruments
chris %>% 
  select(TSH, TSH.Ins) %>% 
  pivot_wider(names_from = TSH.Ins, values_from = TSH, names_prefix = "TSH_") %>% 
  head()

#-------------#
# Some randome ideas for interaction


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


#-----------------------------------------------------#
#------          Thyroid Questionnaire         -------
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
#--------         Kidney Questionnaire         -------
#-----------------------------------------------------#

phenodf$AID <- phenodf$FAM_ID
Kidney$AID  <- as.character(Kidney$AID)
chrisKidQue <- merge(Kidney, phenodf, by = "AID")


chrisKidQue[] <- lapply(chrisKidQue, function(x) if(is.factor(x)) as.character(x) else x)


#-----------------------------------------------------#
#library(gtools)

pdf('Rplot46-PhenoKidneyViolin.pdf', width=15)

for (i in 2:22){
  for (j in 23:39){
    print(
      ggplot(
        data = chrisKidQue, 
        aes(
          chrisKidQue[,i], 
          chrisKidQue[,j]
          )
        ) +
        geom_violin(aes(fill = chrisKidQue[,i]), show.legend=FALSE, color = "steelblue") + 
        geom_boxplot(width = 0.08) + theme_minimal() + #+ ylim(0,18)# + scale_x_discrete(labels=c("Male", "Female"))
        labs(
          x = colnames(chrisKidQue)[i],
          y = colnames(chrisKidQue)[j]
          )
      )
  }
}

dev.off()

#-----------------------------------------------------#
sapply(2:22, function(i) table(chrisKidQue[,i]))



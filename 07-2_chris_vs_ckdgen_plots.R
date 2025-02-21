#=========================================#
#  Visualizing Replication vs. Discovery
#=========================================#


# Started on October 2021
# Last update on February 20, 2025
library(ggfittext)
library(tidyverse)

#-----------------------------------------------------#
#------    Figure 3B: MAF vs. Effect Ratio     -------
#-----------------------------------------------------#

# inputs
egfr_log_res_EA <- "~/projects/kidneyInCHRIS/inputs/11-Dec-23_Suppl._Table_3_163_replicated_SNPs_in_CHRIS_EUA.csv"

# outputs
effects_rangeplot <- "~/projects/kidneyInCHRIS/outputs/21-Feb-25_Figure_2_effect_comparison_leading_variants.png"


# read supplementary Table 3
repSNPs <- data.table::fread(egfr_log_res_EA)

#-------------#
# Beta ratio vs. MAF ratio in CKDGen and CHRIS 
fig_3b <- repSNPs %>% 
  ggplot(aes(x = MAF_Ratio, y = Beta_Ratio, shape = Locus, color = Locus)) + 
  geom_point(alpha = 0.9, size = 3.5) +
  #ggrepel::geom_text_repel(data = tagSNPs, aes(label = Locus), check_overlap = TRUE) +
  #ggfittext::geom_fit_text(grow = TRUE)+
  #geom_text(aes(label = Locus), check_overlap = TRUE) +
  #geom_abline(intercept = 0, slope = +1, color = "grey50", size = .9)+
  scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  scale_color_manual(values= c("maroon1", "darkorchid2", "orange2", "green4", "steelblue2", "darkturquoise",
                               "tomato", "springgreen2", "royalblue2", "gold", "grey50")) +
  #scale_color_brewer(palette = "Dark2") +
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 1)+
  scale_y_continuous(breaks = seq(0,6, 0.5)) +
  scale_x_continuous(breaks = seq(0.75,1.20, .05), limits = c(.75,1.20), expand = c(0,0)) +
  annotate("segment", x = 0.755, xend = 0.995, y = 6, yend = 6, colour = "grey50", size = 1.2,
           arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm"))) +
  annotate("segment", x = 1.005, xend = 1.195, y = 6, yend = 6, colour = "grey50", size = 1.2,
           arrow = arrow(ends = "both", angle = 90, length = unit(0.2, "cm"))) +
  annotate("text", x = 0.85, y = 6.3, size = 4,
           label = "Rarer alleles in CHRIS") + #\nLarger effect in CHRIS
  annotate("text", x = 1.09, y = 6.3, size = 4, 
           label = "Rarer alleles in CKDGen") +#\nLarger effect in CHRIS
  # annotate("text", x = 0.23, y = 5.6, size = 4, label = "CASZ1") +
  # annotate("text", x = 0.7, y = 5.1, size = 4, label = "CASZ1") +
  # annotate("text", x = 0.59, y = 4, size = 4, label = "GAB2") +
  # annotate("text", x = 0.65, y = 3, size = 4, label = "DDX1") +
  # annotate("text", x = 1.15, y = 4, size = 4, label = "IGF1R") +
  # annotate("text", x = 1.2, y = 2.8, size = 4, label = "PIPK1B") +
  labs(x = "\nCHRIS-to-CKDGen MAF ratio\n", 
       y = "\nCHRIS-to-CKDGen effect ratio\n") + 
  #theme_minimal()+ #theme_light(base_size = 10)+ #lims(y = c(1, 6)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
        panel.grid.major.x = element_blank(),
        axis.text  = element_text(size = 11,  face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        #legend.position = "none",
        legend.key.size  = unit(0.99, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.text  = element_text(size = 12, face = "italic"),
        legend.title = element_text(size = 13, face = "bold"))

ggsave("22-Mar-23_MAF Ratio vs Effect Size Ratio in CHRIS and CKDGen.png", 
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")



#-----------------------------------------------------#
#------          Figure 2: Line range           ------
#-----------------------------------------------------#

# The difference of Effect size in CHRIS and CKDGen

plt_effect_comparison <- repSNPs %>%
  filter(SNPid %in% (target_snps %>% str_remove("chr"))) %>%
  select(SNPid, Locus, Beta_CHRIS_ald, Beta_CKDGen_ald, SE_CHRIS, SE_CKDGen) %>% 
  rename(Beta_CHRIS  = Beta_CHRIS_ald, 
         Beta_CKDGen = Beta_CKDGen_ald) %>% 
  pivot_longer(cols = !c(SNPid, Locus),   
               names_to = c("trait", "study"),
               names_pattern = "(.+)_(CHRIS|CKDGen)$",
               values_to = c("score")) %>%
  pivot_wider(names_from = "trait",
              values_from = "score") %>%
  ggplot(aes(x = locus_factor_rsid(Locus),
             y = Beta,
             color = study,
             ymin = Beta - 1.95*SE,
             ymax = Beta + 1.95*SE)) + 
  geom_pointrange(aes(color = study),
                  size = 1.1, fatten = 2.5, alpha = 0.9,
                  #lineend = "round",
                  position = position_dodge(.3)) +
  #scale_color_grey(start=0.55, end=0.25) +
  scale_color_manual(values=c('turquoise3','tomato3'))+
  #scale_size_manual(values=c(2, 4))+
  geom_hline(yintercept = 0, color="black", lty = 1, size = .5) +
  scale_y_continuous(breaks = seq(-.035, .001, .005)) +
  #theme_light(base_size = 10) +
  labs(y = "Effect on ln(eGFRcrea)\n", x = NULL) +
  coord_cartesian(ylim = c(-0.035, 0.001)) +
  theme(legend.title = element_blank(),
        legend.position = c(.9, .365),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'solid', color = "grey80", size = .15),
        axis.text.y  = element_text(size = 8,  face = "bold"),
        axis.text.x  = element_text(size = 6,  face = "bold.italic"),
        axis.title   = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.6, 'cm'),
        legend.text  = element_text(size = 10))

ggsave(effects_rangeplot, plot = plt_effect_comparison, width = 7.5, height = 4.2, dpi = 500, units = "in")

#-----------------------------------------------------#

#filling the locus information missed for some of the replicated variants

tempResult <- tempResult %>%
  mutate(Locus = replace(Locus, which(CHR == '5' & BEG == '39385539' ), "DAB2"))%>%
  mutate(Locus = replace(Locus, which(CHR == '5' & BEG == '39393631' ), "DAB2"))%>%
  mutate(Locus = replace(Locus, which(CHR == '5' & BEG == '39397030' ), "DAB2"))%>%
  mutate(Locus = replace(Locus, which(CHR == '7' & BEG == '77714744' ), "RSBN1L"))%>%
  mutate(Locus = replace(Locus, which(CHR == '7' & BEG == '77736048' ), "RSBN1L"))%>%
  mutate(Locus = replace(Locus, which(CHR == '8' & BEG == '23885208' ), "STC1"))
#-----------------------------------------------------#

p1 <- 
  repSNPs %>% 
  #tempResult %>% 
  #drop_na(any_of(propVars))
  drop_na(any_of(c("Effect.ckdgen",
                   "Effect.CHRIS",
                   "CHRIS.CKDGen.Effect.Ratio"))) %>%
  ggplot(aes(x = Effect.ckdgen,
             y = Effect.CHRIS)) + 
  geom_point(aes(shape = Locus,
                 color = Locus),
             size = 2) +
  scale_shape_manual(values = c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17)) +
  #scale_color_manual(values=c( ' #999999 ' , ' #E69F00 ' , ' #56B4E9 ' ))+
  #scale_size_manual(values=c(1.5, 2, 3))+
  #theme(legend.position="top")
  #geom_quantile()+
  #xlim(-0.0074, -0.002)+
  geom_abline(intercept = 0,
              slope = +1,
              color = "grey50",
              size = 1.5,
              lty = 2)+
  geom_vline(xintercept=0)+#, linetype="dashed", color = "darkgray", size=1.5)+
  geom_hline(yintercept=0)+#, linetype="dashed", color = "darkgray", size=1.5)+
  xlab("Effect CKDGen")+#"Association coefficient for log(eGFRcrea) in CKDGen"
  ylab("Effect CHRIS")+#"Association coefficient for log(eGFRcrea) in CHRIS"
  theme_classic()+
  theme(axis.text  = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.key.size = unit(0.65, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16))

ggsave("effectDiff_eGFRw_slides.png",last_plot(), width = 8, height = 5.5, 
       pointsize = 5, dpi = 300, units = "in")
#-------------#

p2 <- 
  repSNPs %>% 
  #tempResult %>% 
  #drop_na(any_of(propVars))
  drop_na(any_of(c("Effect.ckdgen",
                   "Effect.CHRIS",
                   "CHRIS.CKDGen.Effect.Ratio"))) %>%
  ggplot(aes(x = Effect.ckdgen,
             y = MAF.Diff)) +
  geom_point(aes(color = Locus),
             position  = position_jitter(height = 0L, seed = 1L)) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) + 
  geom_linerange(aes(x = Effect.ckdgen,
                     ymax = MAF.Diff,
                     ymin = 0,
                     color= Locus),
                 position = position_jitter(height = 0L, seed = 1L))+
  xlab("Effect CKDGen") + #"Association coefficient for log(eGFRcrea) in CKDGen"
  ylab("MAF CHRIS - MAF CKDGen") + #"Diffrence of Alleles frequency between CHRIS and CKDGen"
  theme_classic()+
  theme(axis.text  = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.key.size = unit(0.65, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16))
#theme(panel.background = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"),
#      axis.text.x = element_text(#face="bold",color="steelblue",
#  size=8, angle=90))#+
#theme(axis.title.x = element_text(size=12))#,face="italic"
#-------------#

library(gridExtra)
p3 <- grid.arrange(p1, p2, ncol=1)

ggsave("21-Apr-22_CKDGen Betas versus MAF.Diff & CHRIS Betas_4Poster.png", 
       p3, width = 10, height = 6.5, pointsize = 8, dpi = 600, units = "in")
#-----------------------------------------------------#
library(reshape2)
#melt((tempResult[, c("Locus", "MAF.ckdgen", "MAF.log.Std" ,"MAF.Diff")]), id.vars = )
tempResult2 <- data.frame("Study"=rep(c("CKDGen", "CHRIS"), each=nrow(tempResult)), 
                          "MAF"=rbind(as.matrix(tempResult$MAF.ckdgen), as.matrix(tempResult$MAF)),
                          "MAF.Diff"=rep(tempResult$MAF.Diff, each=2))
#-------------#
ggplot(tempResult2, aes(x=MAF.Diff, y=MAF))+ geom_point(aes(color=Study), alpha=0.6, size=3)+
  geom_vline(xintercept=0,linetype="dashed") + geom_hline(yintercept=0, linetype="dashed")+
  theme(panel.background = element_blank(), panel.border = element_blank())
#-------------#

dim(data.frame("Study"=rep(c("CKDGen", "CHRIS"), each=nrow(tempResult))))
#-----------------------------------------------------#
library(tidyverse)

tempResult %>% 
  drop_na(any_of(propVars)) %>%
  ggplot(aes(x=Effect.ckdgen, y=CHRIS.CKDGen.Effect.Ratio))+ 
  geom_point(aes(color=Locus, shape=Locus), size=2)+
  geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
  #geom_text(aes(label= Locus), size = 1.5)+
  scale_shape_manual(values=c(0, 1, 2,3, 7, 8, 9, 15, 16, 17))+ theme_bw()#, 15, 16, 17

ggsave("CKDGen Beta vs CHRIStoCKDGen.Effects Ratio.png",
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-------------#
tempResult %>% 
  drop_na(CHRIS.CKDGen.Effect.Ratio) %>%
  ggplot(aes(MAF.Diff, CHRIS.CKDGen.Effect.Ratio))+ 
  geom_point(aes(shape=Locus, color=Locus), size=2)+
  geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
  scale_shape_manual(values=c(0, 1, 2, 3, 7, 8, 9, 15, 16, 17))+ theme_bw()

ggsave("MAF.Diff vs CHRIStoCKDGen.Effects Ratio.png",
       last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

write.csv(tempResult,"ReplicatedSNPs by eGFRw.log.Res.csv")

#-----------------------------------------------------#
#------       Multi-Anc vs. European-Anc       -------
#-----------------------------------------------------#

#Comapring CKDGen Multi-Ancestry Meta-GWAS vs. European Ancestry

repSNPs %>%
  mutate(SNPid = as.character(noquote(str_split(MARKER_ID_CHRIS,
                                                "_",
                                                simplify=TRUE)[,1]))) %>%
  #filter(Locus == "CASZ1") %>% View
  select(SNPid, Locus, EAF_CKDGen, Beta_CKDGen, SE_CKDGen, Pvalue_CKDGen) %>%
  inner_join(repSNPs_old[c("SNPid", "Locus", "Freq1", "Effect", "StdErr", "P.value")],
             by = c("SNPid", "Locus")) %>%
  mutate(Beta_to_SE_MultiAncestry = Effect / StdErr,
         Beta_to_SE_EuropAncestry = Beta_CKDGen / SE_CKDGen) %>%
  ggplot(aes(x = Beta_to_SE_MultiAncestry, y = Beta_to_SE_EuropAncestry, color = Locus)) +
  #ggplot(aes(x = P.value, y = Pvalue_CKDGen, color = Locus)) +
  #ggplot(aes(x = Freq1, y = EAF_CKDGen, color = Locus)) +
  geom_point() +
  geom_abline(slope = 1) + #geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
  theme_classic()

#-------------#

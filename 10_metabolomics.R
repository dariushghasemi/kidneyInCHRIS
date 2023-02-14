#=========================================#
#     CHRIS Metabolites data analysis
#=========================================#

#Download required packages
install.packages("BiocManager")

BiocManager::install(c("Biobase", "AnnotationFilter", "SummarizedExperiment",
"tibble", "dbplyr", "BiocGenerics", "AnnotationDbi",
"S4Vectors", "robustbase", "BiocStyle", "RMariaDB"), force = TRUE)

library(Biobase)
library(AnnotationFilter)
library(SummarizedExperiment)
library(BiocGenerics)
library(AnnotationDbi)
library(S4Vectors)
library(robustbase)
library(BiocStyle)
library(RMariaDB)
library(usethis)
library(git2r)
library(devtools)


remotes::install_git("https://gitlab.gm.eurac.edu/metabolomics/BioCHRIStes", force=T, git = 'external')
#, auth_token = "glpat-UZ6XJp25AxEH97RDz-2E")

remotes::install_git("https://gitlab.gm.eurac.edu/metabolomics/BioCHRIStes",
                     credentials = git2r::cred_user_pass("dghasemisemeskandeh", "azm000n_amari.IR"))

library(BioCHRIStes)
library(biochristes7500)
browseVignettes("biochristes7500")

library(pander)
data(biochristes7500)
do.call(cbind, metadata(biochristes7500))
#-----------------------------------------------------#

#Metabolites data
rowData(biochristes7500)

#Custimized theme
blank_theme <- theme_minimal()+
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

#frequency table and Pie chart
rowData(biochristes7500) %>% View()
  as_tibble() %>% 
  #count(analyte_name) # 175 metabolites
  count(analyte_class) %>% #   6 categories
  #arrange(desc(n))
  ggplot(aes(x="", y = n, fill = fct_inorder(analyte_class))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  scale_fill_brewer("Blues") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5)) + 
  #geom_text(aes(y = n/3 + c(0, cumsum(n)[-length(n)]),
  #              label = scales::percent(n/100)), size=5) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.key.size = unit(1.3, 'cm'),
          legend.key.width = unit(0.9, 'cm'),
          legend.text  = element_text(size = 14))

ggsave("18-11-2022_Metabolites Categories.png", last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
#-----------------------------------------------------#

# individuals characteristics
colData(biochristes7500)

#individual level data
concentrations(biochristes7500)[1:4, 1:3]
concentrations(biochristes7500, blessing = "none")[1:4, 1:3]

# Retrieve name of the metabolites
rowData(biochristes7500)["C14", ]
rowData(biochristes7500) %>%
  as.data.frame() %>% 
  filter(analyte_name == "alpha-AAA") %>% View

# Box plot
boxplot(split(log2(c14), biochristes7500$plate_name), main = "C14",
        ylab = expression(log[2]~signal), las = 2)

pts <- concentrations(biochristes7500, blessing = "none")["Putrescine", ]
plot(hist(log2(pts), breaks = 128), xlab = expression(log[2]~signal),
     main = "Putrescine")
#-----------------------------------------------------#

# Rank or Inverse normal Transformation
# inormal <- function(x){
#   data <<- (qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
# }

# inverse-normal transformation

do_INT = function(values) {
  qnorm((rank(values, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(values)))
}

#preparing metabolites data for jointing with genotypes
chrisMass <- 
  concentrations(biochristes7500,
                 blessing = "none") %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Metabolite") %>%
  pivot_longer(cols = -Metabolite,
               names_to = "AID") %>% 
  pivot_wider(names_from = Metabolite) %>%
  #Caution: log transformation widely affects the results in step 3 -> no mediator after log trans
  mutate(across(!AID, list("ln" = log, "log2" = log2, "int" = INT), .names = "{.col}_{.fn}")) %>% 
  #mutate(across(!AID, ~ INT)) %>% 
  janitor::clean_names() %>% #janitor::row_to_names(1)
  rename(AID = aid)

#-----------------------------------------------------#
#Saving metabolites names for mediation analysis
metabolites <- chrisMass %>% select(-AID) %>% colnames() #%>% View
#targets    <- vcfmod    %>% select(-AID) %>% colnames()

#-----------------------------------------------------#

#Merge metabolites7500 with CHRIS baseline
vcfMass <-
  chris[c("AID", "Age", "Sex", "eGFRw.log.Res")] %>%
  inner_join(PCs_13K,   by = "AID") %>% #dim()
  inner_join(chrisMass, by = "AID") %>% #dim()
  inner_join(vcfmod,    by = "AID") #%>% dim()

#-----------------------------------------------------#

# Show density plot to compare different transformation

pdf('09-Feb-23_Comparison of metabolites distribution under various transformations.pdf', width=8, height = 10)

map(metabolites, function(mtb) {
  chrisMass %>%
  rename_with(~ paste0(.x, "_org"),
              .cols = any_of(metabolites)) %>%
  pivot_longer(cols = -AID,
               names_to = c("Trait", "Function"),
               names_pattern = "(.+)_(org|ln|log2|int)$",
               values_to = "score") %>%
  filter(str_detect(Trait, mtb)) %>%
  ggplot(aes(x = score, fill = Function), color = "white")+
  geom_density(alpha = 0.7)+
  scale_fill_brewer(palette="Dark2")+
  labs(x = mtb)+
  theme_classic()
})

dev.off()

#-----------------------------------------------------#
#-------------- Step 2: trait as covariate -----------
#-----------------------------------------------------#

step2_to_table <- 
  function(mytrait,
           mytarget,
           data){
    results1 <- lapply(mytrait,
                       function(Trait){
                         map_df(mytarget,
                                function(SNP){
                                  myformula <- as.formula(
    eGFRw.log.Res ~ SNP + Trait + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
                                  m <- lm(myformula, data = data)
                                  s <- coef(m)[2]
                                  p <- summary(m)$coefficients[2, c(1,2,4)] #Beta, se, Pvalue
                                  #t <- tidy(m)[2, c(1,2,4)]#broom
                                  return(p)
                                })
                       })
    results2 <- as.data.frame(do.call(cbind, results1)) #%>% clean_names() #library(janitor)
    names(results2) <- str_replace_all(names(results2),
                                       c("Estimate"   = "Estimate",
                                         "Std. Error" = "SE",
                                         "Pr\\([^\\(]*\\)" = "Pvalue"))
    nloc <- length(mytarget)
    results3 <- cbind(SNPid = names(mytarget),
                      Locus = repSNPs$Locus[1:nloc],
                      results2) %>% as_tibble()
    return(results3)
  }

#---------#
MMA_Step2_raw <-
  step2_to_table(vcfMass[metabolites],
                 vcfMass[targets],
                 vcfMass)
#---------#
#Turning the results to longer format for merging with step3
MMA_results_Step2_long <- 
  MMA_Step2_raw %>% 
  pivot_longer(cols = -c(SNPid, Locus),
               names_to = c("Trait", "value"),
               names_pattern = "(.+).(Estimate|SE|Pvalue)$",
               values_to = c("score")) %>%
  pivot_wider(names_from = "value",
              values_from = "score") %>% 
  group_by(SNPid) %>%
  mutate(outlier = ifelse(is_outlier(Estimate),
                          "Yes",
                          "No")) %>% 
  ungroup()
#---------#
write.csv(MMA_results_Step2_long,
          "29-Nov-2022_MMD_Step2_SNPs adjusted for metabolites_long format.csv",
          row.names = F,
          quote = F)

#---------#
#Heatmap for step 2

library(pheatmap)
png("28-Nov-2022_Heatmap_MMD_Step2_SNPs adjusted for metabolites.png", 
    units="in", res = 300, width=10, height=12)
#pdf('pheatmap2.pdf', width=18, height = 18)

# for HeatMap
MMA_results_Step2_long_HM <-
  repSNPs %>%                    
  select(SNPid, Locus, Beta_CHRIS) %>% 
  inner_join(MMA_Step2_raw,
             by = c("SNPid", "Locus")) %>%
  #filter(SNPid %in% leadingSNPs) %>%
  select(SNPid, Locus, Beta_CHRIS, contains("Estimate")) %>%
  mutate(across(contains("Estimate"), function(x) x / Beta_CHRIS))

pheatmap(MMA_results_Step2_long_HM[,-c(1:3)],
         cluster_cols = F,
         cluster_rows = F, 
         show_rownames = T, 
         labels_row = MMA_Step2_raw$SNPid, 
         border_color = NA, 
         fontsize = 4, 
         angle_col = "45")

dev.off()

#-----------------------------------------------------#
#--------------- Step 3: metabolite as outcome -------
#-----------------------------------------------------#

#Mediation Analysis -> Step 3: Metabolites as the outcome

summary(lm(alpha_aaa ~ `chr1:10599281` + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 , data = vcfMass))

#function for iterating retrieving 
#the regression coefficients

step3_to_table <- 
  function(mytrait, 
           mytarget, 
           myformula,
           data){
    Trait <- names(mytrait)
    SNPid <- names(mytarget)
    res3  <- map2_df(
      .x = mytrait,
      .y = myformula,
      .f = function(trait, formula){
        res1 <- map_df(
          .x = mytarget,
          .f = function(SNP)
          {
            m <- lm(as.formula(formula), data = data)
            p <- summary(m)$coefficients[2, c(1,2,4)]
            return(p)
          }
        )

        colnames(res1) <- c("Estimate",
                            "SE",
                            "Pvalue")
        res2 <-
          cbind(SNPid,
                res1) %>%
          merge(
            repSNPs[c("SNPid", "Locus")],
            .,
            by = "SNPid",
            all = F, sort = F)
        
        return(res2)
      }
    )
    
    res3$pheno <- rep(colnames(mytrait),
                      each = length(unique(res3$SNPid)))
    res4 <- 
      res3 %>%
      pivot_wider(id_cols     = c(Locus, SNPid),
                  names_from  = pheno,
                  values_from = c(Estimate, SE , Pvalue),
                  names_glue  = "{pheno}_{.value}")
    return(res4)
  }

#---------#
# Step 3 using vcfMass_nolog can be found more below

MMA_Step3_raw <- 
  step3_to_table(
    vcfMass[metabolites],
    vcfMass[targets],
    paste("trait ~ SNP + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
    vcfMass) %>% 
  select(contains(c("Locus",   #change the order of the columns
                    "SNPid",
                    metabolites))) #%>% View

#---------#
# Reconstruct to longer format
MMA_results_Step3_long <-
  MMA_Step3_raw %>%
  pivot_longer(cols          = -c(SNPid, Locus),
               names_to      = c("Trait", "value"),
               names_pattern = "(.+)_(Estimate|SE|Pvalue)$",
               values_to     = c("score")) %>%
  pivot_wider(names_from     = "value",
              values_from    = "score") %>%
  #Significant association between variants and metabolites
  mutate(associated = ifelse(Pvalue <= 0.05/1925,
                             "Yes",
                             "No")) #%>% filter(associated == "Yes") %>% View

#---------#
write.csv(MMA_results_Step3_long,
          "29-Nov-2022_MMD_Step3_SNPs associated with metabolites.csv",
          row.names = F, quote = F)
#---------#
# results of no-log-transformed metabolites
MMA_results_Step3_long_nolog <-
  step3_to_table(
  vcfMass_nolog[metabolites],
  vcfMass_nolog[targets],
  paste("trait ~ SNP + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"),
  vcfMass_nolog) %>%
  pivot_longer(cols          = -c(SNPid, Locus),
               names_to      = c("Trait", "value"),
               names_pattern = "(.+)_(Estimate|SE|Pvalue)$",
               values_to     = c("score")) %>%
  pivot_wider(names_from     = "value",
              values_from    = "score") %>%
  mutate(associated = ifelse(Pvalue <= 0.05/1925,
                             "Yes",
                             "No"))

#---------#
# Heatmap of step 3
MMA_results_Step3_long %>%
  mutate(SNP = factor(SNPid,
                      levels = str_sort(unique(SNPid),
                                        numeric = TRUE,
                                        decreasing = TRUE))) %>%
  ggplot(aes(x = pheno, y = SNP, fill = related)) +
  geom_tile(color = "grey80") +
  #scale_fill_gradient(low = "red", high = "white")+
  coord_equal(ratio = 1) + 
  #viridis::scale_fill_viridis(name="Estimate")#(option="magma") +
  #ggthemes::theme_tufte(base_family="Helvetica") +
  #theme_classic()+
  theme(legend.title= element_text(size=7, face="bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size=6, face="bold", angle=60, vjust=1.05, hjust=1),
        axis.text.y = element_text(size=6, face="bold"),
        axis.title  = element_blank())

ggsave("29-Nov-22_MMD_Step3_SNPs associated with metabolites_Pvalue.png", 
       last_plot(), width = 16, height = 12, pointsize = 4, dpi = 300, units = "in")

#---------#
#pretty heatmap
pdf('24-Nov-22_Heatmap of MMD_Step3_SNPs associated with metabolites.pdf', width=8, height = 6)
library(pheatmap)
pheatmap(results_Step3_long$Estimate,
         cluster_cols = T, 
         cluster_rows = T, 
         show_rownames = T, 
         labels_row = results_Step3_long$SNPid,
         border_color = NA, #'Black',
         fontsize_row = 6,
         fontsize_col = 6,
         angle_col = "45")

dev.off()

#-----------------------------------------------------#
#-------------------- Merge 3 Steps ------------------
#-----------------------------------------------------#

#Testing significance of adjusted effect in step 2 for mediator variable in step 3 and step 2
Sobel_test <-
  function(beta_step2, se_step2, beta_step3, se_step3){
    a    = beta_step3
    sa   = se_step3
    b    = beta_step2
    sb   = se_step2
    tmp1 = b^2 * sa^2 + a^2 * sb^2
    zsob = a * b/sqrt(tmp1)
    psob = pnorm(-abs(zsob)) * 2
    p.value = psob
    z.value = zsob
    out = data.frame(cbind(z.value, p.value))
    #names(out) = c("Sobel")
    return(out)
  }


#-----------------------------------------------------#
#-------------------- Merge 3 Steps ------------------
#-----------------------------------------------------#

#MMA_sum3steps_long <-
MMA_results_Step2_long %>%
  inner_join(
    MMA_results_Step3_long,
    by = c("SNPid" = "SNPid",
           "Locus" = "Locus",
           "Trait" = "Trait"),
    suffix = c("_Step2",
               "_Step3")) %>%
  inner_join(
    repSNPs[c("SNPid",
              "Locus",
              "RA_CHRIS_disc",
              "EA_CHRIS_disc",
              "Beta_CHRIS",
              "SE_CHRIS",
              "Pvalue_CHRIS")],
    by = c("SNPid",
           "Locus")) %>%
  mutate(
    EA_OA = paste0(EA_CHRIS_disc, "/", RA_CHRIS_disc),
    outlierRelated = case_when(
    outlier == "Yes" & associated == "Yes" ~ "Yes/Yes",
    outlier == "Yes" & associated != "Yes" ~ "Yes/No",
    outlier != "Yes" & associated == "Yes" ~ "No/Yes",
    outlier != "Yes" & associated != "Yes" ~ "No/No"),
    mediator = ifelse(outlier == "Yes" & associated == "Yes", "Yes", "No"),
    Sobel = Sobel_test(Estimate_Step2, SE_Step2, Estimate_Step3, SE_Step3),
    SNP = factor(SNPid,
                 levels = str_sort(unique(SNPid),
                                   numeric = TRUE,
                                   decreasing = TRUE))) %>%
  rename(
    Estimate_GWAS = Beta_CHRIS,
    SE_GWAS       = SE_CHRIS,
    Pvalue_GWAS   = Pvalue_CHRIS) %>%
  select(
    Locus, SNPid, EA_OA,
    Estimate_GWAS, SE_GWAS, Pvalue_GWAS, everything()) %>%
  filter(
    outlier == "Yes"
    #Trait == "alpha_aaa",
    #Locus == "CASZ1",
    #Pvalue_Step3 < 0.05/11/175,
    #mediator == "Yes"
    ) %>% View
  count(Locus, outlier) %>% View
  #write.csv(., "29-Nov-22_Heatmap_MMD_Step 1&2&3_outlierRelated traits_Mediatory metabolites.csv", row.names = FALSE)
  ggplot(aes(x = trait, y = SNP, fill = outlierRelated)) +
    geom_tile() +
    theme_classic() +
  scale_fill_manual(values=c('white', "grey50", "#FF6666"))+ #'#999999', "#E69F00"
  labs(x    = "",
       y    = "")+
  theme(legend.title = element_text(size=7, face="bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size=4, face="bold", angle=90, vjust=1.05, hjust=1),
        axis.text.y = element_text(size=4, face="bold"),
        axis.title  = element_text(size=12,face="bold"))

ggsave("29-Nov-22_Heatmap_MMD_Step 1&2&3_outlierRelated traits.png", 
       last_plot(), width = 10, height = 9, pointsize = 4, dpi = 300, units = "in")

#---------#
# compare step 3 results for log-transformed vs non-transformed metabolites 
MMA_results_Step3_long_nolog %>%
  inner_join(MMA_results_Step3_long,
             by = c("Locus", "SNPid", "Trait"),
             suffix = c("_nonTrans", "_logTrans")) %>%
  filter(Trait == "alpha_aaa") %>% 
  ggplot(aes(x = Estimate_nonTrans,
             # xmin =  Estimate_nonTrans - SE_nonTrans,
             # xmax =  Estimate_nonTrans + SE_nonTrans,
             y = Estimate_logTrans,
             # ymin =  Estimate_logTrans - SE_logTrans,
             # ymax =  Estimate_logTrans + SE_logTrans,
             color = Locus, shape = Locus))+
  geom_abline(slope = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 0.9, size = 2.5, fatten = 2.2) +
  scale_shape_manual(values = c(19, 17, 18, 17, 19, 15, 17, 19, 17, 19, 15)) +
  scale_color_manual(values = c(
    "maroon1", "darkorchid2","orange2", "green4", "steelblue2", "darkturquoise", "tomato",
    "springgreen2", "royalblue2", "gold", "grey50")) +
  labs(x = "SNPs effect on ln(eGFRcreat) adj for metabolite Alpha_AAA",
       y = "SNPs effect on ln(eGFRcreat) adj for metabolite ln(Alpha_AAA)")+
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.placement = "outside",
        axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=8, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.key.size  = unit(0.99, 'cm'),
        legend.key.width = unit(0.7, 'cm'),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))

ggsave("29-Jan-23_SNPs effect adj for Alpha_AAA metabolite with_out log transformation.png",
       last_plot(), width = 10, height = 7, pointsize = 4, dpi = 300, units = "in")


#---------#
# Checking the distribution
vcfMass %>%
  select(AID, eGFRw.log.Res, alpha_aaa, starts_with("chr1:")) %>%
  mutate(alpha_aaa_log = log(2 + alpha_aaa)) %>%
    # Dosage_SLC34A1 = cut(`chr5:177386403`,
    #                           breaks = c(-Inf, 0.500, 1.500, Inf),
    #                           labels = c("0", "1", "2"))) %>%
  pivot_longer(cols = c(alpha_aaa, alpha_aaa_log),
               names_to  = "metabolite",
               values_to = "value") %>%
  #ggplot(aes(Dosage_SLC34A1, eGFRw.log.Res))+
  #geom_violin(aes(fill = Peptited), position = "dodge")+
  #geom_boxplot(aes(fill = Peptited), width = 0.3, position = position_dodge(.9))+
  ggplot(aes(x = value)) +
  geom_density(aes(fill = metabolite), alpha = 0.6, color = "grey50") + #, show.legend = NULL
  #scale_fill_manual(values = c("violetred1", "turquoise2"), labels = c("No", "Yes")) +
  theme_classic()


ggsave("29-Jan-23_Density plot of Alpha_AAA metabolite with_out log transformation.png", last_plot(), width = 10, height = 7, pointsize = 4, dpi = 300, units = "in")

#---------#



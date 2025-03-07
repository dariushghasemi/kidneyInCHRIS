
# libraries
library(interactions) # for interaction plot in paper
library(jtools)       # for interaction plot in paper
library(cowplot)
source("~/projects/kidneyInCHRIS/scripts/08-1_mediation_data.R") #to load vcfReg_TSHmod


# inputs


# outputs
out_plt_interaction <- "~/projects/kidneyInCHRIS/outputs/07-Mar-25_SNP-TSH interactions rs819196 and rs819185 at STC1.jpg"


# convert quantitative variables from character to numeric

vcfReg_TSHmod <- vcfReg_TSHmod %>% 
  dplyr::mutate(
    across(Age, as.numeric), # convert AGE to numeric
    #TSH_cen = scale(TSH, scale = F), # cetralize TSH
    #eGFRw.log_cen = scale(eGFRw.log, scale = F)
    )


#-----------------------------------------------------#
#-------        Depicting interaction          -------
#-----------------------------------------------------#


# Visualizing SNP-TSH interaction
lapply(
  c("chr8:23885208", "chr8:23894869"), 
  function(SNP){
    
    snp_tidied <- str_replace(SNP, ":", "_")
    x_axis_label <- paste0("15-Feb-23_SNP-TSH interaction effect on eGFRw.log for ", snp_tidied, " at STC1.png")
    
    # prepapre data and plot
    vcfReg_TSHmod %>%
      dplyr::select(eGFRw.log, SNP, TSH, Age, Sex, starts_with("PC")) %>% #starts_with("chr8:")
      filter(!is.na(TSH_cat)) %>%
      dplyr::mutate(
        #TSH_cen = scale(TSH, scale = F),
        TSH_third = case_when(
          TSH < quantile(TSH, 1/3, na.rm = T) ~ "Hyperthyrodism",
          TSH > quantile(TSH, 2/3, na.rm = T) ~ "Hypothyrodism",
          TRUE ~ "Normal TSH"
          ),
        across(c(Age, starts_with("PC")), ~ scale(., scale = F))
        ) %>%
      ggplot(aes(x = !!sym(SNP), y = eGFRw.log, color = TSH_third)) +
      geom_smooth(method = "lm", se = T) + 
      #facet_wrap(~TSH_cat, nrow = 3) +
      scale_color_manual(
        name = "Stratified TSH",
        values = c("steelblue3", "green4", "turquoise2"),
        breaks = c("Hyperthyrodism", "Normal TSH", "Hypothyrodism")) +
      labs(
        x = paste("Dosage of", SNP, "variant at STC1"),
        y = "log(eGFRcreat)"
        ) +
      theme_classic()
    
    # Save the plot
    ggsave(filename = x_axis_label, plot = last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")
    }
  )


#-----------------------------------------------------#
#-------        Fig 4: Interaction plot        -------
#-----------------------------------------------------#


# function to draw theme of the plot
theme_interction <- function(...){
  theme_classic() +
    #jtools::theme_apa()
    theme(
      panel.background = element_rect(fill = "white"),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 12, face = "bold"),
      strip.placement = "outside",
      axis.text    = element_text(size = 12,  face = "bold"),
      axis.title   = element_text(size = 14, face = "bold"),
      axis.ticks = element_line(linewidth = .8),
      axis.ticks.length = unit(2, 'mm'),
      #axis.title.x = ggtext::element_markdown(),
      legend.key.size  = unit(0.8, 'cm'),
      legend.key.width = unit(0.5, 'cm'),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold")
    )
}


#---------#
# Forming principal components term of the model 
PCs <- paste0("PC", 1:10, collapse = " + ")

#---------#
# Emerging SNP-TSH interaction via package
# rs819196 = `chr8:23885208` 
# rs819185 = `chr8:23894869`

# first run model with interaction for each snp
fit_8_23885208 <- lm(paste0("eGFRw.log ~ `chr8:23885208` * TSH + Sex + Age +", PCs), data = vcfReg_TSHmod)
fit_8_23894869 <- lm(paste0("eGFRw.log ~ `chr8:23894869` * TSH + Sex + Age +", PCs), data = vcfReg_TSHmod)

#---------#
# plots

plt_8_23885208 <- interactions::interact_plot(
  model = fit_8_23885208,
  pred = `chr8:23885208`,
  modx = TSH,
  modx.labels = c("Low", "Average", "High"),
  linearity.check = F,
  plot.points = F,
  line.thickness = 2,
  colors = "seagreen",
  #x.label = "Dosage  for  rs819196  at  *STC1*  locus",
  y.label = "ln(eGFRcrea)"
  ) +
  labs(
    x = "\nDosage of rs819196",
    y = "ln(eGFRcrea)\n"
  ) +
  theme_interction()+
  theme(legend.position = c(.2, .8))



plt_8_23894869 <- interactions::interact_plot(
  model = fit_8_23894869,
  pred = `chr8:23894869`,
  modx = TSH,
  #interval = T,
  linearity.check = F,
  plot.points = F,
  line.thickness = 2,
  colors = "seagreen",
  #legend.main = "Custom Legend Title",
  #x.label = "rs819185       dosage",
  y.label = "ln(eGFRcrea)"
) +
  labs(x = "\nDosage of rs819185") +
  theme_interction() +
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  )


#---------#
# Join the plots
plt_joint <- cowplot::plot_grid(plt_8_23885208, plt_8_23894869, ncol = 2)

# save joind plot
ggsave(plt_joint, filename = out_plt_interaction, width = 11, height = 6, dpi = 300, units = "in")





# libraries
library(interactions) # for interaction plot in paper
library(jtools)       # for interaction plot in paper
library(gridExtra)


# inputs


# outputs
out_plt_interaction <- "22-Mar-23_SNP-TSH interactions rs819196 and rs819185 at STC1.png"



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


#---------#
# Emerging SNP-TSH interaction via package `chr8:23885208` `chr8:23894869`
int_SNP2 <- interactions::interact_plot(
    lm(paste0("eGFRw.log ~ `chr8:23894869` * TSH + Sex + Age +", PCs), vcfReg_TSHmod),
    pred = `chr8:23894869`,
    modx = TSH,
    interval = T,
    linearity.check = F,
    plot.points = F,
    colors = "seagreen",
    #legend.main = "Custom Legend Title",
    x.label = "Dosage for rs819185 at *STC1* locus",
    y.label = "ln(eGFRcrea)"
    ) +
  theme_classic() +
  #jtools::theme_apa()
  theme(
    panel.background = element_rect(fill = "white"),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.placement = "outside",
    axis.text    = element_text(size = 8,  face = "bold"),
    axis.title   = element_text(size = 12, face = "bold"),
    axis.title.x = ggtext::element_markdown(),
    legend.position = c(.8, .22),
    legend.key.size  = unit(0.8, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
    )


# 1st way
gridExtra::grid.arrange(int_SNP1, int_SNP2, ncol = 2, labels = c("A", "B"))

# 2nd way
cowplot::plot_grid(int_SNP1, int_SNP2, labels = c("A", "B"), ncol = 2, nrow = 1)

# save joind plot
ggsave(filename = out_plt_interaction, width = 11, height = 5, pointsize = 5, dpi = 300, units = "in")




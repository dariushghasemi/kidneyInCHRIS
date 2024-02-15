
data <- matrix(rnorm(600000), ncol = 6000)

# Run PCA
pca_result <- prcomp(data, scale. = TRUE)

# Calculate cumulative variance explained
cumulative_variance <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
cumulative_variance

# Create a data frame for ggplot
df <- data.frame(Components = 1:length(cumulative_variance), 
                 Cumulative_Variance = cumulative_variance)

# Plot the cumulative variance explained using ggplot
ggplot(df, aes(x = Components, y = Cumulative_Variance)) +
  geom_line() +
  geom_point(color = "steelblue", size = 3) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(0, 100, 200, 500, 1000, 2000, 4000, 6000)) +
  scale_y_continuous(breaks = c(0, .25, .50, .75, .95, 0.99, 1)) +
  labs(x = "Number of Components", y = "Cumulative Variance Explained") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 6))

factoextra::fviz_eig(pca_result, addlabels = F, ncp = 100, ylim = c(0, 15)) + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  theme(title = element_blank())


# compute variance
repSNPs_EA %>%
  mutate(var_CHRIS = (Beta_CHRIS^2)*((2*MAF_CHRIS*(1-MAF_CHRIS))/0.016),
         var_CKDGen = (Beta_CKDGen^2)*((2*EAF_CKDGen*(1-EAF_CKDGen))/0.016))

repSNPs_EA %>%
  ggplot(aes(Beta_CHRIS, Beta_CKDGen)) +
  geom_point()


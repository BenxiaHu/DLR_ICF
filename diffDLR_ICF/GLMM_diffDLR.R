# --------------------------------------------
# Step 0: Install and Load Required Packages
# --------------------------------------------

# Install and load the lme4 package
if (!requireNamespace("lme4", quietly = TRUE)) {
  install.packages("lme4")
}
library(lme4)  # For fitting generalized linear mixed models

# Install and load the lmerTest package
if (!requireNamespace("lmerTest", quietly = TRUE)) {
  install.packages("lmerTest")
}
library(lmerTest)  # Provides p-values for fixed effects

# Install and load the emmeans package
if (!requireNamespace("emmeans", quietly = TRUE)) {
  install.packages("emmeans")
}
library(emmeans)  # For estimated marginal means and contrasts

# Install and load the multcomp package (for multiple comparisons and FDR adjustment)
if (!requireNamespace("multcomp", quietly = TRUE)) {
  install.packages("multcomp")
}
library(multcomp)

# Install and load ggplot2 for plotting (optional)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# --------------------------------------------
# Step 1: Load and Prepare Your Data
# --------------------------------------------

# Load your data (replace 'your_data.csv' with the path to your data file)
data <- read.csv("your_data.csv")

# Ensure that Treatment and Site are factors
data$Treatment <- as.factor(data$Treatment)
data$Site <- as.factor(data$Site)

# View the first few rows of your data
head(data)

# Check the structure of the data
str(data)

# --------------------------------------------
# Step 2: Fit the Generalized Linear Mixed Model (GLMM)
# --------------------------------------------

# Fit the GLMM using glmer from the lme4 package
# Fixed Effects: Treatment (effect of treatment vs. control)
# Random Effects: (Treatment | Site) allows for random intercepts and slopes by Site
# Family: Gaussian (normal distribution) with identity link function
model <- glmer(DLR ~ Treatment + (Treatment | Site), data = data, family = gaussian(link = "identity"))

# Check for convergence warnings
if (length(model@optinfo$conv$lme4$messages) == 0) {
  cat("Model converged successfully.\n")
} else {
  cat("Convergence warnings:\n")
  print(model@optinfo$conv$lme4$messages)
}

# --------------------------------------------
# Step 3: Obtain p-values for Each Site
# --------------------------------------------

# Use emmeans to get estimated marginal means for Treatment within each Site
emm <- emmeans(model, ~ Treatment | Site)

# Perform pairwise comparisons (Treatment vs. Control) within each Site
contrast_results <- contrast(emm, method = "pairwise", adjust = "fdr")

# Convert the results to a data frame
contrast_df <- as.data.frame(contrast_results)

# View the first few rows of the results
head(contrast_df)

# The contrast_df data frame contains:
# - Site: the site identifier
# - contrast: the comparison made (Treatment - Control)
# - estimate: estimated difference between Treatment and Control
# - SE: standard error of the estimate
# - df: degrees of freedom
# - t.ratio: t-statistic
# - p.value: p-value (adjusted for multiple comparisons using FDR)

# --------------------------------------------
# Step 4: Extract and Adjust p-values for Each Site
# --------------------------------------------

# The p-values are already adjusted using FDR in the contrast_results
# If you want to adjust p-values manually, you can do so as follows:
# contrast_df$p.value <- p.adjust(contrast_df$p.value, method = "fdr")

# Rename columns for clarity
site_results <- contrast_df[, c("Site", "estimate", "SE", "df", "t.ratio", "p.value")]
names(site_results) <- c("Site", "Estimate", "Std_Error", "df", "t_value", "p_value")

# View the first few results
head(site_results)

# --------------------------------------------
# Step 5: Interpret the Results
# --------------------------------------------

# Sites where p_value < 0.05 are considered to have a significant treatment effect
significant_sites <- site_results[site_results$p_value < 0.05, ]

# View significant sites
print("Significant Sites (p < 0.05):")
print(significant_sites)

# --------------------------------------------
# Step 6: Save or Export the Results
# --------------------------------------------

# Save the site_results data frame to a CSV file
write.csv(site_results, "site_specific_p_values.csv", row.names = FALSE)

# --------------------------------------------
# Step 7: Additional Diagnostics (Optional)
# --------------------------------------------

# Plot the estimated treatment effects by site
ggplot(site_results, aes(x = as.numeric(Site), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - 1.96 * Std_Error, ymax = Estimate + 1.96 * Std_Error), width = 0.2) +
  labs(title = "Site-Specific Treatment Effects",
       x = "Site",
       y = "Estimated Treatment Effect") +
  theme_minimal()

# --------------------------------------------
# Step 8: Model Diagnostics
# --------------------------------------------

# Residuals vs. Fitted plot
plot(fitted(model), resid(model),
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red", lty = 2)

# Q-Q plot of residuals
qqnorm(resid(model),
       main = "Normal Q-Q Plot of Residuals")
qqline(resid(model), col = "red", lty = 2)

# --------------------------------------------
# End of Updated Code
# --------------------------------------------

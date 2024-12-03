# --------------------------------------------
# Step 0: Install and Load Required Packages
# --------------------------------------------

# Install and load the lme4 package
if (!requireNamespace("lme4", quietly = TRUE)) {
  install.packages("lme4")
}
library(lme4)  # Used for fitting linear mixed models

# Install and load the lmerTest package
if (!requireNamespace("lmerTest", quietly = TRUE)) {
  install.packages("lmerTest")
}
library(lmerTest)  # Provides p-values for fixed effects using Satterthwaite's method

# Install and load the multcomp package (for multiple comparisons and FDR adjustment)
if (!requireNamespace("multcomp", quietly = TRUE)) {
  install.packages("multcomp")
}
library(multcomp)

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
# Step 2: Fit the Linear Mixed Model (LMM)
# --------------------------------------------

# Fit the LMM using lmer from the lmerTest package
# Fixed Effects: Treatment (effect of treatment vs. control)
# Random Effects: (1 | Site) models random intercepts for each site
model <- lmer(DLR ~ Treatment + (Treatment | Site), data = data)

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

# Adjust p-values using FDR (already adjusted in contrast_results)
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
library(ggplot2)

ggplot(site_results, aes(x = as.numeric(Site), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - 1.96 * Std_Error, ymax = Estimate + 1.96 * Std_Error), width = 0.2) +
  labs(title = "Site-Specific Treatment Effects",
       x = "Site",
       y = "Estimated Treatment Effect") +
  theme_minimal()

# --------------------------------------------
# End of Updated Code
# --------------------------------------------
library(glmmTMB)
library(dplyr)
library(stats)
library(multcomp)
library(ggplot2)
library(tidyr)

set.seed(42)  # For reproducibility

# 2. Simulate Example Data
n_bins <- 2000
# True means for D and L
mu_D_true <- 2000
mu_L_true <- 2000

# Dispersion parameters (theta) for NB
theta_D <- 10
theta_L <- 20

# Function to simulate NB counts
simulate_nb <- function(mu, theta, size) {
  # In glmmTMB, 'size' corresponds to the dispersion parameter
  rnbinom(n_bins, size = theta, mu = mu)
}

# Simulate counts
D_counts <- simulate_nb(mu_D_true, theta_D, n_bins)
L_counts <- simulate_nb(mu_L_true, theta_L, n_bins)

# Create DataFrame
data <- data.frame(
  Bin = paste0("Bin_", 1:n_bins),
  D = D_counts,
  L = L_counts
)

# View first few rows
head(data)

# 3. Reshape Data for Modeling
# Convert to long format
data_long <- data %>%
  pivot_longer(cols = c(D, L), names_to = "Interaction", values_to = "Count")

# Add a column to indicate the interaction type as a factor
data_long$Interaction <- factor(data_long$Interaction, levels = c("L", "D"))

# View the transformed data
head(data_long)

# 4. Define the Model
## a Generalized Linear Mixed Model (GLMM) with an NB distribution to account for overdispersion. The model will include:
## Fixed Effect: Interaction type (D vs. L).
## Random Effect: Bin identifier to account for the pairing within bins.

# Fit the GLMM
model <- glmmTMB(
  Count ~ Interaction + (1 | Bin),
  family = nbinom2,  # Negative binomial family with variance = mu + mu^2 / theta
  data = data_long
)

# Summarize the model
summary(model)

# 5. Perform Hypothesis Testing
# a. Calculate DLRs and Their Confidence Intervals
# Extract fixed effect coefficients
fixed_effects <- fixef(model)$cond
interaction_coef <- fixed_effects["InteractionD"]

# Calculate the log D/L ratio
# log(D/L) = InteractionD coefficient
# Exponentiate to get D/L ratio on the original scale
global_DLR <- exp(interaction_coef)
print(paste("Global D/L Ratio:", round(global_DLR, 2)))


# b. Estimate Bin-Specific DLRs
# Extract random effects (bin-specific deviations)
random_effects <- ranef(model)$Bin

# Compute bin-specific log(D/L) = InteractionD + bin deviation
# Since the model is Count ~ Interaction + (1 | Bin),
# and Interaction is a fixed effect, bin deviations capture deviations from the fixed effect

# Merge random effects with bin identifiers
random_effects_df <- data.frame(
  Bin = rownames(random_effects),
  deviation = random_effects[, "Intercept"]
)

# Merge with the original data
data_with_deviation <- data_long %>%
  pivot_wider(names_from = Interaction, values_from = Count) %>%
  left_join(random_effects_df, by = "Bin") %>%
  mutate(
    log_DLR = log(D / L) + deviation  # Adjust log(D/L) with bin deviation
  ) %>%
  mutate(
    DLR = D / L
  )

# View first few rows
head(data_with_deviation)


# c. Define Thresholds for Significant DLR Peaks
# Calculate standard error of the global DLR
# The variance of the fixed effect can be extracted from the model summary
model_summary <- summary(model)
se_interaction <- coef(model_summary)$cond["InteractionD", "Std. Error"]

# Define confidence interval (e.g., 95%)
alpha <- 0.05
z_score <- qnorm(1 - alpha / 2)

# Compute lower and upper bounds on the log scale
lower_log <- interaction_coef - z_score * se_interaction
upper_log <- interaction_coef + z_score * se_interaction

# Exponentiate to get bounds on the D/L ratio scale
lower_DLR <- exp(lower_log)
upper_DLR <- exp(upper_log)

print(paste("95% CI for Global D/L Ratio: [", round(lower_DLR, 2), ",", round(upper_DLR, 2), "]", sep = ""))

# d. Classify Each Bin Based on DLR
# Define classification based on DLR thresholds
data_with_deviation <- data_with_deviation %>%
  mutate(
    DLR_classification = case_when(
      DLR < lower_DLR ~ "Local DLR",
      DLR > upper_DLR ~ "Distal DLR",
      TRUE ~ "Neutral"
    )
  )

# View classifications
head(data_with_deviation)


# e. Visualize the Results
# Plot DLR classifications
ggplot(data_with_deviation, aes(x = Bin, y = DLR, color = DLR_classification)) +
  geom_point() +
  geom_hline(yintercept = lower_DLR, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = upper_DLR, linetype = "dashed", color = "blue") +
  theme_bw() +
  labs(title = "Distal-to-Local Ratio (DLR) per Genomic Bin",
       x = "Genomic Bin",
       y = "DLR (D/L)",
       color = "DLR Classification") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Extract bin-specific deviations
bin_effects <- ranef(model)$Bin$Intercept

# Calculate standard error for the fixed effect
se_fixed <- coef(summary(model))$cond["InteractionD", "Std. Error"]

# Calculate bin-specific log(D/L)
bin_log_DLR <- fixed_effects["InteractionD"] + bin_effects

# Calculate standard error for bin-specific log(D/L)
# Assuming independence between fixed effect and random effect
# The variance of bin_log_DLR is se_fixed^2 + Var(random_effects)
# Estimate Var(random_effects) from model
var_bin <- VarCorr(model)$Bin$cond["Intercept"]

se_bin <- sqrt(se_fixed^2 + var_bin)

# Calculate z-scores
z_scores <- bin_log_DLR / se_bin

# Calculate two-tailed p-values
p_values <- 2 * (1 - pnorm(abs(z_scores)))

# Add p-values to the data
data_with_deviation$p_value <- p_values

# Adjust p-values using BH method
data_with_deviation$adj_p_value <- p.adjust(data_with_deviation$p_value, method = "BH")

# Classify based on adjusted p-values
data_with_deviation <- data_with_deviation %>%
  mutate(
    DLR_classification_adj = case_when(
      adj_p_value < 0.05 & DLR > upper_DLR ~ "Distal DLR",
      adj_p_value < 0.05 & DLR < lower_DLR ~ "Local DLR",
      TRUE ~ "Neutral"
    )
  )

# View significant peaks
significant_peaks <- data_with_deviation %>%
  filter(DLR_classification_adj != "Neutral")

print(significant_peaks[, c("Bin", "D", "L", "DLR", "DLR_classification_adj")])

# 6. Multiple Testing Correction
# a. Calculate p-values for Each Bin's DLR
# Extract bin-specific deviations
bin_effects <- ranef(model)$Bin$Intercept

# Calculate standard error for the fixed effect
se_fixed <- coef(summary(model))$cond["InteractionD", "Std. Error"]

# Calculate bin-specific log(D/L)
bin_log_DLR <- fixed_effects["InteractionD"] + bin_effects

# Calculate standard error for bin-specific log(D/L)
# Assuming independence between fixed effect and random effect
# The variance of bin_log_DLR is se_fixed^2 + Var(random_effects)
# Estimate Var(random_effects) from model
var_bin <- VarCorr(model)$Bin$cond["Intercept"]

se_bin <- sqrt(se_fixed^2 + var_bin)

# Calculate z-scores
z_scores <- bin_log_DLR / se_bin

# Calculate two-tailed p-values
p_values <- 2 * (1 - pnorm(abs(z_scores)))

# Add p-values to the data
data_with_deviation$p_value <- p_values

# Adjust p-values using BH method
data_with_deviation$adj_p_value <- p.adjust(data_with_deviation$p_value, method = "BH")

# Classify based on adjusted p-values
data_with_deviation <- data_with_deviation %>%
  mutate(
    DLR_classification_adj = case_when(
      adj_p_value < 0.05 & DLR > upper_DLR ~ "Distal DLR",
      adj_p_value < 0.05 & DLR < lower_DLR ~ "Local DLR",
      TRUE ~ "Neutral"
    )
  )

# View significant peaks
significant_peaks <- data_with_deviation %>%
  filter(DLR_classification_adj != "Neutral")

print(significant_peaks[, c("Bin", "D", "L", "DLR", "DLR_classification_adj")])

# b. Summary of Significant DLR Peaks
# Count of each classification
data_with_deviation %>%
  group_by(DLR_classification_adj) %>%
  summarise(Count = n())












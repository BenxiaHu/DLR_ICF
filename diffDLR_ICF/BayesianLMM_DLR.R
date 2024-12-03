# Install brms if not already installed
if (!require(brms)) {
  install.packages("brms")
  library(brms)
}

# Load other necessary packages
library(ggplot2)
library(dplyr)

# Load your data (replace 'your_data.csv' with your actual data file)
data <- read.csv("your_data.csv")

# Ensure that Treatment and Site are factors
data$Treatment <- as.factor(data$Treatment)
data$Site <- as.factor(data$Site)

# View the first few rows
head(data)

# Summary statistics
summary(data)

# Explore the distribution of DLR
ggplot(data, aes(x = DLR)) +
  geom_histogram(binwidth = 0.5, fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of DLR Values", x = "DLR", y = "Frequency")

## 2.3. Estimating Hyperparameters from the Data
### 2.3.1. Estimating the Variance Components

#### Use the method of moments or a simple linear model to get initial estimates.
# Calculate overall mean and variance of DLR
mean_DLR <- mean(data$DLR)
var_DLR <- var(data$DLR)

# Estimate within-site variance
within_site_var <- data %>%
  group_by(Site) %>%
  summarize(var_within = var(DLR)) %>%
  pull(var_within) %>%
  mean(na.rm = TRUE)

# Estimate between-site variance
between_site_var <- var_DLR - within_site_var

# Within-Site Variance: Average variance of DLR within each site.
# Between-Site Variance: Difference between total variance and within-site variance.

### 2.3.2. Setting Priors Based on Estimates
# Set priors using estimated variances
beta_prior_sd <- sqrt(var_DLR)       # For fixed effects
sigma_u_prior_sd <- sqrt(between_site_var)  # For random effects
sigma_prior_sd <- sqrt(within_site_var)     # For residual variance

# Define priors
priors <- c(
  set_prior(paste0("normal(0, ", beta_prior_sd, ")"), class = "b"),
  set_prior(paste0("normal(0, ", sigma_u_prior_sd, ")"), class = "sd", group = "Site", coef = "Intercept", lb = 0),
  set_prior(paste0("normal(0, ", sigma_prior_sd, ")"), class = "sigma", lb = 0)
)

### Priors Based on Data: Use the estimated variances to set the standard deviations of the priors.
### Lower Bounds (lb = 0): Ensure standard deviations are positive.

### 2.4. Fitting the Bayesian LMM

# Define the model formula
formula <- bf(DLR ~ Treatment + (1 | Site))

# Fit the Bayesian LMM
bayes_lmm <- brm(
  formula = formula,
  data = data,
  family = gaussian(),
  prior = priors,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 123
)

### Model Specification: Treatment as a fixed effect, Site as a random effect.
### Iterations and Warmup: Sufficient iterations for convergence.
### Cores: Parallel computation for efficiency.

### 2.5. Summarizing and Interpreting the Results

# Summarize the model
summary(bayes_lmm)

##### Fixed Effects (Population-Level Effects): Estimates of the intercept and treatment effect.
##### Random Effects (Group-Level Effects): Estimates of the standard deviation for Site (sd(Site)(Intercept)).
##### Residual Standard Deviation (sigma): Estimate of residual variability.

#### Extracting Hyperparameters:
# Extract fixed effects estimates
fixed_effects <- fixef(bayes_lmm)

# Extract random effects standard deviation
random_effects_sd <- VarCorr(bayes_lmm)$Site$sd

# Extract residual standard deviation
residual_sd <- sigma(bayes_lmm)

# Print estimates
print("Fixed Effects:")
print(fixed_effects)

print("Random Effects Standard Deviation (Site):")
print(random_effects_sd)

print("Residual Standard Deviation:")
print(residual_sd)

#### fixef(bayes_lmm): Extracts fixed effects estimates (β).
#### VarCorr(bayes_lmm): Provides variance components, including σu
#### sigma(bayes_lmm): Residual standard deviation (σ).

#### 2.6. Detecting Site-Specific Effects
#### To detect significant treatment effects at each site, we'll extend the model to include random slopes.

#### 2.6.1. Modifying the Model to Include Random Slopes
# Define the model with random slopes
formula_random_slopes <- bf(DLR ~ Treatment + (Treatment | Site))

# Fit the model
bayes_lmm_random_slopes <- brm(
  formula = formula_random_slopes,
  data = data,
  family = gaussian(),
  prior = priors,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = 4,
  seed = 123
)

### Random Slopes: Allows the treatment effect to vary by site.
### Increased Iterations: More iterations may be necessary for model convergence due to added complexity.

#### 2.6.2. Extracting Site-Specific Treatment Effects
# Extract site-specific coefficients
site_effects <- coef(bayes_lmm_random_slopes)$Site

# Create a data frame of site-specific treatment effects
site_treatment_effects <- data.frame(
  Site = rownames(site_effects),
  Intercept = site_effects[, "Estimate", "Intercept"],
  Treatment_Effect = site_effects[, "Estimate", "TreatmentTreatment"],
  Treatment_Effect_SE = site_effects[, "Est.Error", "TreatmentTreatment"],
  Treatment_Effect_lower = site_effects[, "Q2.5", "TreatmentTreatment"],
  Treatment_Effect_upper = site_effects[, "Q97.5", "TreatmentTreatment"]
)

# View the first few site-specific effects
head(site_treatment_effects)

### coef(bayes_lmm_random_slopes): Extracts estimated coefficients for each site.
### Credible Intervals: Provided by Q2.5 and Q97.5 for 95% credible intervals.

#### 2.6.3. Identifying Significant Sites
# Determine which sites have significant treatment effects
site_treatment_effects$Significant <- with(site_treatment_effects, 
                                           ifelse(Treatment_Effect_lower > 0 | Treatment_Effect_upper < 0, "Yes", "No"))

# View significant sites
significant_sites <- site_treatment_effects[site_treatment_effects$Significant == "Yes", ]
print(significant_sites)

### Significance Criteria: A site is significant if the 95% credible interval does not include zero.
### Significant Column: Indicates whether the treatment effect at each site is significant.

#### 2.7. Visualizing Site-Specific Effects
# Plot site-specific treatment effects with credible intervals
ggplot(site_treatment_effects, aes(x = reorder(Site, Treatment_Effect), y = Treatment_Effect)) +
  geom_point() +
  geom_errorbar(aes(ymin = Treatment_Effect_lower, ymax = Treatment_Effect_upper), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Site-Specific Treatment Effects", x = "Site", y = "Treatment Effect") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

#### Visualization: Helps identify sites with significant positive or negative treatment effects.
#### Horizontal Line at Zero: Indicates no treatment effect; sites with intervals not crossing zero are significant.

#### 2.8. Model Diagnostics
# Plot diagnostics
plot(bayes_lmm_random_slopes)

# Check R-hat values
rhat_values <- rhat(bayes_lmm_random_slopes)
print(rhat_values)

# Posterior predictive checks
pp_check(bayes_lmm_random_slopes)

#### Trace Plots and Density Plots: Assess convergence and mixing of chains.
#### R-hat Values: Values close to 1 indicate convergence.
#### Posterior Predictive Checks: Evaluate model fit.

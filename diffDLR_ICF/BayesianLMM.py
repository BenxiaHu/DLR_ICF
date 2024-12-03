# Install packages if necessary
# !pip install pymc3 arviz pandas numpy matplotlib

import pymc3 as pm
import arviz as az
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#  3.2. Loading and Exploring Your DLR Data
# Load your data
data = pd.read_csv("your_data.csv")

# Ensure that Treatment and Site are categorical variables
data['Treatment'] = data['Treatment'].astype('category')
data['Site'] = data['Site'].astype('category')

# Encode categorical variables
data['Treatment_code'] = data['Treatment'].cat.codes  # Control=0, Treatment=1
data['Site_code'] = data['Site'].cat.codes

# View the first few rows
print(data.head())

# Summary statistics
print(data.describe())

# Plot histogram of DLR values
plt.hist(data['DLR'], bins=30, color='blue', alpha=0.7)
plt.title('Distribution of DLR Values')
plt.xlabel('DLR')
plt.ylabel('Frequency')
plt.show()


### 3.3. Estimating Hyperparameters from the Data
### 3.3.1. Estimating the Variance Components
# Calculate overall mean and variance of DLR
mean_DLR = data['DLR'].mean()
var_DLR = data['DLR'].var()

# Estimate within-site variance
within_site_var = data.groupby('Site')['DLR'].var().mean()

# Estimate between-site variance
between_site_var = var_DLR - within_site_var

### Variance Estimates: Provide initial values for setting priors.

### 3.4. Setting Priors Based on Estimates
# Set priors using estimated variances
beta_prior_sd = np.sqrt(var_DLR)           # For fixed effects
sigma_u_prior_sd = np.sqrt(between_site_var)  # For random effects
sigma_prior_sd = np.sqrt(within_site_var)     # For residual variance
# Priors Based on Data: Align priors with observed variability.

## 3.5. Defining the Bayesian LMM
# Extract variables
DLR = data['DLR'].values
Treatment = data['Treatment_code'].values
Site = data['Site_code'].values
n_sites = len(np.unique(Site))

with pm.Model() as bayes_lmm:
    # Priors for fixed effects
    beta_0 = pm.Normal('beta_0', mu=0, sigma=beta_prior_sd)
    beta_Treatment = pm.Normal('beta_Treatment', mu=0, sigma=beta_prior_sd)
    
    # Priors for random intercepts (Site)
    sigma_u = pm.HalfNormal('sigma_u', sigma=sigma_u_prior_sd)
    u = pm.Normal('u', mu=0, sigma=sigma_u, shape=n_sites)
    
    # Residual standard deviation
    sigma = pm.HalfNormal('sigma', sigma=sigma_prior_sd)
    
    # Expected value
    mu = beta_0 + beta_Treatment * Treatment + u[Site]
    
    # Likelihood
    Y_obs = pm.Normal('Y_obs', mu=mu, sigma=sigma, observed=DLR)

### Model Components:
#### Fixed Effects: 𝛽0 and 𝛽Treatment
#### Random Effects: Site-specific intercepts 𝑢
#### Residuals: Modeled with 𝜎.

## 3.6. Fitting the Bayesian LMM
with bayes_lmm:
    trace = pm.sample(4000, tune=2000, chains=4, cores=4, random_seed=123, target_accept=0.95)

# Sampling Parameters: Sufficient iterations and tuning steps for convergence.

## 3.7. Summarizing and Interpreting the Results
# Summarize the posterior distributions
summary = az.summary(trace, var_names=['beta_0', 'beta_Treatment', 'sigma_u', 'sigma'], hdi_prob=0.95)
print(summary)
## beta_Treatment: Estimate of the treatment effect.
## sigma_u: Estimate of the standard deviation of site random effects.
## sigma: Estimate of residual standard deviation.

## 3.8. Detecting Site-Specific Effects
## To detect significant treatment effects at each site, we'll include random slopes.

### 3.8.1. Modifying the Model to Include Random Slopes
with pm.Model() as bayes_lmm_random_slopes:
    # Priors for fixed intercept
    beta_0 = pm.Normal('beta_0', mu=0, sigma=beta_prior_sd)
    
    # Priors for random intercepts and slopes
    sigma_u0 = pm.HalfNormal('sigma_u0', sigma=sigma_u_prior_sd)
    sigma_u1 = pm.HalfNormal('sigma_u1', sigma=sigma_u_prior_sd)
    
    # Random intercepts and slopes
    u0 = pm.Normal('u0', mu=0, sigma=sigma_u0, shape=n_sites)
    u1 = pm.Normal('u1', mu=0, sigma=sigma_u1, shape=n_sites)
    
    # Residual standard deviation
    sigma = pm.HalfNormal('sigma', sigma=sigma_prior_sd)
    
    # Expected value
    mu = beta_0 + u0[Site] + (beta_Treatment + u1[Site]) * Treatment
    
    # Likelihood
    Y_obs = pm.Normal('Y_obs', mu=mu, sigma=sigma, observed=DLR)

### Random Slopes (u1): Allow treatment effect to vary by site.

### 3.8.2. Fitting the Model
with bayes_lmm_random_slopes:
    trace_random_slopes = pm.sample(6000, tune=3000, chains=4, cores=4, random_seed=123, target_accept=0.95)

# Increased Iterations: Necessary due to model complexity.

## 3.8.3. Extracting Site-Specific Treatment Effects
# Extract site-specific treatment effects
u1_samples = trace_random_slopes.posterior['u1']

# Compute mean and credible intervals for each site's treatment effect
site_treatment_effects = pd.DataFrame({
    'Site': data['Site'].cat.categories,
    'Treatment_Effect_Mean': u1_samples.mean(dim=['chain', 'draw']).values,
    'Treatment_Effect_SD': u1_samples.std(dim=['chain', 'draw']).values
})

# Compute 95% credible intervals
hdi = az.hdi(u1_samples, hdi_prob=0.95)
site_treatment_effects['Treatment_Effect_lower'] = hdi.sel(hdi='lower').values
site_treatment_effects['Treatment_Effect_upper'] = hdi.sel(hdi='upper').values

# View the first few site-specific effects
print(site_treatment_effects.head())

## u1_samples: Posterior samples of site-specific treatment effects.
## Credible Intervals: Computed using az.hdi().

### 3.8.4. Identifying Significant Sites
# Determine which sites have significant treatment effects
site_treatment_effects['Significant'] = ((site_treatment_effects['Treatment_Effect_lower'] > 0) |
                                         (site_treatment_effects['Treatment_Effect_upper'] < 0))

# View significant sites
significant_sites = site_treatment_effects[site_treatment_effects['Significant']]
print(significant_sites)

# 3.9. Visualizing Site-Specific Effects
# Plot site-specific treatment effects with credible intervals
plt.figure(figsize=(10, 8))
plt.errorbar(
    x=site_treatment_effects['Treatment_Effect_Mean'],
    y=site_treatment_effects['Site'],
    xerr=[
        site_treatment_effects['Treatment_Effect_Mean'] - site_treatment_effects['Treatment_Effect_lower'],
        site_treatment_effects['Treatment_Effect_upper'] - site_treatment_effects['Treatment_Effect_Mean']
    ],
    fmt='o', ecolor='gray', capsize=3
)
plt.axvline(x=0, color='red', linestyle='--')
plt.xlabel('Treatment Effect')
plt.ylabel('Site')
plt.title('Site-Specific Treatment Effects')
plt.tight_layout()
plt.show()

#  3.10. Model Diagnostics
# Trace plots
az.plot_trace(trace_random_slopes, var_names=['beta_0', 'beta_Treatment', 'sigma_u0', 'sigma_u1', 'sigma'])

# R-hat values
rhat = az.rhat(trace_random_slopes)
print(rhat)

# Posterior predictive checks
with bayes_lmm_random_slopes:
    ppc = pm.sample_posterior_predictive(trace_random_slopes, samples=1000)
az.plot_ppc(az.from_pymc3(posterior_predictive=ppc, model=bayes_lmm_random_slopes))

# Trace Plots: Assess convergence and mixing.
# R-hat Values: Should be close to 1.
# Posterior Predictive Checks: Evaluate model fit.





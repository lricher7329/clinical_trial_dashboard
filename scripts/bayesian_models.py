#!/usr/bin/env Rscript

# Bayesian models for clinical trial analysis

# Install packages if needed
# install.packages(c("rstan", "rstanarm", "brms", "bayesplot", "loo", "tidyverse"))

# Load libraries
library(rstan)
library(rstanarm)
library(brms)
library(bayesplot)
library(loo)
library(tidyverse)

# Set seed for reproducibility
set.seed(42)

# Read data
trial_data <- read.csv("../data/trial_data.csv")
biomarker_data <- read.csv("../data/biomarker_data.csv")

# Convert categorical variables to factors
trial_data$treatment <- as.factor(trial_data$treatment)
trial_data$sex <- as.factor(trial_data$sex)
biomarker_data$treatment <- as.factor(biomarker_data$treatment)
biomarker_data$timepoint <- factor(biomarker_data$timepoint, 
                                  levels = c("Baseline", "Week 6", "Week 12"))

# Create age groups for subgroup analysis
trial_data$age_group <- cut(trial_data$age, 
                           breaks = c(0, 60, 70, 100), 
                           labels = c("<60", "60-70", ">70"))

# Create output directory
dir.create("../results", showWarnings = FALSE)

#' Fit Bayesian model for primary outcome
#' @return A list with the model and summary statistics
fit_primary_outcome_model <- function() {
  cat("Fitting Bayesian model for primary outcome...\n")
  
  # Fit model using rstanarm
  primary_model <- stan_glm(
    primary_outcome ~ treatment + age + sex,
    data = trial_data,
    family = gaussian(),
    prior = normal(0, 2.5),
    prior_intercept = normal(5, 5),
    seed = 123,
    refresh = 0  # Suppress sampling progress output
  )
  
  # Extract posterior samples
  posterior_samples <- as.data.frame(primary_model)
  
  # Calculate summary statistics
  treatment_effect <- posterior_samples$treatmentTreatment
  mean_effect <- mean(treatment_effect)
  ci_effect <- quantile(treatment_effect, c(0.025, 0.975))
  prob_positive <- mean(treatment_effect > 0)
  prob_clinically_significant <- mean(treatment_effect > 1)
  
  # Create summary table
  summary_stats <- data.frame(
    Parameter = "Treatment Effect",
    Mean = mean_effect,
    Lower_CI = ci_effect[1],
    Upper_CI = ci_effect[2],
    Prob_Positive = prob_positive,
    Prob_Significant = prob_clinically_significant
  )
  
  # Save model
  saveRDS(primary_model, "../results/primary_outcome_model.rds")
  
  # Save summary stats
  write.csv(summary_stats, "../results/primary_outcome_summary.csv", row.names = FALSE)
  
  # Create posterior distribution plot
  pdf("../results/primary_outcome_posterior.pdf", width = 8, height = 6)
  mcmc_areas(posterior_samples, pars = "treatmentTreatment", prob = 0.95) +
    labs(title = "Posterior Distribution of Treatment Effect",
         subtitle = "Primary Outcome | 95% Credible Interval")
  dev.off()
  
  # Return results
  return(list(
    model = primary_model,
    summary = summary_stats
  ))
}

#' Fit Bayesian model for biomarker data
#' @return A list with the model and summary statistics
fit_biomarker_model <- function() {
  cat("Fitting Bayesian model for biomarker data...\n")
  
  # Prepare data for brms
  # Convert timepoint to numeric for ease of modeling
  biomarker_data$time_numeric <- as.numeric(biomarker_data$timepoint) - 1  # 0, 1, 2
  
  # Fit longitudinal mixed effects model
  biomarker_model <- brm(
    biomarker_value ~ treatment * time_numeric + (time_numeric | patient_id),
    data = biomarker_data,
    family = gaussian(),
    chains = 2,  # Use fewer chains for speed
    iter = 1000,  # Fewer iterations for speed
    seed = 456,
    refresh = 0
  )
  
  # Extract posterior samples
  posterior_samples <- as.data.frame(biomarker_model)
  
  # Calculate treatment effects at different timepoints
  treatment_effect_baseline <- posterior_samples$b_treatmentTreatment
  treatment_effect_week12 <- posterior_samples$b_treatmentTreatment + 
                           2 * posterior_samples$`b_treatmentTreatment:time_numeric`
  
  # Create summary table
  summary_stats <- data.frame(
    Timepoint = c("Baseline", "Week 12"),
    Mean_Effect = c(mean(treatment_effect_baseline), mean(treatment_effect_week12)),
    Lower_CI = c(quantile(treatment_effect_baseline, 0.025), quantile(treatment_effect_week12, 0.025)),
    Upper_CI = c(quantile(treatment_effect_baseline, 0.975), quantile(treatment_effect_week12, 0.975)),
    Prob_Positive = c(mean(treatment_effect_baseline > 0), mean(treatment_effect_week12 > 0))
  )
  
  # Save model (may be large)
  saveRDS(biomarker_model, "../results/biomarker_model.rds")
  
  # Save summary stats
  write.csv(summary_stats, "../results/biomarker_summary.csv", row.names = FALSE)
  
  # Create conditional effects plot
  pdf("../results/biomarker_conditional_effects.pdf", width = 10, height = 6)
  plot(conditional_effects(biomarker_model, "treatment:time_numeric"), 
       points = TRUE, point_args = list(alpha = 0.3))
  dev.off()
  
  # Return results
  return(list(
    model = biomarker_model,
    summary = summary_stats
  ))
}

#' Run subgroup analysis
#' @return A data frame with subgroup results
run_subgroup_analysis <- function() {
  cat("Running subgroup analysis...\n")
  
  subgroups <- list(
    all = list(filter = rep(TRUE, nrow(trial_data))),
    male = list(filter = trial_data$sex == "Male"),
    female = list(filter = trial_data$sex == "Female"),
    young = list(filter = trial_data$age_group == "<60"),
    middle = list(filter = trial_data$age_group == "60-70"),
    old = list(filter = trial_data$age_group == ">70")
  )
  
  results <- data.frame(
    Subgroup = character(),
    N = integer(),
    Effect_Size = numeric(),
    Lower_CI = numeric(),
    Upper_CI = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(subgroups)) {
    # Subset data
    subset_data <- trial_data[subgroups[[name]]$filter, ]
    
    # Skip if too few observations
    if (nrow(subset_data) < 10) {
      cat("Skipping subgroup", name, "due to insufficient data\n")
      next
    }
    
    # Fit model
    model <- stan_glm(
      primary_outcome ~ treatment,
      data = subset_data,
      family = gaussian(),
      prior = normal(0, 2.5),
      prior_intercept = normal(5, 5),
      seed = 123,
      refresh = 0
    )
    
    # Extract results
    posterior <- as.data.frame(model)
    effect_size <- mean(posterior$treatmentTreatment)
    ci <- quantile(posterior$treatmentTreatment, c(0.025, 0.975))
    p_value <- mean(posterior$treatmentTreatment <= 0)
    
    # Add to results
    results <- rbind(results, data.frame(
      Subgroup = name,
      N = nrow(subset_data),
      Effect_Size = effect_size,
      Lower_CI = ci[1],
      Upper_CI = ci[2],
      P_Value = p_value
    ))
  }
  
  # Save results
  write.csv(results, "../results/subgroup_analysis.csv", row.names = FALSE)
  
  # Create forest plot
  pdf("../results/subgroup_forest_plot.pdf", width = 8, height = 6)
  ggplot(results, aes(x = Effect_Size, y = Subgroup)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Treatment Effect by Subgroup",
         x = "Effect Size",
         y = "Subgroup") +
    theme_minimal()
  dev.off()
  
  return(results)
}

# Run all analyses
if (!interactive()) {
  cat("Running Bayesian analyses for clinical trial data...\n")
  
  primary_results <- fit_primary_outcome_model()
  biomarker_results <- fit_biomarker_model()
  subgroup_results <- run_subgroup_analysis()
  
  cat("Analysis complete. Results saved to ../results/\n")
}

################################################################################################################
############################# R script to simulate four stand-alone 2-arm non-inferiority trial#################
################################################################################################################

# Load required libraries
library(here)
library (tidyverse)
library(sandwich)

# Source helper and model files
# Scenario A
source(here("config_trct_a.R"))

# Scenario B
source(here("config_trct_b.R"))

# Source all functions
source(here("R", "traditionalrct", "recruitment_trct.R"))
source(here("R", "traditionalrct", "outcomes_trct.R"))
source(here("R", "traditionalrct", "type1error_trct.R"))
source(here("R", "traditionalrct", "wrapper_trct.R"))
source(here("R", "traditionalrct", "simulate_trct.R"))

  # Run simulation
res <- simulate_trials_and_results_multi(
  nSims = nSims,
  nPatients = nPatients,
  death0 = death0,
  death1 = death1,
  ni_margin = ni_margin,
  mean_sites_per_month = mean_sites_per_month,
  max_sites = max_sites,
  mean_patients_per_site_per_year = mean_patients_per_site_per_year,
  death_alt = death_alt,
  starts_months = starts_months
)  

################################################################################################################

# Individual trial summaries
batch_summary_df <- res$simulation_results %>%
  group_by(Trial) %>%
  summarise(
    mean_risk_control = mean(risk_control, na.rm = TRUE),
    mean_risk_treatment = mean(risk_treatment, na.rm = TRUE),
    mean_rd = mean(rd, na.rm = TRUE),
    mean_lcl = mean(lcl, na.rm = TRUE),
    mean_ucl = mean(ucl, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE),
    
    relative_bias = (mean_rd - unique(true_rd)) / unique(true_rd),
    relative_bias_percent = ((mean_rd - unique(true_rd)) / unique(true_rd)) * 100,
    
    empirical_se = sd(rd, na.rm = TRUE),
    mean_model_se = mean(se_vec, na.rm = TRUE),
    relative_error = (mean_model_se - empirical_se) / empirical_se,
    relative_error_pct = relative_error * 100,
    
    power = mean(success, na.rm = TRUE),
    type1_error_rate = mean(type1_error, na.rm = TRUE),
    avg_months = mean(total_months, na.rm = TRUE),
    median_duration = median(total_months, na.rm = TRUE),
    average_total_patients = mean(total_patients, na.rm = TRUE),
    median_total_patients = median(total_patients, na.rm = TRUE),
    average_control_patients = mean(control_patients, na.rm = TRUE),
    median_control_patients = median(control_patients, na.rm = TRUE),
    avg_sites = mean(total_sites, na.rm = TRUE),
    mse_rd = mean((rd - (mean_risk_treatment - mean_risk_control))^2, na.rm = TRUE),
    coverage_prob = mean((lcl <= (mean_risk_treatment - mean_risk_control)) &
                           (ucl >= (mean_risk_treatment - mean_risk_control)), na.rm = TRUE) * 100,
    empirical_se = sd(rd, na.rm = TRUE),
    mean_model_se = mean(se_vec, na.rm = TRUE),
    
    # Start month quantiles
    median_start_month = median(Start_Month, na.rm = TRUE),
    q1_start_month     = as.numeric(quantile(Start_Month, 0.25, na.rm = TRUE, type = 7)),
    q3_start_month     = as.numeric(quantile(Start_Month, 0.75, na.rm = TRUE, type = 7)),
    
    # End month quantiles
    median_end_month   = median(End_Month, na.rm = TRUE),
    q1_end_month       = as.numeric(quantile(End_Month, 0.25, na.rm = TRUE, type = 7)),
    q3_end_month       = as.numeric(quantile(End_Month, 0.75, na.rm = TRUE, type = 7)),
    
    .groups = "drop"
  )

cat("\nSummary statistics for each trial:\n")
print(batch_summary_df, n = Inf, width = Inf)   # all rows & all columns

# Overall combined trial summaries
res$summary_stats

# Trial duration
median_duration_per_trial <- res$simulation_results %>%
  mutate(duration_months = End_Month - Start_Month + 1) %>%   # inclusive duration
  group_by(Trial) %>%
  summarise(
    median_duration = median(duration_months, na.rm = TRUE),
    q25_duration = quantile(duration_months, 0.25, na.rm = TRUE),
    q75_duration = quantile(duration_months, 0.75, na.rm = TRUE),
    mean_duration = mean(duration_months, na.rm = TRUE),
    .groups = "drop"
  )

median_duration_per_trial

# Program duration
prog_quartiles <- quantile(res$program_summary$Program_Duration_Months, 
                           probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
prog_summary_stats <- data.frame(
  q25 = prog_quartiles[1],
  median = prog_quartiles[2],
  q75 = prog_quartiles[3]
)

prog_summary_stats

# Staggered trial power (at least one effective trial declared efficacious)
staggered_power_any <- res$simulation_results %>%
  group_by(Sim) %>%                               # one row per simulation replicate
  summarise(any_success = any(success == 1), .groups = "drop") %>%
  summarise(power_any = mean(any_success)) %>%
  pull(power_any)

staggered_power_any

check_success <- res$simulation_results %>%
  group_by(Sim) %>%
  summarise(
    n_success = sum(success == 1),
    any_success = any(success == 1),
    all_success = all(success == 1),
    .groups = "drop"
  )

head(check_success, 10)

# Count number of successful trials per simulation
success_counts <- res$simulation_results %>%
  group_by(Sim) %>%
  summarise(n_success = sum(success == 1), .groups = "drop")

# Frequency table
table(success_counts$n_success)

# Proportions
prop.table(table(success_counts$n_success))

# Staggered trial type I error rate (at least one ineffective trial incorrectly declared efficacious)
# 1. Count number of Type I errors per simulation
type1_counts <- res$simulation_results %>%
  group_by(Sim) %>%
  summarise(n_type1 = sum(type1_error == 1, na.rm = TRUE),
            any_type1 = any(type1_error == 1), .groups = "drop")

# 2. Family-wise Type I error rate (at least one false positive)
staggered_type1_any <- mean(type1_counts$any_type1)
staggered_type1_any

# Frequency table
table(type1_counts$n_type1)

# Proportion table
prop.table(table(type1_counts$n_type1))

# Total control patients is the sum of the average control patients across all simulations
sum_avg_control_patients <- res$simulation_results %>%
  group_by(Sim) %>%
  summarise(total_control_patients = sum(control_patients, na.rm = TRUE)) %>%
  summarise(mean_total_control_patients = mean(total_control_patients)) %>%
  pull(mean_total_control_patients)

sum_avg_control_patients

# Total patients is the sum of the average patients across all simulations
sum_avg_total_patients <- res$simulation_results %>%
  group_by(Sim) %>%
  summarise(total_patients = sum(total_patients, na.rm = TRUE)) %>%
  summarise(mean_total_patients = mean(total_patients)) %>%
  pull(mean_total_patients)

sum_avg_total_patients

#####################################################
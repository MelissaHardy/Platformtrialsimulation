################################################################################################################
############################# R Script to run the Staggered Platform Trial simulation #########################
################################################################################################################
# This code will simulate 5-arm non-inferiority randomised controlled staggered platform trial using 
# a predefined allocation of treatment 0 ("best available therapy") versus treatments 1-4. There are three options for random allocation
# Equal allocation mode (1:1:1:1:1), fixed allocation mode (2:1:1:1:1) and an adaptive sqrt allocation mode (k):1:1:1:1
# where k is the number of arms concurrently open. The primary outcome is "all-cause mortality at 28 days"
# Patients receiving treatment 0 and treatment 1-4 will have a probability of mortality as outlined in the configuration
# Analysis will be performed with a estimated absolute risk difference in the proportion of participants
# with death between the treatment arm and control with no covariate adjustment
# We will simulate 1000 / 10,000 RCT's and store the risk difference, 95% confidence interval, and success vs. the noninferiority 
# margin for each trial and type I error rate
# We will estimate the total duration of the trial based on the number of sites activated per month and
# the number of patients recruited per year

# Load required libraries
library(here)
library(tidyverse)
library(sandwich)

# Source helper and model files

# Scenario A
source(here("config_pt_a.R"))

# Scenario A with the fixed allocation 2:1:1:1:1
source(here("config_pt_a_21.R"))

# Scenario B
source(here("config_pt_b.R"))

# Scenario B with the fixed allocation 2:1:1:1:1
source(here("config_pt_b_21.R"))

# Source staggered platform trial modules
source(here("R", "staggeredpt_cc", "recruitment_staggeredpt_cc.R"))
source(here("R", "staggeredpt_cc", "outcomes_staggeredpt_cc.R"))
source(here("R", "staggeredpt_cc", "type1error_staggeredpt_cc.R"))
source(here("R", "staggeredpt_cc", "wrapper_staggeredpt_cc.R"))
source(here("R", "staggeredpt_cc", "simulate_staggeredpt_cc.R"))

# ---- 1. Equal allocation mode (1:1:1:1:1) ----

# set seed
set.seed(1234)

results_sptf1 <- simulate_trials_spt(
  nSims                        = nSims,
  n_arms                       = n_arms,
  patients_per_experimental_arm = patients_per_experimental_arm,
  death_probs                   = death_probs,
  death_probs_inferior          = death_probs_inferior,
  death_alt                     = death_alt,
  ni_margin                     = ni_margin,
  stagger_months                = stagger_months,
  mean_sites_per_month          = mean_sites_per_month,
  max_sites                     = max_sites,
  mean_patients_per_site_per_year = mean_patients_per_site_per_year,
  true_rd = true_rd,
  allocation_mode = "fixed",
  allocation_ratio = c(1, 1, 1, 1, 1),
  expected_control_arm = expected_control_arm
)

# Print the simulation summary results for each arm
print(results_sptf1$summary_stats)

# Print the simulation summary results the whole platform
print(results_sptf1$trial_level_means)

# Print the risk by arm summary
print(results_sptf1$risk_by_arm_summary)

# Print the number of simulations with a control arm less than the expected
print(results_sptf1$low_control_summary)

################################################################################

# Save results
save(results_sptf1, file = "results/sim_results_staggeredpt_cc_fixed1_a.RData")

# Save random seed
attr(results_sptf1, "seed") <- .Random.seed
saveRDS(.Random.seed, file = "results/sim_seed_staggeredpt_cc_fixed1_a.rds")

# load("sim_seed.RData")       # or
# .Random.seed <- readRDS("sim_seed.rds")

################################################################################

# ---- 2. Fixed allocation mode (2:1:1:1:1) ----

set.seed(1234)

results_sptf2 <- simulate_trials_spt(
  nSims = nSims,
  n_arms = n_arms,
  patients_per_experimental_arm = patients_per_experimental_arm,
  death_probs = death_probs,
  death_probs_inferior = death_probs_inferior,
  death_alt = death_alt,
  ni_margin = ni_margin,
  stagger_months = stagger_months,
  mean_sites_per_month = mean_sites_per_month,
  max_sites = max_sites,
  mean_patients_per_site_per_year = mean_patients_per_site_per_year,
  allocation_mode = "fixed",
  allocation_ratio = c(2, 1, 1, 1, 1),
  true_rd = true_rd,
  expected_control_arm = expected_control_arm
)

# Print the simulation summary results for each arm
print(results_sptf2$summary_stats)

# Print the simulation summary results the whole platform
print(results_sptf2$trial_level_means)

# Print the risk by arm summary
print(results_sptf2$risk_by_arm_summary)

# Print the number of simulations with a control arm less than the expected
print(results_sptf2$low_control_summary)

################################################################################

# Save results
save(results_sptf2, file = here("results", "sim_results_results_staggered_fixed2_a.RData"))

# Save random seed
attr(results_sptf2, "seed") <- .Random.seed
saveRDS(.Random.seed, file = "results/sim_seed_staggeredpt_cc_fixed2_a.rds")

# load("sim_seed.RData")       # or
# .Random.seed <- readRDS("sim_seed.rds")

# ---- 3. sqrt(k) allocation mode ----

set.seed(1234)

results_sqrtk <- simulate_trials_spt(
  nSims                        = nSims,
  n_arms                       = n_arms,
  patients_per_experimental_arm = patients_per_experimental_arm,
  death_probs                   = death_probs,
  death_probs_inferior          = death_probs_inferior,
  death_alt                     = death_alt,
  ni_margin                     = ni_margin,
  stagger_months                = stagger_months,
  mean_sites_per_month          = mean_sites_per_month,
  max_sites                     = max_sites,
  mean_patients_per_site_per_year = mean_patients_per_site_per_year,
  allocation_mode = "sqrtk",
  true_rd = true_rd,
  expected_control_arm = expected_control_arm
)

# Print the simulation summary results for each arm
print(results_sqrtk$summary_stats)

# Print the simulation summary results the whole platform
print(results_sqrtk$trial_level_means)

# Print the risk by arm summary
print(results_sqrtk$risk_by_arm_summary)

# Print the number of simulations with a control arm less than the expected
print(results_sqrtk$low_control_summary)

################################################################################

save(results_sqrtk, file = here("results", "sim_results_staggered_sqrtk_b.RData"))

# Save random seed
attr(results_sqrtk, "seed") <- .Random.seed
saveRDS(.Random.seed, file = "results/sim_seed_staggeredpt_cc_sqrtk_b.rds")


################################################################################################################

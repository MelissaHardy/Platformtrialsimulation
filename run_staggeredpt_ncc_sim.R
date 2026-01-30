################################################################################################################
############################# R Script to run the Staggered Platform Trial simulation #########################
################################################################################################################
# This code will simulate 5-arm non-inferiority randomised controlled staggered platform trial using 
# a non-concurrent control with treatment 0 ("best available therapy") versus treatment 1 - 4. There random allocation will be 1:1:1:1:1 until each arm reaches capacity
# The primary outcome is "all-cause mortality at 28 days"
# Patients receiving treatment 0 and treatment 1-4 will have a probability of mortality as outlined in the configuration
# Analysis will be performed with a estimated absolute risk difference in the proportion of participants
# with death between treatment 1-4 and control with no covariate adjustment
# We will run 1000 / 10,000 simulated RCT's and store the risk difference, 95% confidence interval, and success vs. the noninferiority margin for each trial
# We will also estimate the total duration of the trial based on the number of sites activated per month and
# the number of patients recruited per year

# Load required libraries
library(here)
library(tidyverse)
library(sandwich)

# Source helper and model files
source(here("config_pt_a.R"))

# Source helper and model files
source(here("config_pt_b.R"))

# Source staggered platform trial modules
source(here("R", "staggeredpt_ncc", "recruitment_staggeredpt_ncc.R"))
source(here("R", "staggeredpt_ncc", "outcomes_staggeredpt_ncc.R"))
source(here("R", "staggeredpt_ncc", "type1error_staggeredpt_ncc.R"))
source(here("R", "staggeredpt_ncc", "wrapper_staggeredpt_ncc.R"))
source(here("R", "staggeredpt_ncc", "simulate_staggeredpt_ncc.R"))

# set seed
set.seed(12345)


results_sptncc <- simulate_trials_sptncc(
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
  mean_patients_per_site_per_year = mean_patients_per_site_per_year
)

# Print the simulation summary results for each arm
print(results_sptncc$summary_stats)

# Print the simulation summary results the whole platform
print(results_sptncc$trial_level_means)

recruitment_df <- results_sptncc$recruitment_details[[1]]

################################################################################################################

# Save results
save(results_sptncc, file = "results/sim_results_staggeredpt_ncc_b.RData")

# Save random seed
attr(results_sptncc, "seed") <- .Random.seed
saveRDS(.Random.seed, file = "results/sim_seed_staggeredpt_ncc_b.rds")

# load("sim_seed.RData")       # or
# .Random.seed <- readRDS("sim_seed.rds")

################################################################################################################

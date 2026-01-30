################################################################################################################
############################# R Script to run the Parallel Platform Trial simulation #########################
################################################################################################################
# This code will simulate a 5-arm non-inferiority randomised controlled platform trial using 1:1 allocation 
# of participants to treatment 0 ("best available therapy") versus treatments 1 - 4
# The primary outcome is "all-cause mortality at 28 days" after randomisation
# Patients receiving treatment 0 and treatments 1-4 will have a probability of mortality as outlined in the configuration
# Analysis will be performed with a estimated absolute risk difference in the proportion of participants
# with death between the treatments 1-4 and control with no covariate adjustment
# We will run 1000 / 10,000 simulated RCT's and store the risk difference, 95% confidence interval, and success vs. the noninferiority margin for each trial
# We will also estimate the total duration of the trial based on the number of sites activated per month and
# the number of patients recruited per year

# Load required libraries
library(here)
library(tidyverse)
library(sandwich)

# Source helper and model files
# Scenario A
source(here("config_pt_a.R"))

# Scenario B
source(here("config_pt_b.R"))

# Source parallel platform trial modules
source(here("R", "parallel_pt", "recruitment_parallelpt.R"))
source(here("R", "parallel_pt", "outcomes_parallelpt.R"))
source(here("R", "parallel_pt", "type1error_parallelpt.R"))
source(here("R", "parallel_pt", "wrapper_parallelpt.R"))
source(here("R", "parallel_pt", "simulate_parallelpt.R"))

set.seed(1234)

results_ppt <- simulate_trials_ppt(nSims,
                               n_arms,
                               patients_per_experimental_arm,
                               death0, 
                               death_probs, 
                               ni_margin,
                               mean_sites_per_month, 
                               max_sites, 
                               mean_patients_per_site_per_year, 
                               death_alt,
                               true_rd)


# Print the simulation summary results for each arm
print(results_ppt$summary_stats)

# Print the simulation summary results for the platform
print(results_ppt$trial_level_means)

################################################################################################################

# Save results
save(results_ppt, file = "results/sim_results_parallelpt_b.RData")

# Save random seed
attr(results_ppt, "seed") <- .Random.seed
saveRDS(.Random.seed, file = "results/sim_seed_parallelpt_b.rds")

# load("sim_seed.RData")       # or
# .Random.seed <- readRDS("sim_seed.rds")

################################################################################################################


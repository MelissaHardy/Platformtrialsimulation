################################################################################################################
###### R Script to simulate the type 1 error rate from a Staggered Platform Trial with 5 arms###################
################################################################################################################
# This code will simulate the type 1 error rate from a 5-arm non-inferiority randomised controlled platform trial 
# the predefined randomisation allocation to treatment 0 ("best available therapy") versus treatment 1 versus 
# treatment 2 versus treatment 3 and treatment 4 dependent on the arms that are actively recruiting
# Patients receiving treatment 1-4 will have an inferior mortality risk as outlined below
# Non-inferiority margin will be predefined
# Analysis will be performed with a estimated absolute risk difference in the proportion of participants
# with death between the investigational agent and control with no covariate adjustment
# We will run x simulated RCT's and store the risk difference, 95% confidence interval, and success vs. 
# the noninferiority margin for each trial

simulate_type1_error_spt <- function(recruitment_data,
                                     death_probs_inferior,
                                     ni_margin,
                                     stagger_months,
                                     patients_per_experimental_arm) {
  
  # Step 2: Run a single trial simulation under the null
  sim_result <- simulate_trial_outcomes_spt(
    recruitment_data = recruitment_data,
    death_probs = death_probs_inferior,
    ni_margin = ni_margin,
    stagger_months = stagger_months,
    patients_per_experimental_arm = patients_per_experimental_arm
  )
  
  # Step 3: Extract relevant results
  type1_errors <- sim_result$success              # 1 = falsely declared non-inferiority
  ucls <- sim_result$ucl                          # Upper confidence limits
  comparison_names <- sim_result$comparisons      # e.g. "Arm1 vs Control"
  trialdata_inf <- sim_result$trial_data          # Patient-level outcomes
  
  # Step 4: Return results
  return(list(
    ucl_inf = ucls,
    type1_error = type1_errors,
    comparisons = comparison_names,
    trialdata_inf = trialdata_inf
  ))
}

######################################################

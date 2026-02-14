################################################################################################################
######R Script for simulating type 1 error rate for a Staggered Platform Trial with a non-concurrent control###
################################################################################################################
# This code will simulate the type 1 error rate with a 5-arm non-inferiority randomised controlled platform trial 
# using a non-concurrent control and an inferior mortality rate for the investigational treatment arms.

simulate_type1_error_sptncc <- function(recruitment_data,
                                     death_probs_inferior,
                                     ni_margin,
                                     stagger_months,
                                     patients_per_experimental_arm) {
  
  # Step 2: Run a single trial simulation under the null
  sim_result <- simulate_trial_outcomes_sptncc(
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

################################################################################################################

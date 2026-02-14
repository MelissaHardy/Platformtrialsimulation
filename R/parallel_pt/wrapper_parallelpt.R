################################################################################################################
####### R script for a wrapper function that simulates a Parallel Platform Trial with 5 arms####################
################################################################################################################
# This code will wrap the recruirment, outcomes and type 1 error functions to simulate a 5-arm non-inferiority 
# randomised controlled parallel platform trial

simulate_full_trial_ppt <- function(n_arms,
                                patients_per_experimental_arm,
                                death0, 
                                death_probs, 
                                ni_margin,
                                mean_sites_per_month, 
                                max_sites, 
                                mean_patients_per_site_per_year,
                                death_alt, true_rd) {
  
  # Step 1: Compute total sample size
  nPatients <- n_arms * patients_per_experimental_arm
  
  # Step 2: Recruitment simulation
  recruitment <- simulate_recruitment_ppt(nPatients, 
                                          mean_sites_per_month, 
                                          max_sites, 
                                          mean_patients_per_site_per_year)

  # Step 3: Trial outcomes simulation
  outcomes <- simulate_trial_outcomes_ppt(n_arms, 
                                          patients_per_experimental_arm, 
                                          death_probs, 
                                          ni_margin)
  
  # Step 4: Type I error simulation
  type1 <- simulate_type1_error_ppt(outcomes$trialdata, 
                                    death0, 
                                    death_alt, 
                                    ni_margin)
  
  # Aggregate and return results
  list(
    recruitment = recruitment,
    outcomes = outcomes,
    type1_error = type1,
    site_accrual = recruitment$site_accrual
  )
}

################################################################################################################


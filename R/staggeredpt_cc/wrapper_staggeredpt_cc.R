################################################################################################################
######R script for a wrapper function that simulates a Staggered Platform Trial ################################
################################################################################################################
# This code will wrap the recruirment, outcomes and type 1 error functions to simulate a 5-arm non-inferiority 
# randomised controlled platform trial

simulate_full_trial_spt <- function(n_arms,
                                    patients_per_experimental_arm,
                                    death_probs,
                                    death0,
                                    death_alt,
                                    death_probs_inferior,
                                    ni_margin,
                                    stagger_months,
                                    mean_sites_per_month,
                                    max_sites,
                                    mean_patients_per_site_per_year,
                                    allocation_mode = c("fixed", "sqrtk"),
                                    allocation_ratio = NULL) {
  
  # Step 1: Simulate recruitment
  recruitment <- simulate_recruitment_spt(
    n_arms = n_arms,
    patients_per_experimental_arm = patients_per_experimental_arm,
    mean_sites_per_month = mean_sites_per_month,
    max_sites = max_sites,
    mean_patients_per_site_per_year = mean_patients_per_site_per_year,
    stagger_months = stagger_months,
    allocation_mode = allocation_mode,
    allocation_ratio = allocation_ratio
  )
  
  recruitment_data <- recruitment$recruitment_details
  
  # Step 2: Simulate trial outcomes under alternative hypothesis
  outcomes <- simulate_trial_outcomes_spt(
    recruitment_data = recruitment$recruitment_details,
    death_probs = death_probs,
    ni_margin = ni_margin,
    stagger_months = stagger_months,
    patients_per_experimental_arm = patients_per_experimental_arm
    )
  
  # Step 3: Simulate Type I error (under null: all treatments worse)
  type1 <- simulate_type1_error_spt(
    recruitment_data = recruitment$recruitment_details,
    death_probs_inferior = death_probs_inferior,
    ni_margin = ni_margin,
    stagger_months = stagger_months,
    patients_per_experimental_arm = patients_per_experimental_arm
  )
  
  # Return combined results
  list(
    recruitment = recruitment,
    outcomes = outcomes,
    type1_error = type1,
    site_accrual = recruitment$site_accrual
  )
}

################################################################################################################

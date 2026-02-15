################################################################################################################
############################# R script for the master function to run a full simulation trial############################
################################################################################################################

simulate_full_trial <- function(nPatients, death0, death1, ni_margin,
                                death_alt) {
 
  # Step 1: Trial outcomes simulation
  outcomes <- simulate_trial_outcomes(nPatients, death0, death1, ni_margin)
  
  # Step 2: Type I error simulation
  type1 <- simulate_type1_error(outcomes$trialdata, death0, death_alt, ni_margin)
  
  # Aggregate and return results
  list(
    outcomes = outcomes,
    type1_error = type1
    )
}

######################################################


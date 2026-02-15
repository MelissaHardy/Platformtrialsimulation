###############################################################################################
######## R script for simulating outcomes from a stand-alone 2-arm non-inferiority trial########
###############################################################################################

# This code will simulate outcomes from a 2-arm non-inferiority randomised controlled trial 
# using 1:1 allocation of participants to treatment 0 ("best available therapy") versus 
# treatment 1 ("new treatment")
# The primary outcome is "all-cause mortality at 28 days" after randomisation
# Patients receiving treatment 0 and treatment 1 will have a probability of mortality 
# as outlined in the parameters
# Non-inferiority margin will be outlined in the parameters
# Analysis will be performed with a estimated absolute risk difference in the proportion of participants
# with death between the investigational agent and control with no covariate adjustment
# We will run x simulated RCT's and store the risk difference, 95% confidence interval, 
# and success vs. the noninferiority margin for each trial

simulate_trial_outcomes <- function(nPatients, death0, death1, ni_margin) {
  # Number of patients per group
  n_per_arm <- nPatients / 2
  
  # Generate trial data: assign patient IDs, treatment, and outcomes
  pid <- 1:nPatients
  
  # Assign fixed number to each group: 0 = control, 1 = treatment
  treatment <- rep(c(0, 1), each = n_per_arm)
  
  # Shuffle treatment assignment randomly
  treatment <- sample(treatment)
  
  # Assign the probability of the outcome for each patient based on treatment allocation
  deathprob <- ifelse(treatment == 0, death0, death1)
  
  # Simulate each patient's outcome as a random draw from a binomial distribution
  death <- rbinom(nPatients, 1, deathprob)
  
  # Create a data frame to store trial data
  trialdata <- data.frame(pid = pid, treatment = treatment, death = death)
  
  # Fit a generalised linear regression model to the trial's simulated data
  model <- glm(death ~ treatment, family = gaussian(link = "identity"), data = trialdata)
  
  # Extract estimate and robust SE
  est <- coef(model)["treatment"]
  se <- sqrt(vcovHC(model, type = "HC0")["treatment", "treatment"])  # robust SE
  
  # Extract trial results
  rd <- round(est, 4)                        # Risk difference
  se <- round(se, 4)                         # Standard error
  lcl <- round(rd - 1.96 * se, 4)            # Lower limit of 95% CI
  ucl <- round(rd + 1.96 * se, 4)            # Upper limit of 95% CI
  success <- ifelse(ucl < ni_margin, 1, 0)   # Non-inferiority success (1 = success, 0 = failure)
  
  # Calculate mean predicted probabilities for each group
  predicted_probabilities <- predict(model, type = "response")
  risk_control <- round(mean(predicted_probabilities[treatment == 0]), 4)
  risk_treatment <- round(mean(predicted_probabilities[treatment == 1]), 4)
  
  list(
    trialdata = trialdata,
    rd = rd,
    lcl = lcl,
    ucl = ucl,
    se = se,
    success = success,
    risk_control = risk_control,
    risk_treatment = risk_treatment,
    se = se
  )
}

##############################################################
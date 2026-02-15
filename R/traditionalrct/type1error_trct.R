####################################################################################################
### R script for simulating outcomes and type 1 error rate from a 2-arm non-inferiority trial#######
####################################################################################################

# This code will simulate recruitment for a 2-arm non-inferiority randomised controlled trial 
# using 1:1 allocation of participants to treatment 0 ("best available therapy") versus 
# treatment 1 ("new treatment")
# Patients receiving treatment 0 and treatment 1
# This function will estimate the type 1 error rate with an inferior mortality rate

######################################################
# Function to simulate type I error scenario
######################################################

simulate_type1_error <- function(trialdata, death0, death_alt, ni_margin) {
 
   # Use the same treatment assignments as in trialdata
  treatment <- trialdata$treatment
  pid <- trialdata$pid
  
  # Adjust outcome probabilities to simulate an inferior treatment
  deathprob_inf <- ifelse(treatment == 0, death0, death_alt)
  death_inf <- rbinom(length(pid), 1, deathprob_inf)
  trialdata_inf <- data.frame(pid = pid, treatment = treatment, death = death_inf)
  
  # Fit the model for type I error analysis
  model_inf <- glm(death ~ treatment, family = gaussian(link = "identity"), data = trialdata_inf)
  
  # Extract estimate and robust SE
  est <- coef(model_inf)["treatment"]
  se <- sqrt(vcovHC(model_inf, type = "HC0")["treatment", "treatment"])  # robust SE
  
  rd <- round(est, 4)                        # Risk difference
  ucl_inf <- round(rd + 1.96 * se, 4)
  
  type1_error <- ifelse(ucl_inf < ni_margin, 1, 0)
  
  list(
    ucl_inf = ucl_inf,
    type1_error = type1_error,
    trialdata_inf = trialdata_inf
  )
}

##############################################################
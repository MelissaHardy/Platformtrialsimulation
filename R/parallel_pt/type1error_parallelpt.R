################################################################################################################
#################### R script to simulate the type 1 error rate for a Parallel Platform Trial with 5 arms######
################################################################################################################
# This code will simulate the type 1 error rate for a 5-arm non-inferiority randomised controlled platform trial 
# using 1:1 allocation of participants to treatment 0 ("best available therapy") versus treatment 1 versus treatment 
# 2 versus treatment 3 and treatment 4
# It will calculate the probability of incorrectly concluding non-inferiority when the experimental treatments 
# are actually inferior to control
# Patients receiving treatment 1-4 will have a inferior probability of mortality risk

simulate_type1_error_ppt <- function(trialdata, death0, death_alt, ni_margin) {
  
  # Use the same treatment assignments as in trialdata
  treatment <- trialdata$treatment
  pid <- trialdata$pid
  
  # Assign true probabilities: control (0) = death0, arms 1–4 = inferior death_alt
  deathprob_inf <- ifelse(treatment == 0, death0, death_alt)
  death_inf <- rbinom(length(pid), 1, deathprob_inf)
  trialdata_inf <- data.frame(pid = pid,
                              treatment = relevel(factor(treatment), ref = "0"),  # ensure control is ref
                              death = death_inf)
  
  # Fit the model for type I error analysis
  model_inf <- glm(death ~ treatment, family = gaussian(link = "identity"), data = trialdata_inf)
  beta       <- stats::coef(model_inf)
  vcov_rob   <- sandwich::vcovHC(model_inf, type = "HC0")
  coef_names <- names(beta)
  
  # Indices for treatment-vs-control coefficients (exclude intercept)
  treat_idx  <- which(coef_names != "(Intercept)")
  
  # Risk differences vs Treatment0 and robust SEs
  rd     <- beta[treat_idx]
  se <- sqrt(diag(vcov_rob))[treat_idx]
  
  ucl_inf  <- rd + 1.96 * se
  type1_error <- as.integer(ucl_inf < ni_margin)   # 1 = false NI claim
  
  # Names straight from the model (e.g., "treatment1","treatment2",...)
  comparisons <- coef_names[treat_idx]
  
  list(
    ucl_inf = ucl_inf,
    type1_error = type1_error,
    comparisons = comparisons,
    trialdata_inf = trialdata_inf
  )
}

##########################################################################################################


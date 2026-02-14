################################################################################################################
############################# R script to simulate outcomes of a Parallel Platform Trial with 5 arms###########
################################################################################################################
# This code will simulate a 5-arm non-inferiority randomised controlled platform trial using 1:1 allocation of participants
# to treatment 0 ("best available therapy") versus treatment 1 versus treatment 2 versus treatment 3 and treatment 4
# The primary outcome is "all-cause mortality at 28 days" after randomisation
# Patients receiving treatment 0 and treatment 1-4 will have a probability of mortality as outlined in the design parameters
# Analysis will be performed with a estimated absolute risk difference in the proportion of participants
# with death between the investigational agent and control with no covariate adjustment
# We will run x simulated RCT's and store the risk difference, 95% confidence interval, and success vs. the noninferiority margin for each trial

simulate_trial_outcomes_ppt <- function(n_arms, 
                                        patients_per_experimental_arm, 
                                        death_probs, 
                                        ni_margin) {
  
  # Step 1: Set up and randomisation
  nPatients <- n_arms * patients_per_experimental_arm
  
  # Equal randomisation
  treatment <- sample(rep(0:(n_arms - 1), each = patients_per_experimental_arm))
  pid         <- seq_len(nPatients)
  
  # Step 2: Simulate death outcomes based on arm-specific probabilities
  death_prob <- death_probs[treatment + 1]  # arm-specific probabilities
  death      <- rbinom(nPatients, size = 1, prob = death_prob)
  
  # Make treatment a FACTOR with "0" as the first level (control)
  treatment <- factor(treatment, levels = 0:(n_arms - 1))
  
  trialdata <- data.frame(pid = pid,
                          treatment = treatment,
                          death = death)
  
  # Step 3: Fit a linear probability model (treatment as factor)
  model <- glm(death ~ treatment, family = gaussian(link = "identity"), data = trialdata)
  
  # Step 4: Extract risk difference estimates for Arms 1–4 vs Arm 0 and robust SE
  beta       <- stats::coef(model)
  vcov_rob   <- sandwich::vcovHC(model, type = "HC0")
  coef_names <- names(beta)
  treat_idx  <- which(coef_names != "(Intercept)")
  
  # treatment-vs-control coefficients (exclude intercept)
  comparisons <- coef_names[treat_idx]   # e.g. "treatment1", "treatment2", ...
  
  # Risk differences vs Treatment0 and robust SEs
  rd     <- beta[treat_idx]
  se_rob <- sqrt(diag(vcov_rob))[treat_idx]
  
  # 95% CIs and NI success per arm
  lcl     <- rd - 1.96 * se_rob
  ucl     <- rd + 1.96 * se_rob
  success <- as.integer(ucl < ni_margin)   # declare NI if upper CI < margin
  
  # Step 5: Predicted risk by arm
  predicted <- predict(model, type = "response")
  risk_by_arm <- round(tapply(predicted, trialdata$treatment, mean), 4)
  
  # Patients per arm
  patients_by_arm <- as.integer(table(trialdata$treatment))
  names(patients_by_arm) <- levels(trialdata$treatment)
  
  list(
    trialdata = trialdata,
    comparisons = comparisons,
    rd = rd,
    se = se_rob,         # robust (HC0) SEs
    lcl = lcl,
    ucl = ucl,
    success = success,
    risk_by_arm    = risk_by_arm,
    patients_by_arm = rep(patients_per_experimental_arm, n_arms)
  )
}

################################################################################################################

################################################################################################################
############################# R script for simulating outcomes from a Staggered Platform Trial with 5 arms#####
################################################################################################################
# This code will simulate a 5-arm non-inferiority randomised controlled staggered platform trial the randomisation
# allocation specified in the recruitment function. It will randomise treatment 0 ("best available therapy") 
# versus treatment 1 versus treatment 2 versus treatment 3 and treatment 4 dependent on the arms that are actively 
# recruiting
# The primary outcome is "all-cause mortality at 28 days" after randomisation
# Patients receiving treatment 0 and treatment 1-4 will have a probability of mortality as outlined below
# Non-inferiority margin will be outlined below
# Analysis will be performed with a estimated absolute risk difference in the proportion of participants
# with death between the investigational agent and control with no covariate adjustment
# We will run x simulated RCT's and store the risk difference, 95% confidence interval, and success vs. the noninferiority margin for each trial

simulate_trial_outcomes_spt <- function(recruitment_data, 
                                        death_probs, 
                                        ni_margin, 
                                        stagger_months, 
                                        patients_per_experimental_arm) {
  
  n_arms <- length(death_probs)  # e.g., 5 for Arm0 to Arm4
  
  # Step 1: Expand recruitment data to patient-level format
  trial_data <- recruitment_data %>%
    select(Month, matches("^Arm\\d+$")) %>%
    pivot_longer(cols = -Month, names_to = "arm", values_to = "monthly_count") %>%
    filter(monthly_count > 0) %>%
    uncount(weights = monthly_count) %>%
    mutate(
      arm = as.numeric(gsub("Arm", "", arm)),  # Arm0 = 0, etc.
      pid = row_number()
    )
  
  # Step 2: Simulate death outcomes
  trial_data$death <- rbinom(nrow(trial_data), 1, prob = death_probs[trial_data$arm + 1])
  
  # Step 3: Initialize result containers
  rd <- se <- lcl <- ucl <- success <- rep(NA, n_arms - 1)
  comparisons <- paste0("Arm", 1:(n_arms - 1), " vs Control")
  risk_by_arm_list <- vector("list", n_arms - 1)
  
  # Step 4: Loop through experimental arms
  for (j in 1:(n_arms - 1)) {
    arm_index <- j
    arm_col <- paste0("Arm", arm_index)
    
    start_month <- stagger_months * max(0, arm_index - 1)
    
    arm_data <- trial_data %>%
      filter(arm == arm_index & Month >= start_month) %>%
      arrange(Month) %>%
      mutate(cumulative = cumsum(rep(1, n())))
    
    cutoff_month <- arm_data %>%
      filter(cumulative >= patients_per_experimental_arm) %>%
      slice(1) %>%
      pull(Month)
    
    eligible_months <- start_month:cutoff_month
    
    # Subset patients from concurrent control and treatment arm
    sub <- trial_data %>%
      filter(Month %in% eligible_months & arm %in% c(0, arm_index)) %>%
      mutate(treatment = factor(ifelse(arm == 0, "Control", paste0("Arm", arm_index)),
                                levels = c("Control", paste0("Arm", arm_index))))
    
    # Skip if both arms not present
    if (length(unique(sub$treatment)) < 2) next
    
    # Fit linear probability model
    model <- glm(death ~ treatment, data = sub, family = gaussian(link = "identity"))
    
    # Extract risk difference estimates for Arms 1–4 vs Arm 0
    beta <- coef(model)
    treat_coef_name <- setdiff(names(beta), "(Intercept)")[1]
    
    # Robust SE (HC0)
    vcov_rob   <- sandwich::vcovHC(model, type = "HC0")
    
    # Estimate, robust SE, CI, NI decision
    rd_j <- unname(beta[treat_coef_name])
    se_j <- sqrt(vcov_rob[treat_coef_name, treat_coef_name])
    lcl_j <- rd_j - 1.96 * se_j
    ucl_j <- rd_j + 1.96 * se_j
    success_j <- as.integer(ucl_j < ni_margin)
    
      # Store outputs
      rd[arm_index] <- round(rd_j, 4)
      se[arm_index] <- round(se_j, 4)
      lcl[arm_index] <- round(lcl_j, 4)
      ucl[arm_index] <- round(ucl_j, 4)
      success[arm_index] <- success_j
      
      # Predicted probabilities
      predicted <- predict(model, type = "response")
      risk_estimates <- tapply(predicted, sub$treatment, mean)
      n_patients <- table(sub$treatment)
      
      risk_by_arm_list[[j]] <- data.frame(
        Arm = names(risk_estimates),
        Predicted_Prob = round(as.numeric(risk_estimates), 4),
        N = as.integer(n_patients[names(risk_estimates)])
      )
    }
  
  # Step 5: Return results
  return(list(
    trial_data = trial_data,
    rd = rd,
    se = se,
    lcl = lcl,
    ucl = ucl,
    success = success,
    comparisons = comparisons,
    risk_by_arm = risk_by_arm_list
  ))
}

##############################################################

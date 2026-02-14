################################################################################################################
######################### R script to simulate multiple Parallel Platform Trials with 5 arms ##################
################################################################################################################
# This code will simulate multiple 5-arm non-inferiority randomised controlled platform trials and aggregate
# the summary statistics and save the lists

simulate_trials_ppt <- function(nSims, 
                            n_arms,
                            patients_per_experimental_arm,
                            death0, 
                            death_probs, 
                            ni_margin,
                            mean_sites_per_month, 
                            max_sites, 
                            mean_patients_per_site_per_year, 
                            death_alt,
                            true_rd) {
  
  nPatients <- n_arms * patients_per_experimental_arm
  
  # Initialize vectors and lists to store aggregated results
  trialnum <- numeric(nSims)
  total_sites <- numeric(nSims)
  total_patients <- numeric(nSims)
  control_patients <- numeric(nSims)
  total_months <- numeric(nSims)
  type1_error <- vector("list", nSims)
  recruitment_details <- vector("list", nSims)
  trial_data_list <- vector("list", nSims)  # List to store full trial data
  trial_data_inf_list <- vector("list", nSims)  # List to store type 1 error trial data
  site_accrual_list <- vector("list", nSims)
  
  # Initialize lists for outcome results (vector per arm)
  rd_list <- list()
  lcl_list <- list()
  ucl_list <- list()
  ci_width <- list()
  se_list <- list()
  success_list <- list()
  risk_by_arm_list <- list()
  
  # Loop over simulations
  for (i in seq_len(nSims)) {
    sim <- simulate_full_trial_ppt(n_arms, 
                                   patients_per_experimental_arm, 
                                   death0, death_probs, 
                                   ni_margin,
                                   mean_sites_per_month, 
                                   max_sites, 
                                   mean_patients_per_site_per_year, 
                                   death_alt, 
                                   true_rd)
    
    trialnum[i] <- i
    total_sites[i] <- sim$recruitment$total_sites
    total_patients[i] <- sim$recruitment$total_patients
    control_patients[i] <- sum(sim$outcomes$trialdata$treatment == 0)
    total_months[i] <- sim$recruitment$total_months
    recruitment_details[[i]] <- sim$recruitment$recruitment_details
    trial_data_list[[i]] <- sim$outcomes$trialdata
    trial_data_inf_list[[i]] <- sim$type1_error$trialdata_inf
    type1_error[[i]] <- sim$type1_error$type1_error
    site_accrual_list[[i]] <- sim$site_accrual
    
    # Per-arm results
    rd_list[[i]] <- sim$outcomes$rd
    lcl_list[[i]] <- sim$outcomes$lcl
    ucl_list[[i]] <- sim$outcomes$ucl
    ci_width[[i]] <- sim$outcomes$ucl - sim$outcomes$lcl
    se_list[[i]] <- sim$outcomes$se
    success_list[[i]] <- sim$outcomes$success
    risk_by_arm_list[[i]] <- sim$outcomes$risk_by_arm
  }
  
  # Convert list of vectors to data frames for summary
  rd_df <- do.call(rbind, rd_list)
  lcl_df <- do.call(rbind, lcl_list)
  ucl_df <- do.call(rbind, ucl_list)
  ci_width_df <- do.call(rbind, ci_width)
  se_df <- do.call(rbind, se_list)
  success_df <- do.call(rbind, success_list)
  type1_df <- do.call(rbind, type1_error)
  risk_df    <- do.call(rbind, risk_by_arm_list)  # nSims x n_arms
  
  # Empirical standard error for each treatment arm - standard deviation across simulations
  empirical_se <- apply(rd_df, 2, sd)  # vector of length (n_arms - 1)
  
  # Overall mean empirical SE across all treatment arms
  overall_empirical_se <- mean(empirical_se)
  
  # Expand single true_rd into a vector if needed
  if (length(true_rd) == 1) {
    true_rd <- rep(true_rd, n_arms - 1)
  }
  
  # Compute coverage probability for each treatment arm
  coverage_by_arm <- colMeans((lcl_df <= true_rd) & (ucl_df >= true_rd)) * 100

  # Coverage across platform: at least one arm per trial with CI covering the true RD
  platform_coverage <- mean(apply((lcl_df <= true_rd) & (ucl_df >= true_rd), 1, any)) * 100
  
  # Calculate per-simulation power
  arm_success_rate <- rowMeans(success_df)  # gives a vector of length nSims
  mean_platform_power <- mean(arm_success_rate)
  
  # Platform-level power (at least one effective arm declared efficacious)
  overall_platform_power_any <- mean(apply(success_df, 1, function(x) any(x == 1)))
  
  # Calculate per-simulation type 1 error rate
  arm_error_rate <- rowMeans(type1_df)  # gives a vector of length nSims
  mean_platform_errorrate <- mean(arm_error_rate)
  
  # Platform-level error rate (at least one ineffective arm incorrectly declared efficacious)
  overall_platform_type1_any <- mean(apply(type1_df, 1, function(x) any(x == 1)))
  
  # Labels
  arm_labels <- paste0("Treatment", 0:(n_arms - 1))  # including control
  
  # Initialize tidy data frame
  long_summary <- data.frame(
    sim = rep(1:nSims, each = n_arms),
    arm = rep(arm_labels, times = nSims),
    risk = as.vector(t(risk_df)),
    rd = NA,
    lcl = NA,
    ucl = NA,
    ci_width = NA,
    success = NA,
    type1error = NA
  )
  
  # Fill in RD, CI, and success for treatment arms (Arm1 to Arm4)
  for (j in 1:(n_arms - 1)) {
    row_indices <- seq(j + 1, nrow(long_summary), by = n_arms)  # skip Arm0
    long_summary$rd[row_indices]      <- as.vector(rd_df[, j])
    long_summary$lcl[row_indices]     <- as.vector(lcl_df[, j])
    long_summary$ucl[row_indices]     <- as.vector(ucl_df[, j])
    long_summary$ci_width[row_indices]     <- as.vector(ci_width_df[, j])
    long_summary$success[row_indices] <- as.vector(success_df[, j])
    long_summary$type1error[row_indices] <- as.vector(type1_df[, j])
  }
  
  # Calculate summary statistics per arm
  mean_rd <- colMeans(rd_df)
  mean_lcl <- colMeans(lcl_df)
  mean_ucl <- colMeans(ucl_df)
  mean_ci_width <- colMeans(ci_width_df)
  mean_se <- colMeans(se_df)
  power <- colMeans(success_df)
  type1_error_rate <- colMeans(type1_df)
  mean_model_se <- colMeans(se_df) # Mean of SEs estimated by the model
  
  # Calculate relative bias (proportion and percentage) - quantifies the average deviation of your 
  # estimated parameter from the true parameter as a proportion of the true parameter.
  # If on average, the relative bias is -10%, the estimator underestimates the true value by 10%.
  relative_bias <- (mean_rd - true_rd) / true_rd
  relative_bias_pct <- relative_bias * 100
  
  # Calculate the mean squared error which quantifies the average squared difference between the 
  # estimates produced by your model (or simulation) and the true parameter value.
  # A lower MSE indicates that the estimator is, on average, closer to the true parameter value, combining both its bias and variability.

  # Per-arm MSE: (Bias^2 + Variance) for each column of rd_df
  mse_per_arm <- mapply(function(est, tr) {
    b <- mean(est) - tr         # bias for this arm
    v <- var(est)               # variance across simulations
    b^2 + v
  }, as.data.frame(rd_df), true_rd)
  
  mean_risk <- colMeans(risk_df)
  
  # Ratio of model SE error vs empirical SE
  relative_error_se <- (mean_model_se - empirical_se) / empirical_se
  relative_error_pct_se <- relative_error_se * 100
  
  summary_stats <- data.frame(
    treatment = paste0("Treatment", 0:(n_arms - 1)),
    mean_risk = mean_risk,
    mean_rd = c(NA, mean_rd),              # NA for control
    mean_lcl = c(NA, mean_lcl),
    mean_ucl = c(NA, mean_ucl),
    mean_ci_width = c(NA, mean_ci_width),
    power = c(NA, power),
    type1_error_rate = c(NA, type1_error_rate),
    relative_bias = c(NA, relative_bias),
    relative_bias_pct = c(NA, relative_bias_pct),
    mse_rd = c(NA, mse_per_arm)
    )
  
  # Add coverage to summary stats
  summary_stats$coverage <- c(NA, coverage_by_arm)
  summary_stats$empirical_se <- c(NA, empirical_se)  # NA for control
  summary_stats$mean_model_se <- c(NA, mean_model_se)
  summary_stats$relative_error_se <- c(NA, relative_error_se)
  summary_stats$relative_error_pct_se <- c(NA, relative_error_pct_se)
  
  list(
    simulation_results = list(
      rd = rd_df,
      lcl = lcl_df,
      ucl = ucl_df,
      ci_width = ci_width_df,
      se = se_df,
      success = success_df,
      type1_error = type1_df),
    
    recruitment_details = recruitment_details,
    trial_data = trial_data_list,
    trial_data_inf = trial_data_inf_list,
    summary_stats = summary_stats,
    long_summary = long_summary,
    site_accrual = site_accrual_list,
    
    # Whole-trial level summaries
    trial_level_summary = data.frame(
      sim = 1:nSims,
      total_patients = total_patients,
      control_patients = control_patients,
      total_sites = total_sites,
      total_months = total_months,
      proportion_successful_arms = arm_success_rate,
      proportion_error_arms = arm_error_rate
    ),
    
    # Summary values
    trial_level_means = data.frame(
      mean_total_patients       = mean(total_patients),
      median_total_patients     = median(total_patients),
      q1_total_patients         = quantile(total_patients, 0.25),
      q3_total_patients         = quantile(total_patients, 0.75),
      
      mean_control_patients     = mean(control_patients),
      median_control_patients   = median(control_patients),
      q1_control_patients       = quantile(control_patients, 0.25),
      q3_control_patients       = quantile(control_patients, 0.75),
      
      mean_total_sites          = mean(total_sites),
      median_total_sites        = median(total_sites),
      q1_total_sites            = quantile(total_sites, 0.25),
      q3_total_sites            = quantile(total_sites, 0.75),
      
      mean_total_months         = mean(total_months),
      median_total_months       = median(total_months),
      q1_total_months           = quantile(total_months, 0.25),
      q3_total_months           = quantile(total_months, 0.75),
      
      mean_overall_platform_power = mean_platform_power,
      overall_platform_power_any  = overall_platform_power_any,
      mean_overall_error_rate     = mean_platform_errorrate,
      overall_platform_type1_any  = overall_platform_type1_any,
      
      mse_rd = mean(mse_per_arm),
      relative_bias = mean(relative_bias),
      relative_bias_pct = mean(abs(relative_bias_pct)),
      coverage_prob = mean(coverage_by_arm),
      platform_coverage_prob = platform_coverage,
      mean_overall_empirical_se = overall_empirical_se,
      mean_overall_model_se     = mean(mean_model_se),
      mean_relative_error_se    = mean(relative_error_se),
      mean_relative_error_pct_se = mean(abs(relative_error_pct_se))
    ))
}

######################################################

################################################################################################################
######################### R script to simulate multiple Staggered Platform Trials with 5 arms ##################
################################################################################################################
# This code will simulate multiple 5-arm non-inferiority randomised controlled platform trials and aggregate
# the summary statistics and save the lists

simulate_trials_spt <- function(nSims, 
                                n_arms,
                                patients_per_experimental_arm,
                                death_probs,
                                death_probs_inferior,
                                death_alt,
                                ni_margin,
                                stagger_months,
                                mean_sites_per_month,
                                max_sites,
                                mean_patients_per_site_per_year,
                                allocation_mode = c("fixed", "sqrtk"),
                                allocation_ratio = NULL,
                                true_rd, 
                                expected_control_arm
) {
  
  allocation_mode <- match.arg(allocation_mode)
  
  # -----------------------------------------------------
  # 1. Initialize storage for trial‐level summaries
  # -----------------------------------------------------
  trialnum            <- numeric(nSims)
  total_sites         <- numeric(nSims)
  total_patients      <- numeric(nSims)
  control_patients    <- numeric(nSims)
  total_months        <- numeric(nSims)
  recruitment_details <- vector("list", nSims)
  trial_data_list     <- vector("list", nSims)
  trial_data_inf_list <- vector("list", nSims)
  site_accrual_list <- vector("list", nSims)
  
  # -----------------------------------------------------
  # 2. Preallocate lists for per‐arm, per‐simulation outcomes
  #    Each element is a numeric vector of length (n_arms - 1)
  # -----------------------------------------------------
  rd_list            <- vector("list", nSims)
  lcl_list           <- vector("list", nSims)
  ucl_list           <- vector("list", nSims)
  ci_width_list      <- vector("list", nSims)
  se_list            <- vector("list", nSims)
  success_list       <- vector("list", nSims)
  type1_error_list   <- vector("list", nSims)
  
  # -----------------------------------------------------
  # 3. Preallocate list for risk_by_arm (each element = small data.frame)
  # -----------------------------------------------------
  risk_by_arm_all <- vector("list", nSims)
  
  # -----------------------------------------------------
  # 4. Main simulation loop (using your exact call)
  # -----------------------------------------------------
  for (i in seq_len(nSims)) {
    sim <- simulate_full_trial_spt(
      n_arms                         = n_arms,
      patients_per_experimental_arm  = patients_per_experimental_arm,
      death_probs                    = death_probs,
      death_probs_inferior           = death_probs_inferior,
      death_alt                      = death_alt,
      ni_margin                      = ni_margin,
      stagger_months                 = stagger_months,
      mean_sites_per_month           = mean_sites_per_month,
      max_sites                      = max_sites,
      mean_patients_per_site_per_year = mean_patients_per_site_per_year,
      allocation_mode                = allocation_mode,
      allocation_ratio               = allocation_ratio
    )
    
    # 4a. Trial‐level identifiers
    trialnum[i]       <- i
    total_sites[i]    <- sim$recruitment$total_sites
    total_patients[i] <- sim$recruitment$total_patients
    control_patients[i] <- sum(sim$outcomes$trial_data$arm == 0)
    total_months[i]   <- sim$recruitment$total_months
    
    # 4b. Store recruitment details, site-accrual and patient‐level data
    recruitment_details[[i]] <- sim$recruitment$recruitment_details
    trial_data_list[[i]]     <- sim$outcomes$trial_data
    trial_data_inf_list[[i]] <- sim$type1_error$trialdata_inf
    site_accrual_list[[i]] <- sim$site_accrual
    
    # 4c. Store per‐arm RD / CI / SE / success (length = n_arms - 1)
    rd_list[[i]]       <- sim$outcomes$rd
    lcl_list[[i]]      <- sim$outcomes$lcl
    ucl_list[[i]]      <- sim$outcomes$ucl
    se_list[[i]]       <- sim$outcomes$se
    ci_width_list[[i]] <- sim$outcomes$ucl - sim$outcomes$lcl
    success_list[[i]]  <- sim$outcomes$success
    
    # 4d. Store type I error vector (length = n_arms - 1)
    type1_error_list[[i]] <- sim$type1_error$type1_error
    
    # 4e. Build a small data.frame of risk_by_arm for this simulation
    this_risk_df            <- do.call(rbind, sim$outcomes$risk_by_arm)
    this_risk_df$Comparison <- rep(sim$outcomes$comparisons, each = 2)
    this_risk_df$sim        <- i
    risk_by_arm_all[[i]]    <- this_risk_df
    
  }
  
  # -----------------------------------------------------
  # 5. Convert lists of vectors into matrices/data.frames
  # -----------------------------------------------------
  rd_df         <- do.call(rbind, rd_list)            # nSims × (n_arms - 1)
  lcl_df        <- do.call(rbind, lcl_list)
  ucl_df        <- do.call(rbind, ucl_list)
  ci_width_df   <- do.call(rbind, ci_width_list)
  se_df         <- do.call(rbind, se_list)
  success_df    <- do.call(rbind, success_list)
  type1_df      <- do.call(rbind, type1_error_list) 
  
  risk_by_arm_full <- do.call(rbind, risk_by_arm_all)
  risk_by_arm_full <- risk_by_arm_full[, c("sim", "Comparison", "Arm", "N", "Predicted_Prob")]
  
  # Create a summary of mean predicted risk and mean N per arm per comparison
  risk_by_arm_summary <- aggregate(
    cbind(Predicted_Prob, N) ~ Arm + Comparison,
    data = risk_by_arm_full,
    FUN = mean
  )
  
  # Step 5b: Identify control arms with < patients_per_experimental_arm
  risk_by_arm_full$low_control_flag <- ifelse(
    risk_by_arm_full$Arm == "Control" & 
      risk_by_arm_full$N < expected_control_arm,
    1, 0
  )
  
  # Summarize number and proportion of low-control simulations per comparison
  low_control_summary <- aggregate(
    low_control_flag ~ Comparison,
    data = risk_by_arm_full,
    FUN = function(x) c(count = sum(x), prop = mean(x))
  )
  
  # Flatten the matrix-style result into a data frame
  low_control_summary <- do.call(data.frame, low_control_summary)
  colnames(low_control_summary)[2:3] <- c("control_sample<expected", "prop_low_control")
  
  # Compute mean N when control was low
  control_subset <- subset(risk_by_arm_full, Arm == "Control" & N < expected_control_arm)
  
  if (nrow(control_subset) > 0) {
    mean_low_control_N <- aggregate(
      N ~ Comparison,
      data = control_subset,
      FUN = mean
    )
    colnames(mean_low_control_N)[2] <- "mean_N_if_low_control"
  } else {
    mean_low_control_N <- data.frame(
      Comparison = unique(risk_by_arm_full$Comparison),
      mean_N_if_low_control = "none"  # Replace NA with "none"
    )
  }
  
  # Merge into summary
  low_control_summary <- merge(
    low_control_summary, 
    mean_low_control_N, 
    by = "Comparison", 
    all.x = TRUE
  )
  
  # Step 5c: Calculate per-simulation/overall power and type 1 error rate
  
  arm_success_rate <- rowMeans(success_df)
  mean_platform_power <- mean(arm_success_rate)
  
  # Platform-level power (at least one effective arm declared efficacious)
  overall_platform_power_any <- mean(apply(success_df, 1, function(x) any(x == 1)))
  
  arm_error_rate <- rowMeans(type1_df)  # gives a vector of length nSims
  mean_platform_errorrate <- mean(arm_error_rate)
 
  # Platform-level error rate (at least one ineffective arm incorrectly declared efficacious)
  overall_platform_type1_any <- mean(apply(type1_df, 1, function(x) any(x == 1)))
  
  # Step 5d: Calculate other performance measures
  
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
  
  # -----------------------------------------------------
  # 6. Build “long_summary”: one row per sim × arm
  # -----------------------------------------------------
  arm_labels_full <- paste0("Treatment", 0:(n_arms - 1))
  long_summary <- data.frame(
    sim          = rep(1:nSims, each = n_arms),
    arm          = rep(arm_labels_full, times = nSims),
    rd           = NA,
    lcl          = NA,
    ucl          = NA,
    ci_width     = NA,
    success      = NA,
    type1error = NA  
  )
  
  # Fill RD/CI/success for experimental arms
  for (j in seq_len(n_arms - 1)) {
    idx   <- which(long_summary$arm == paste0("Treatment", j))
    long_summary$rd[idx]       <- rd_df[, j]
    long_summary$lcl[idx]      <- lcl_df[, j]
    long_summary$ucl[idx]      <- ucl_df[, j]
    long_summary$ci_width[idx] <- ci_width_df[, j]
    long_summary$success[idx]  <- success_df[, j]
    long_summary$type1error[idx]  <- type1_df[, j]
  }
  
  # -----------------------------------------------------
  # 7. Compute per‐arm summary statistics
  # -----------------------------------------------------
  
  mean_rd        <- colMeans(rd_df)
  mean_lcl       <- colMeans(lcl_df)
  mean_ucl       <- colMeans(ucl_df)
  mean_ci_width <- colMeans(ci_width_df)
  mean_se        <- colMeans(se_df)
  mean_power     <- colMeans(success_df)
  mean_type1     <- colMeans(type1_df)
  mean_model_se <- colMeans(se_df) # Mean of SEs estimated by the model
  
  # 7b. Calculate further performance calculations
  
  # Calculate relative bias (proportion and percentage) - quantifies the average deviation of your 
  # estimated parameter from the true parameter as a proportion of the true parameter.
  # If on average, the relative bias is -10%, the estimator underestimates the true value by 10%.
  bias_per_arm            <- mean_rd - true_rd
  relative_bias_per_arm   <- bias_per_arm / true_rd
  relative_bias_pct_per_arm <- 100 * relative_bias_per_arm
  
  # Calculate the mean squared error which quantifies the average squared difference between the 
  # estimates produced by your model (or simulation) and the true parameter value.
  # A lower MSE indicates that the estimator is, on average, closer to the true parameter value, combining both its bias and variability.
  
  # Per-arm MSE: (Bias^2 + Variance) for each column of rd_df
  mse_per_arm <- mapply(function(est, tr) {
    b <- mean(est) - tr         # bias for this arm
    v <- var(est)               # variance across simulations
    b^2 + v
  }, as.data.frame(rd_df), true_rd)
  
  # Ratio of model SE error vs empirical SE
  relative_error_se <- (mean_model_se - empirical_se) / empirical_se
  relative_error_pct_se <- relative_error_se * 100
  
  # -----------------------------------------------------
  # 8. Build `summary_stats` (one row per arm)
  # -----------------------------------------------------
  summary_stats <- data.frame(
    treatment         = arm_labels_full,
    mean_rd           = c(NA, mean_rd),
    mean_lcl          = c(NA, mean_lcl),
    mean_ucl          = c(NA, mean_ucl),
    mean_ci_width = c(NA, mean_ci_width),
    mean_se           = c(NA, mean_se),
    power             = c(NA, mean_power),
    type1_error_rate  = c(NA, mean_type1),
    mse_rd            = c(NA, mse_per_arm),
    relative_bias        = c(NA, relative_bias_per_arm),
    relative_bias_pct    = c(NA, relative_bias_pct_per_arm)
  )
  
  # Add coverage to summary stats
  summary_stats$coverage <- c(NA, coverage_by_arm)
  summary_stats$empirical_se <- c(NA, empirical_se)  # NA for control
  summary_stats$relative_error_se <- c(NA, relative_error_se)
  summary_stats$relative_error_pct_se <- c(NA, relative_error_pct_se)
  summary_stats$mean_model_se <- c(NA, mean_model_se)
  
  # -----------------------------------------------------
  # 9. Build a trial‐level summary (one row per simulation)
  # -----------------------------------------------------
  trial_level_summary <- data.frame(
    sim                  = 1:nSims,
    total_patients       = total_patients,
    control_patients     = control_patients,
    total_sites          = total_sites,
    total_months         = total_months,
    proportion_successful_arms = arm_success_rate,
    proportion_error_arms = arm_error_rate)
  
  # -----------------------------------------------------
  # 10. Aggregated summaries across all simulations
  # -----------------------------------------------------
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
    relative_bias = mean(relative_bias_per_arm),
    relative_bias_pct = mean(abs(relative_bias_pct_per_arm)),
    coverage_prob = mean(coverage_by_arm),
    platform_coverage_prob = platform_coverage,
    mean_overall_empirical_se = overall_empirical_se,
    mean_overall_model_se     = mean(mean_model_se),
    mean_relative_error_se    = mean(relative_error_se),
    mean_relative_error_pct_se = mean(abs(relative_error_pct_se))
  )
  
  # -----------------------------------------------------
  # 11. Return all results
  # -----------------------------------------------------
  return(list(
    simulation_results = list(
      rd         = rd_df,
      lcl        = lcl_df,
      ucl        = ucl_df,
      ci_width = ci_width_df,
      se         = se_df,
      success    = success_df,
      type1_error = type1_df),
    
    recruitment_details = recruitment_details,
    trial_data          = trial_data_list,
    trial_data_inf      = trial_data_inf_list,
    risk_by_arm_full    = risk_by_arm_full,
    risk_by_arm_summary = risk_by_arm_summary,
    low_control_summary = low_control_summary,
    summary_stats       = summary_stats,
    long_summary        = long_summary,
    trial_level_summary = trial_level_summary,
    trial_level_means   = trial_level_means,
    site_accrual        = site_accrual_list
  ))
}

################################################################################################################


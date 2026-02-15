################################################################################################################
####################### R script Simulation wrapper to run multiple simulations and aggregate results###########
################################################################################################################

simulate_trials_and_results_multi <- function(nSims, 
                                        nPatients, 
                                        death0, 
                                        death1, 
                                        ni_margin,
                                       mean_sites_per_month, 
                                       max_sites, 
                                       mean_patients_per_site_per_year,
                                       death_alt,
                                       starts_months = starts_months
) {
  
  # 1) Run shared-calendar recruitment once for all sims/trials
  recr_all <- simulate_4_staggered_competing(
    nSims                 = nSims,
    nPatients             = nPatients,
    mean_sites_per_month  = mean_sites_per_month,
    max_sites             = max_sites,
    mean_patients_per_site_per_year = mean_patients_per_site_per_year,
    starts_months         = starts_months,
    nTrials = nTrials,
    site_concurrency_cap = site_concurrency_cap,  # max concurrent trials per site
    global_nsites = max_sites                     # size of the global site universe
  )
  
  program_summary <- do.call(rbind, lapply(recr_all$sims, function(x) x$program_summary))
  
  # 2) Pre-allocate per (Sim x Trial) storage
  nRows <- nSims * nTrials
  Sim <- integer(nRows)
  Trial <- character(nRows)
  rd <- numeric(nRows); lcl <- numeric(nRows); ucl <- numeric(nRows)
  se_vec <- numeric(nRows); ci_width <- numeric(nRows); success <- numeric(nRows)
  risk_treatment <- numeric(nRows); risk_control <- numeric(nRows)
  total_sites <- integer(nRows); total_patients <- integer(nRows)
  control_patients <- integer(nRows); total_months <- integer(nRows)
  type1_error <- numeric(nRows)
  Start_Month <- integer(nRows);  End_Month   <- integer(nRows)
  
  recruitment_details_list <- vector("list", nRows)
  trial_data_list <- vector("list", nRows)
  trial_data_inf_list <- vector("list", nRows)
  site_accrual_list <- vector("list", nRows)
  
  # 3) For each Sim x trial: outcomes + type I error (single-trial)
  # Simulates trial 1 x sims for outcomes and type I error then moves onto the next trial
  
  r <- 1L
  for (s in seq_len(nSims)) {                  # Outer loop: simulation replicate
    per_trial <- recr_all$sims[[s]]$per_trial  # recruitment list for all trials in sim x: "Trial_1", "Trial_2", ...
    
    for (t in seq_len(nTrials)) {              # Inner loop: trial index within that replicate
      tr_name <- paste0("Trial_", t)
      rec_obj <- per_trial[[tr_name]]          # recruitment log for one trial
      
      # Simulate outcomes and type I error
      outcomes <- simulate_trial_outcomes(nPatients, death0, death1, ni_margin)
      type1    <- simulate_type1_error(outcomes$trialdata, death0, death_alt, ni_margin)
      
      # Fill row r
      Sim[r] <- s
      Trial[r] <- tr_name
      
      rd[r] <- outcomes$rd
      lcl[r] <- outcomes$lcl
      ucl[r] <- outcomes$ucl
      ci_width[r] <- ucl[r] - lcl[r]
      se_vec[r] <- outcomes$se
      success[r] <- outcomes$success
      risk_control[r] <- outcomes$risk_control
      risk_treatment[r] <- outcomes$risk_treatment
      
      total_sites[r] <- rec_obj$total_sites
      total_patients[r] <- nPatients
      control_patients[r] <- sum(outcomes$trialdata$treatment == 0)
      total_months[r] <- rec_obj$total_months
      type1_error[r] <- type1$type1_error
      
      # Carry logs with Sim/Trial tags
      if (!is.null(rec_obj$recruitment_details) && nrow(rec_obj$recruitment_details)) {
        tmp_rd <- rec_obj$recruitment_details
        tmp_rd$Sim <- s; tmp_rd$Trial <- tr_name
        recruitment_details_list[[r]] <- tmp_rd
      } else {
        recruitment_details_list[[r]] <- data.frame(Sim=integer(0), Trial=character(0))
      }
      
      if (!is.null(rec_obj$site_accrual) && nrow(rec_obj$site_accrual)) {
        tmp_sa <- rec_obj$site_accrual
        tmp_sa$Sim <- s; tmp_sa$Trial <- tr_name
        site_accrual_list[[r]] <- tmp_sa
      } else {
        site_accrual_list[[r]] <- data.frame(Sim=integer(0), Trial=character(0))
      }
      
      td <- outcomes$trialdata; td$Sim <- s; td$Trial <- tr_name
      trial_data_list[[r]] <- td
      
      tdi <- type1$trialdata_inf; tdi$Sim <- s; tdi$Trial <- tr_name
      trial_data_inf_list[[r]] <- tdi
      
      r <- r + 1L
    }
  }
  
  # --- 4) Bind per (Sim x Trial) results ---
  
  simulation_results <- data.frame(
    Sim, Trial,
    risk_treatment, risk_control, rd, lcl, ucl, ci_width, se_vec, success,
    total_months, total_sites, total_patients, control_patients, type1_error,
    row.names = NULL
  )
  
  # Draw start and end month for each trial and simulation
  se_tab <- recr_all$simulation_results[, c("Sim","Trial","Start_Month","End_Month")]
  
  # Merge into your main trial-level results
  simulation_results <- merge(
    simulation_results, 
    se_tab,
    by = c("Sim", "Trial"),
    all.x = TRUE
  )
  
  # --- 5) Summary stats across all Sim x Trial rows ---
  true_rd <- death1 - death0
  mean_risk_control <- mean(simulation_results$risk_control)
  mean_risk_treatment <- mean(simulation_results$risk_treatment)
  mean_rd <- mean(simulation_results$rd)
  mean_lcl <- mean(simulation_results$lcl)
  mean_ucl <- mean(simulation_results$ucl)
  mean_ci_width <- mean(simulation_results$ci_width)
  avg_sites = mean(simulation_results$total_sites)
  power <- mean(simulation_results$success)          # Proportion of simulations with success = 1
  type1_error_rate <- mean(simulation_results$type1_error)  # Proportion with type 1 error
  average_total_patients <- mean(simulation_results$total_patients)
  median_total_patients <- median(simulation_results$total_patients)
  average_control_patients <- mean(simulation_results$control_patients)
  median_control_patients <- median(simulation_results$control_patients)
  mean_program_duration <- mean(program_summary$Program_Duration_Months)
  median_program_duration <- median(program_summary$Program_Duration_Months)
  
  # Calculate relative bias (proportion and percentage) - quantifies the average deviation of your 
  # estimated parameter from the true parameter as a proportion of the true parameter.
  # If on average, the relative bias is -10%, the estimator underestimates the true value by 10%.
  relative_bias <- (mean_rd - true_rd) / true_rd
  relative_bias_pct <- relative_bias * 100
  
  # Calculate the mean squared error which quantifies the average squared difference between the 
  # estimates produced by your model (or simulation) and the true parameter value.
  # A lower MSE indicates that the estimator is, on average, closer to the true parameter value, combining both its bias and variability.
  mse_rd <- mean((simulation_results$rd - true_rd)^2)
  
  # Coverage probability is the proportion of simulations where the true risk difference lies between the lower and upper confidence limits
  coverage <- mean((simulation_results$lcl <= true_rd) & (simulation_results$ucl >= true_rd)) * 100
  
  # Calculate the Empirical Standard Error (ESE) for the risk difference estimates
  # The empirical standard error reflects the actual variability of your estimator,
  # you can compare it to the model-based standard error to check how well your model’s 
  # theoretical uncertainty approximates the actual variability seen in the simulation.
  empirical_se <- sd(simulation_results$rd)
  
  # Calculate the Mean Model-Based Standard Error
  mean_model_se <- mean(simulation_results$se_vec)
  
  # Compute Relative Error in Model-Based SE:A normalized measure to show the discrepancy between 
  # the average model-based uncertainty and the empirical variability.This value tells you the proportion 
  # by which the model-based SE overestimates (if positive) or underestimates (if negative) the actual observed variability.
  relative_error <- (mean_model_se - empirical_se) / empirical_se
  relative_error_pct <- relative_error * 100
  
  summary_stats <- data.frame(
    mean_risk_control = mean_risk_control,
    mean_risk_treatment = mean_risk_treatment,
    mean_rd = mean_rd,
    mean_lcl = mean_lcl,
    mean_ucl = mean_ucl,
    mean_ci_width = mean_ci_width,
    power = power,
    type1_error_rate = type1_error_rate,
    mean_program_duration = mean_program_duration,
    median_program_duration = median_program_duration,
    average_total_patients = average_total_patients,
    median_total_patients = median_total_patients,
    average_control_patients = average_control_patients,
    median_control_patients = median_control_patients,
    avg_sites = avg_sites,
    relative_bias = relative_bias,
    relative_bias_pct = relative_bias_pct,
    mse_rd = mse_rd,
    coverage_prob = coverage,
    empirical_se = empirical_se,
    mean_model_se = mean_model_se,
    relative_error_model_se = relative_error,
    relative_error_model_se_pct = relative_error_pct
    )  

  # --- 6) Return combined outputs ---
  list(
    simulation_results = simulation_results,
    recruitment_details = do.call(rbind, recruitment_details_list),     # combined monthly logs
    trial_data = trial_data_list,
    trial_data_inf = trial_data_inf_list,
    summary_stats = summary_stats,
    site_accrual = do.call(rbind, site_accrual_list),
    program_summary = program_summary
  )
}

######################################################

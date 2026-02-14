################################################################################################################
####R Script for simulating recruitment for a Staggered Platform Trial with 5 arms and a non-concurrent control#
################################################################################################################
# This code will simulate recruitment for a 5-arm non-inferiority randomised controlled platform trial using 
# non-concurrent control arm which will continue to enrol until it reaches its required per arm sample size
# versus treatment 1 versus treatment 2 versus treatment 3 and treatment 4 with a non-concurrent control
# This function will estimate the total duration of the trial based on the number of sites activated per month and
# the number of patients recruited per year

simulate_recruitment_sptncc <- function(n_arms,
                                     patients_per_experimental_arm,
                                     mean_sites_per_month,
                                     max_sites,
                                     mean_patients_per_site_per_year,
                                     stagger_months) {
  
  # Arm names from 0 to (n_arms - 1)
  arm_names <- paste0("Arm", 0:(n_arms - 1))
  
  # Arm 0 (control) and Arm 1 open at month 0, others are staggered
  arm_open_month <- c(0, 0, stagger_months * (1:(n_arms - 2)))
  
  # Recruitment targets: all arms capped equally
  arm_target <- rep(patients_per_experimental_arm, n_arms)

  arm_recruited <- rep(0, n_arms)
  
  # Initialize trial state
  sites <- 0
  month <- 0
  
  # Trackers
  patients_per_month <- c()
  sites_per_month <- c()
  month_numbers <- c()
  recruited_by_arm <- list()
  
  # Initialize site-level accrual data frame
  site_accrual <- data.frame()
  
  # Loop until all experimental arms are full
  while (any(arm_recruited[-1] < patients_per_experimental_arm)) {
    
    # Add new sites
    if (sites < max_sites) {
      new_sites <- rpois(1, mean_sites_per_month)
      sites <- min(sites + new_sites, max_sites)
    }
    
    # Simulate monthly patient recruitment
    patients_by_site <- rpois(sites, mean_patients_per_site_per_year / 12)
    monthly_patients <- sum(patients_by_site)
    
    # Determine open arms (must be open and not full)
    open_arms <- which(arm_open_month <= month & arm_recruited < arm_target)
    remaining_capacity <- arm_target[open_arms] - arm_recruited[open_arms]
    total_remaining_capacity <- sum(remaining_capacity)
    
    # Cap total recruitment if needed
    if (monthly_patients > total_remaining_capacity) {
      monthly_patients <- total_remaining_capacity
    }
    
    # Allocate patients
    allocation <- rep(0, n_arms)
    
    if (length(open_arms) > 0 && monthly_patients > 0) {
      alloc_probs <- rep(1, length(open_arms)) / length(open_arms)
      raw_alloc <- rmultinom(1, monthly_patients, prob = alloc_probs)
      capped_alloc <- pmin(as.vector(raw_alloc), remaining_capacity)
      
      # Redistribute if shortfall
      shortfall <- monthly_patients - sum(capped_alloc)
      if (shortfall > 0) {
        still_open <- which(remaining_capacity - capped_alloc > 0)
        if (length(still_open) > 0) {
          redistrib <- rmultinom(1, shortfall, prob = alloc_probs[still_open])
          capped_alloc[still_open] <- capped_alloc[still_open] + as.vector(redistrib)
        }
      }
      
      allocation[open_arms] <- capped_alloc
      arm_recruited <- arm_recruited + allocation
    }
    
    # Track site-level accrual
    site_df <- data.frame(
      Month = rep(month, sites),
      Site_ID = seq_len(sites),
      Patients = patients_by_site
    )
    site_accrual <- rbind(site_accrual, site_df)
    
    # Track data
    patients_per_month <- c(patients_per_month, monthly_patients)
    sites_per_month <- c(sites_per_month, sites)
    month_numbers <- c(month_numbers, month)
    recruited_by_arm[[length(recruited_by_arm) + 1L]] <- allocation
    
    # Increment month ---
    month <- month + 1L
  }
  
  # Format recruitment data
  recruited_df <- do.call(rbind, recruited_by_arm)
  colnames(recruited_df) <- arm_names
  
  recruitment_details <- data.frame(
    Month = month_numbers,
    Active_Sites = sites_per_month,
    Total_Recruited = patients_per_month
  )
  recruitment_details <- cbind(recruitment_details, recruited_df)
  recruitment_details$Cumulative_Recruited <- cumsum(recruitment_details$Total_Recruited)
  
  for (i in seq_len(n_arms)) {
    col <- arm_names[i]
    recruitment_details[[paste0("Cumulative_", col)]] <- cumsum(recruitment_details[[col]])
  }
  
  list(
    recruitment_details = recruitment_details,
    total_sites = sites,
    total_patients = sum(arm_recruited),
    total_months = month,
    patients_per_arm = arm_recruited,
    site_accrual = site_accrual
  )
}

##############################################################

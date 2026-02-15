################################################################################################################
##### R Script for simulating recruitment for a Staggered Platform Trial with 5 arms and a concurrent control####
################################################################################################################
# This code will simulate recruitment for a 5-arm non-inferiority randomised controlled platform trial using 
# three options of random allocation to participants with a # 1:1:1:1:1, 2:1:1:1:1 or √k:1:1:1:1 allocation
# to treatment 0 ("best available therapy") versus treatment 1 versus treatment 2 versus treatment 3 and 
# treatment 4 dependent on the arms that are actively recruiting
# This fuinction will estimate the total duration of the trial based on the number of sites activated per month and
# the number of patients recruited per year

simulate_recruitment_spt <- function(n_arms,
                                     patients_per_experimental_arm,
                                     mean_sites_per_month,
                                     max_sites,
                                     mean_patients_per_site_per_year,
                                     stagger_months,
                                     allocation_mode = c("fixed", "sqrtk"),
                                     allocation_ratio = NULL) {
  
  allocation_mode <- match.arg(allocation_mode)
  
  # ---- Arm labels ----
  arm_names <- paste0("Arm", 0:(n_arms - 1))
  
  # ---- Opening schedule: Arm0=control, Arm1 open at month 0; others staggered ----
  arm_open_month <- c(0, 0, stagger_months * (1:(n_arms - 2)))
  
  # Recruitment targets: control (Arm0) uncapped, others capped
  arm_target <- rep(NA_real_, n_arms)                # NA = uncapped
  arm_target[-1] <- patients_per_experimental_arm    # arms 2..n capped; arm 1 (index 1) is control and NA
  
  arm_recruited <- rep(0, n_arms)
  
  # Initialize trial state
  sites <- 0
  month <- 0
  
  # Trackers
  patients_per_month <- integer(0)
  sites_per_month <- integer(0)
  month_numbers <- integer(0)
  recruited_by_arm <- list()
  site_accrual <- data.frame()

  control_arm   <- 1L  # Arm0 = index 1
  
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
    
    # Determine open arms (must be open and not full) with uncapped control
    open_by_time  <- which(arm_open_month <= month)
    # experimental arms that still have capacity (ignore control NA in arm_target)
    exp_can_recruit <- which((seq_len(n_arms) != control_arm) & (arm_recruited < arm_target))
    # final open arms = control if open-by-time + experimental with capacity and open-by-time
    open_arms <- union(intersect(open_by_time, control_arm),
                       intersect(open_by_time, exp_can_recruit))
    
    # Remaining capacity; treat control as Inf for capping ONLY
    remaining_capacity <- arm_target[open_arms] - arm_recruited[open_arms]
    if (length(open_arms) > 0L) {
      is_control_in_open <- open_arms == control_arm
      if (any(is_control_in_open)) {
        remaining_capacity[is_control_in_open] <- Inf
      }
    }
    total_remaining_capacity <- sum(remaining_capacity)
    
    # Cap total recruitment if needed
    if (monthly_patients > total_remaining_capacity) {
      monthly_patients <- total_remaining_capacity
    }
    
    # Allocate patients
    allocation <- rep(0, n_arms)
    
    if (length(open_arms) > 0 && monthly_patients > 0) {
      
      if (allocation_mode == "fixed") {
        if (is.null(allocation_ratio)) {
          stop("allocation_ratio must be provided when allocation_mode is 'fixed'")
        }
        alloc_weights <- allocation_ratio[open_arms]
        alloc_probs <- alloc_weights / sum(alloc_weights)
      } 
      else if (allocation_mode == "sqrtk") {
        control_arm <- 1L  # Arm0 = index 1
        control_in_open <- control_arm %in% open_arms
        experimental_open <- setdiff(open_arms, control_arm)
        k <- length(experimental_open)
        
        alloc_weights <- numeric(length(open_arms))
        names(alloc_weights) <- arm_names[open_arms]
        
        if (control_in_open) {
          alloc_weights[as.character(arm_names[control_arm])] <- if (k > 0) sqrt(k) else 1
        }
        if (k > 0) {
          alloc_weights[as.character(arm_names[experimental_open])] <- 1
        }
        
        total_weight <- sum(alloc_weights)
        
        if (total_weight == 0 || any(is.na(alloc_weights))) {
          warning(paste("Allocation fallback to equal allocation at month", month))
          alloc_probs <- rep(1 / length(open_arms), length(open_arms))
        } else {
          alloc_probs <- alloc_weights / total_weight
        }
      }
      
      raw_alloc <- as.vector(rmultinom(1, monthly_patients, prob = alloc_probs))
      capped_alloc <- pmin(raw_alloc, remaining_capacity)
      
      if (sum(capped_alloc) < monthly_patients) {
        shortfall <- monthly_patients - sum(capped_alloc)
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
    recruited_by_arm[[length(recruited_by_arm) + 1L]] <- allocation  # don't index by month (month can be 0)
    
    # Increment month at the *end*
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
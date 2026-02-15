###############################################################################################
##R script for simulating trial recruitment for 4 stand-alone 2-arm non-inferiority trial######
###############################################################################################

# This code will simulate recruitment for four staggered stand-alone 2-arm randomised controlled trials 
# This simulation will estimate the total duration of the trial based on the number of sites activated per month, 
# the number of patients recruited per site per year, the maximum sites for all trials and sites can only run a maximum
# of 2 trials at one time. If only one trial is running it gets the full allocation of patients.
# If 2 trials overlap, each gets ½ the usual per-site rate.

# Step 1 Initialise a single-trial recruitment state object that captures the 
# trials recruitment from day zero. This sets up the variables and data frames 
# needed to tract site activation and patient accrual over time.
trct_init <- function(nPatients,
                      mean_sites_per_month,
                      max_sites,
                      mean_patients_per_site_per_year) {
  list(
    nPatients = nPatients,           # Total number of patients needed for this trial
    sites = 0L,                      # Site count
    sites_vec = integer(0),          # GLOBAL site IDs assigned to this trial
    patients = 0L,                   # Patient count
    month = 0L,                      # Month count
    done = FALSE,                    # Initialised as FALSE and becomes TRUE once recruitment is complete
    mean_sites_per_month = mean_sites_per_month,   # Average number of new sites activated per month
    max_sites = max_sites,                         # Max sites for each trial and all trials
    base_rate_month = mean_patients_per_site_per_year / 12,    # Monthly rate of patients per site
    patients_per_month = integer(0),               # Patient count
    sites_per_month = integer(0),                  # Site count
    month_numbers = integer(0),                    # Month by month identifiers
    site_accrual = data.frame()                    # Site accrual data
  )
}

# Step 2: Advance recruitment of a trial by ONE month
# Multiplies the per-site rate by `site_weights`, and log a calendar month.
# - site_weights: A vector specifying the recruitment rate for each site for this month.
# If a site participates in two trials, the allocation is reduced to 0.5 and shared between the two sites
# - current_global_month: the current calendar month in the shared simulation timeline to align all staggered trials
# start month
# Returns the updated state list for that trial with month, patient counts and logs until done

trct_step <- function(state, site_weights = NULL, current_global_month) {
  if (state$done) return(state)
  
  # Advance the trial’s local month counter
  state$month <- state$month + 1L
  
  # Use the orchestrator’s calendar label directly for all trials
  calendar_month <- as.integer(current_global_month)
  
  # number of currently assigned (GLOBAL) sites
  n_sites <- length(state$sites_vec)
  state$sites <- as.integer(n_sites)
  
  ## --- Per-site accrual with site-specific weights ---
  # Draws number of patients per site based on the number of active trials
  
  if (n_sites > 0L) {
    # default to equal weight 1 if not provided / wrong length
    if (is.null(site_weights) || length(site_weights) != n_sites) {
      site_weights <- rep(1, n_sites)
    }
    # Each site gets its own λ depending on concurrency
    lambda_vec <- state$base_rate_month * site_weights
    
    patients_by_site <- rpois(n_sites, lambda_vec)
  } else {
    patients_by_site <- integer(0)
  }
  
  ## --- Capping at target sample size ---
  # Ensures that recruitment for the trial stops when the target is reached
  # It detects an overshoot of patients allocated and corrects it
  patients_recruited <- sum(patients_by_site)
  if (state$patients + patients_recruited > state$nPatients) {
    excess <- (state$patients + patients_recruited) - state$nPatients
    # Reduce one-by-one from sites with >0 until on target
    while (excess > 0L && sum(patients_by_site) > 0L) {
      idx <- which(patients_by_site > 0L)[1]
      patients_by_site[idx] <- patients_by_site[idx] - 1L
      excess <- excess - 1L
    }
    patients_recruited <- sum(patients_by_site)
  }
  
  ## --- Log site-level accrual for this calendar month ---
  if (n_sites > 0L) {
    site_df <- data.frame(
      Month    = rep(calendar_month, n_sites),
      Site_ID  = state$sites_vec,
      Patients = patients_by_site
    )
    state$site_accrual <- rbind(state$site_accrual, site_df)
  }
  
  ## --- Update monthly totals and flags ---
  state$patients             <- state$patients + patients_recruited
  state$patients_per_month   <- c(state$patients_per_month, patients_recruited)
  state$sites_per_month      <- c(state$sites_per_month, n_sites)
  state$month_numbers        <- c(state$month_numbers, calendar_month)
  
  if (state$patients >= state$nPatients) state$done <- TRUE
  state
}

# Step 3 - Return relevant data to the overall dataframe
finalize_recruitment_state <- function(state) {
  recruitment_details <- data.frame(
    Month = state$month_numbers,
    Active_Sites = state$sites_per_month,
    Patients = state$patients_per_month
  )
  start_m <- if (length(state$month_numbers)) min(state$month_numbers) else NA_integer_
  end_m   <- if (length(state$month_numbers)) max(state$month_numbers) else NA_integer_
  dur_m   <- if (length(state$month_numbers)) (end_m - start_m + 1L) else 0L
  
  list(
    recruitment_details = recruitment_details,
    total_sites    = if (length(state$sites_per_month)) tail(state$sites_per_month, 1) else 0L,
    total_patients = state$patients,
    start_month    = start_m,
    end_month      = end_m,
    total_months   = dur_m,
    site_accrual   = state$site_accrual
  )
}

# Overall simulation function for 4 staggered stand-alone trials that recruits on a shared
# global calendar, enforces a per-site concurrency cap on the number of trials and competes for 
# the same global site pool. It splits each sites accrual capacity across the trials that are
# concurrently running. Returns per-simulation, per-trial summaries and detailed monthly site logs

simulate_4_staggered_competing <- function(
    nSims,
    nPatients,
    mean_sites_per_month,
    max_sites,
    mean_patients_per_site_per_year,
    starts_months = starts_months,
    nTrials = nTrials,
    site_concurrency_cap = site_concurrency_cap,  # max concurrent trials per site
    global_nsites = max_sites                     # size of the global site universe
) {
  starts_months <- as.integer(starts_months)
  sims_out <- vector("list", nSims)
 
# Create one state per trial - each state logs patients, sites etc   
  for (s in seq_len(nSims)) {
    # init one state per trial
    states <- lapply(seq_len(nTrials), function(i) trct_init(
      nPatients, mean_sites_per_month, max_sites, mean_patients_per_site_per_year
    ))
    
    # GLOBAL SITE POOL for this simulation
    # active_count = how many trials this site is currently on
    global_pool <- data.frame(
      Site_ID = seq_len(global_nsites),
      active_count = integer(global_nsites),                 # how many trials the site is on this month
      base_rate_year = rep(mean_patients_per_site_per_year, global_nsites),  # e.g., 2
      stringsAsFactors = FALSE
    )
    # Tracking log for each simulation, month, how sites and accrual are distributed
    # Starts at null and is appended each month as trials activate or update weighting
    alloc_log <- data.frame(Sim=integer(), Month=integer(),
                            Trial=character(), K=integer(),
                            Sites_Added=integer(), Mean_Site_Weight=numeric(),
                            stringsAsFactors = FALSE)
    
    # tracks monthly progressing of which trials are active, have finished and when to stop 
    # the loop 
    was_done <- rep(FALSE, nTrials)
    
    month <- 0L
    while (!all(vapply(states, function(z) z$done, logical(1)))) {
      active <- which(month >= starts_months &
                        !vapply(states, `[[`, logical(1), "done"))
      K <- length(active)
      
      if (K > 0L) {
        
        ## ---- Site activation under concurrency cap (per active trial) ----
        # Randomly draws new sites and only selects sites not reached the cap.
        # Assigns up to the max sites per trial and globally
        # Logs number of sites added each month
        
        for (i in active) {
          desired <- rpois(1, states[[i]]$mean_sites_per_month)
          current_sites <- length(states[[i]]$sites_vec)
          remaining_trial <- max(0L, states[[i]]$max_sites - current_sites)
          
          n_add <- 0L
          if (remaining_trial > 0L && desired > 0L) {
            # available global sites: not already in this trial, and below cap
            avail <- global_pool$Site_ID[global_pool$active_count < site_concurrency_cap &
                                           !(global_pool$Site_ID %in% states[[i]]$sites_vec)]
            
            n_add <- min(desired, remaining_trial, length(avail))
            if (n_add > 0L) {
              # choose first n_add (or sample(avail, n_add))
              new_ids <- head(avail, n_add)
              # assign to trial
              states[[i]]$sites_vec <- c(states[[i]]$sites_vec, new_ids)
              states[[i]]$sites <- length(states[[i]]$sites_vec)
              # bump concurrency counters
              global_pool$active_count[match(new_ids, global_pool$Site_ID)] <-
                global_pool$active_count[match(new_ids, global_pool$Site_ID)] + 1L
            }
            }
          # compute current mean site weight for logging (after activation)
          if (length(states[[i]]$sites_vec)) {
            ac <- global_pool$active_count[match(states[[i]]$sites_vec, global_pool$Site_ID)]
            mean_w <- mean(1 / pmax(1L, ac))
          } else {
            mean_w <- NA_real_
          }
          
          # log activation for this month/trial
          alloc_log <- rbind(
            alloc_log,
            data.frame(Sim = s, Month = month, Trial = paste0("Trial_", i),
                       K = K, Sites_Added = n_add, Mean_Site_Weight = mean_w)
          )
        }
        ## ---- Step active trials with site-specific weights ----
        # Builds site weights for currently assigned sites based on number of trials at
        # each site and randomly draws patient accrual per site and caps at target sample size
        
        for (i in active) {
          site_ids <- states[[i]]$sites_vec
          if (length(site_ids) > 0L) {
            ac <- global_pool$active_count[match(site_ids, global_pool$Site_ID)]
            site_w <- 1 / pmax(1L, ac)    # with cap=2 → weights are 1.0 or 0.5
          } else {
            site_w <- NULL
          }
          
          states[[i]] <- trct_step(
            states[[i]],
            site_weights = site_w,        # << key change
            current_global_month = month
          )
        }
          # ---- Release sites from any trials that JUST finished this month ----
          # Detects trials newly finished each month and frees them for other trials
          now_done <- vapply(states, `[[`, logical(1), "done")
          newly_finished <- which(now_done & !was_done)
          if (length(newly_finished)) {
            for (i in newly_finished) {
              if (length(states[[i]]$sites_vec)) {
                idx <- match(states[[i]]$sites_vec, global_pool$Site_ID)
                global_pool$active_count[idx] <- pmax(0L, global_pool$active_count[idx] - 1L)
              }
              # NOTE: we DO NOT clear states[[i]]$sites_vec so logs remain consistent;
              # the pool is freed for future assignment.
            }
          }
          was_done <- now_done
        }
      
      month <- month + 1L
    }
# Build a per-trial summary 
    trial_ids <- paste0("Trial_", seq_len(nTrials))
    
    pt <- setNames(lapply(states, finalize_recruitment_state), trial_ids)
    
    sim_summary <- data.frame(
      Sim = s,
      Trial = trial_ids,
      Start_Month = vapply(pt, function(x) x$start_month, integer(1)),
      End_Month   = vapply(pt, function(x) x$end_month, integer(1)),
      Duration_Months = vapply(pt, function(x) x$total_months, integer(1)),
      Total_Sites_At_End = vapply(states, function(z) as.integer(if (length(z$sites_per_month)) tail(z$sites_per_month, 1) else 0L), integer(1)),
      Total_Patients = vapply(states, function(z) as.integer(z$patients), integer(1)),
      stringsAsFactors = FALSE
    )
    
    program_start <- min(sim_summary$Start_Month, na.rm = TRUE)
    program_end   <- max(sim_summary$End_Month,   na.rm = TRUE)
    program_duration <- program_end - program_start + 1L
    
    rec_logs <- do.call(rbind, Map(function(id, st) {
      if (length(st$month_numbers)) {
        data.frame(Sim = s, Trial = id, Month = st$month_numbers,
                   Active_Sites = st$sites_per_month, Patients = st$patients_per_month)
      } else data.frame(Sim = integer(0), Trial = character(0),
                        Month = integer(0), Active_Sites = integer(0), Patients = integer(0))
    }, trial_ids, states))
    
    sims_out[[s]] <- list(
      simulation_results = sim_summary,
      recruitment_details = rec_logs,
      per_trial = pt,
      alloc_log = alloc_log,
      program_summary = data.frame(
        Sim = s,
        Program_Start = program_start,
        Program_End = program_end,
        Program_Duration_Months = program_duration
    )
    )
  }
  
  list(
    simulation_results  = do.call(rbind, lapply(sims_out, `[[`, "simulation_results")),
    recruitment_details = do.call(rbind, lapply(sims_out, `[[`, "recruitment_details")),
    sims = sims_out
  )
}

##############################################################


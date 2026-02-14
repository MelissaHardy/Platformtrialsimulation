################################################################################################################
############################ R script to simulate a Parallel Platform Trial with 5 arms#########################
################################################################################################################
# This code will simulate recruitment of a 5-arm non-inferiority randomised controlled platform trial using 1:1 
# allocation of participants to treatment 0 ("best available therapy") versus treatment 1 versus treatment 2 
# versus treatment 3 and treatment 4

# This function will estimate the total duration of the trial based on the number of sites activated per month and
# the number of patients recruited per year

simulate_recruitment_ppt <- function(nPatients, mean_sites_per_month, max_sites, mean_patients_per_site_per_year) {
  sites <- 0
  patients <- 0
  month <- 0
  
  # Vectors to track recruitment over time
  patients_per_month <- c()
  sites_per_month <- c()
  month_numbers <- c()
  
  # Initialize site-level accrual data frame
  site_accrual <- data.frame()
  
  while (patients < nPatients) {
    month <- month + 1
    
    # Activate new sites (ensuring we don't exceed max_sites)
    if (sites < max_sites) {
      sites_activated <- rpois(1, mean_sites_per_month)
      if (sites + sites_activated > max_sites) {
        sites_activated <- max_sites - sites
      }
      sites <- sites + sites_activated
    }
    
    # Simulate patients by site
    patients_by_site <- rpois(sites, mean_patients_per_site_per_year / 12)
    
    # Cap if total exceeds required
    patients_recruited <- sum(patients_by_site)
    if (patients + patients_recruited > nPatients) {
      excess <- (patients + patients_recruited) - nPatients
      # Reduce patient count proportionally from sites with recruitment
      while (excess > 0 && sum(patients_by_site) > 0) {
        idx <- which(patients_by_site > 0)[1]
        patients_by_site[idx] <- patients_by_site[idx] - 1
        excess <- excess - 1
      }
      patients_recruited <- sum(patients_by_site)
    }
    
    # Track site-level accrual
    site_df <- data.frame(
      Month = rep(month, sites),
      Site_ID = seq_len(sites),
      Patients = patients_by_site
    )
    site_accrual <- rbind(site_accrual, site_df)
    
    patients <- patients + patients_recruited
    patients_per_month <- c(patients_per_month, patients_recruited)
    sites_per_month <- c(sites_per_month, sites)
    month_numbers <- c(month_numbers, month)
  }
  
  recruitment_details <- data.frame(
    Month = month_numbers,
    Active_Sites = sites_per_month,
    Patients = patients_per_month
  )
  
  list(
    recruitment_details = recruitment_details,
    total_sites = sites,
    total_patients = patients,
    total_months = month,
    site_accrual = site_accrual
  )
}

################################################################################################################
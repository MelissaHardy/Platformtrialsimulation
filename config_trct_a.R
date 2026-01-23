##################################################################################
#### Trial Design Parameters for traditional stand-alone trials for Scenario A####
##################################################################################

# The type 1 error is 0.025
# The power is 80%
# Sample size was predetermined

nPatients <- 120 # number of patients in each RCT
death0 <- 0.41 # Mortality rate for the control arm
death1 <- 0.36 # Mortality rate for the treatment arms
rue_rd <- death1 - death0 # True risk difference
ni_margin <- 0.20 # non-inferiority margin for Scenario A
death_alt <- death0 + ni_margin  # Experimental treatment is inferior (worse outcome rate) for the type 1 error rate
mean_sites_per_month <- 10 # mean number of sites activated per month
max_sites <- 200 # maximum number of sites for the trial
mean_patients_per_site_per_year <- 2 # mean number of patients recruited per site per year
nSims <- 10000 # number of trials to simulate
starts_months <- as.integer(c(0, 6, 12, 18)) # Required for staggered trials
nTrials <- 4
site_concurrency_cap <- 2 # Number of trials each site could concurrently participate in

# set seed
set.seed(1234)

##################################################################################
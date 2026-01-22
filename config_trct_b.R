
#### Trial Design Parameters ####
# The type 1 error is 0.025
# The power is 80%
# Total sample size was predetermined
# The total duration required the mean number of patients recruited per site per year

nPatients <- 186 # number of patients in each RCT
death0 <- 0.41 # mortality rate for patients receiving 'treatment 0' BAT
death1 <- 0.36 # mortality rate for patients receiving 'treatment 1' new investigational treatment
true_rd <- death1 - death0 # True risk difference
ni_margin <- 0.15 # non-inferiority margin
death_alt <- death0 + ni_margin  # Experimental treatment is inferior (worse outcome rate)
mean_sites_per_month <- 10 # mean number of sites activated per month
max_sites <- 200 # maximum number of sites for the trial
mean_patients_per_site_per_year <- 2 # mean number of patients recruited per site per year
nSims <- 10000 # number of trials to simulate

starts_months <- as.integer(c(0, 6, 12, 18)) #Required for the recruitment simulation
nTrials <- 4
site_concurrency_cap <- 2

# set seed
set.seed(1234)

##########################################

#### Trial Design Parameters ####
# The type 1 error is 0.025
# The power is 80%
# Total sample size was predetermined
# The total duration required the mean number of patients recruited per site per year

n_arms <- 5 # Number of treatment arms including BAT

# Total number of patients in each arm 
patients_per_experimental_arm <- 71
death0 <- 0.41
death1 <- 0.36
death_probs = c(0.41, 0.36, 0.36, 0.36, 0.36)  # mortality rate for the 5 arms
true_rd <- death1 - death0 # True risk difference

# non-inferiority margin - need to change from 0.20 to 0.15 for second study
ni_margin <- 0.15  
death_alt <- death0 + ni_margin  # Experimental treatment is inferior (worse outcome rate)
mean_sites_per_month <- 10 # mean number of sites activated per month
max_sites <- 200 # maximimum number of sites for the trial
mean_patients_per_site_per_year <- 2 # mean number of patients recruited per site per year
stagger_months <- 6 # Each new trial is staggered by 6 months
arm_open_months <- c(0, 6, 12, 18)
arm_open_month <- c(0, 0, stagger_months * 1:(n_arms - 2))
expected_control_arm <- 141


# Inferior death probabilities vector: control remains at 0.41, others worse
death_probs_inferior <- c(death0, rep(death_alt, n_arms - 1))

nSims <- 10000

# set seed
set.seed(1234)

##########################################
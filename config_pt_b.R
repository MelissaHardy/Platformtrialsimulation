#####################################################################
#### Trial Design Parameters for the platform trials (Scenario B)####
#####################################################################

# The type 1 error is 0.025
# The power is 80%
# Sample size was predetermined

# Total number of patients the platform trial
nPatients <- 465 

n_arms <- 5 # Number of treatment arms including control arm

patients_per_experimental_arm <- 93 # Total number of patients in each arm in Scenario B
death0 <- 0.41 # Mortality rate for the control arm
death1 <- 0.36 # Mortality rate for the treatment arms
death_probs = c(0.41, 0.36, 0.36, 0.36, 0.36)  # Mortality rate for the control and 4 treatment arms
true_rd <- death1 - death0 # True risk difference

ni_margin <- 0.15  # non-inferiority margin fro scenario B
death_alt <- death0 + ni_margin  # Experimental treatment is inferior (worse outcome rate) for the type 1 error rate
mean_sites_per_month <- 10 # mean number of sites activated per month
max_sites <- 200 # maximum number of sites for the trial
mean_patients_per_site_per_year <- 2 # mean number of patients recruited per site per year
stagger_months <- 6 # Each new trial is staggered by 6 months
arm_open_months <- c(0, 6, 12, 18)
arm_open_month <- c(0, 0, stagger_months * 1:(n_arms - 2))
expected_control_arm <- 93


# Inferior death probabilities vector: control remains at 0.41, others worse
death_probs_inferior <- c(death0, rep(death_alt, n_arms - 1))

nSims <- 10000

# set seed
set.seed(1234)

#####################################################################
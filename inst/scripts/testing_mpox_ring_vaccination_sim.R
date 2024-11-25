source("R/mpox_ring_vaccination_offspring_function.R")
seed <- 120
population <- 10^6
initial_immune <- 0
check_final_size <- 10000
load("data/south_kivu_final.rda")
synthetic_household_df <- data.table::as.data.table(south_kivu_final)
seeding_cases <- 5
t0 <- 0
infection_to_onset <- function(n) { rep(5, n) }
prob_quarantine <- 0.4
onset_to_quarantine <- function(n) { rep(2, n) }
vaccine_logistical_delay <- 2
vaccine_protection_delay <- 2
mn_offspring_sexual_SW <- 5
disp_offspring_sexual_SW <- 100
mn_offspring_sexual_PBS <- 5
disp_offspring_sexual_PBS <- 10
mn_offspring_sexual_genPop <- 5
disp_offspring_sexual_genPop <- 10
sexual_transmission_occupation_matrix <- matrix(data = c(0.00, 0.10, 0.90,
                                                         0.20, 0.10, 0.70,
                                                         0.00, 0.90, 0.10), nrow = 3, ncol = 3, byrow = TRUE)
mn_offspring_hh <- 7
disp_offspring_hh <- 100
mn_offspring_community <- 5
disp_offspring_community <- 100
community_transmission_age_matrix <- matrix(data = c(0.34, 0.33, 0.33,
                                                     0.34, 0.33, 0.33,
                                                     0.34, 0.33, 0.33), nrow = 3, ncol = 3, byrow = TRUE)
vaccine_efficacy_transmission <- 0.5
generation_time <- function(n) { rep(10, n)}
vaccine_coverage <- 0.5
vaccine_efficacy_infection <- 0.8
prob_quarantine <- 0.8
onset_to_quarantine <- function(n) { rep(0.1, n) }
quarantine_efficacy_household <- 0.77
quarantine_efficacy_else <- 0.6
vaccine_start <- 0
prop_asymptomatic <- 0.1

x <- basic_ring_vaccination_sim(mn_offspring_sexual_SW = mn_offspring_sexual_SW,
                                disp_offspring_sexual_SW = disp_offspring_sexual_SW,
                                mn_offspring_sexual_PBS = mn_offspring_sexual_PBS,
                                disp_offspring_sexual_PBS = disp_offspring_sexual_PBS,
                                mn_offspring_sexual_genPop = mn_offspring_sexual_genPop,
                                disp_offspring_sexual_genPop = disp_offspring_sexual_genPop,
                                sexual_transmission_occupation_matrix = sexual_transmission_occupation_matrix,
                                mn_offspring_hh = mn_offspring_hh,
                                disp_offspring_hh = disp_offspring_hh,
                                mn_offspring_community = mn_offspring_community,
                                disp_offspring_community = disp_offspring_community,
                                community_transmission_age_matrix = community_transmission_age_matrix,
                                prop_asymptomatic = prop_asymptomatic,
                                generation_time = generation_time,
                                infection_to_onset = infection_to_onset,
                                vaccine_start = vaccine_start,
                                vaccine_coverage = vaccine_coverage,
                                vaccine_efficacy_infection = vaccine_efficacy_infection,
                                vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                vaccine_logistical_delay = vaccine_logistical_delay,
                                vaccine_protection_delay = vaccine_protection_delay,
                                prob_quarantine = prob_quarantine,
                                onset_to_quarantine = onset_to_quarantine,
                                quarantine_efficacy_household = quarantine_efficacy_household,
                                quarantine_efficacy_else = quarantine_efficacy_else,
                                synthetic_household_df = synthetic_household_df,
                                t0 = t0,
                                population = population,
                                check_final_size = check_final_size,
                                initial_immune = initial_immune,
                                seeding_cases = seeding_cases,
                                seed = seed)

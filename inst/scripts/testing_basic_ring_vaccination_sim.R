## Testing the ring vaccination function
mn_offspring <- 6         # mean of the offspring distribution
disp_offspring <- 1.1     # overdispersion of the offspring distribution (1 = Poisson, but has to be > 1 for Negative Binomial, which we use)
t0 <- 0                   # calendar time to start at
population <- 10^8        # population size (this is largely irrelevant, only used in the susceptible depletion term for modifying the offspring distribution)
check_final_size <- 10000 # maximum number of infections to simulate
initial_immune <- 0       # number of individuals who are initially immune
seeding_cases <- 10       # number of initially infected individuals
seed <- 10                # stochastic seed for reproducibility

## Sourcing the function (in functions/basic_ring_vaccination_sim.R)
source("R/basic_ring_vaccination_sim.R")

## Key features of the model are:
#### Ring vaccination (of both contacts AND contacts of contacts) - initiated following symptom onset in infections
#### Quarantine - initiated following symptom onset in infections

## Example 1: No ring-vaccination, quarantine and symptom onset happens very early relative to generation time, so we expect it to be effective
### Should expect R to be reduced to (prob_quarantine * mn_offspring) + (prob_quarantine * (1 - quarantine_efficacy) * mn_offspring)
example1 <- basic_ring_vaccination_sim(mn_offspring = mn_offspring,
                                       disp_offspring = disp_offspring,
                                       t0 = t0,
                                       population = population,
                                       check_final_size = check_final_size,
                                       initial_immune = initial_immune,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = 0,
                                       generation_time = function(n) { rgamma(n, shape = 50, rate = 2) },
                                       infection_to_onset = function(n) { rgamma(n, shape = 2, rate = 2) },
                                       vaccine_start = 100000,
                                       vaccine_coverage = 0,
                                       vaccine_efficacy_infection = 0,
                                       vaccine_efficacy_transmission = 0,
                                       vaccine_logistical_delay = 10000,
                                       vaccine_protection_delay = 10000,
                                       prob_quarantine = 0.5,
                                       onset_to_quarantine = function(n) { rgamma(n, shape = 5, rate = 2) },
                                       quarantine_efficacy = 0.8,
                                       seed = seed)
mean(example1$n_offspring[!is.na(example1$time_infection_absolute)], na.rm = TRUE)
mean(example1$n_offspring_new[!is.na(example1$time_infection_absolute)], na.rm = TRUE)
mean(example1$n_offspring_quarantine[!is.na(example1$time_infection_absolute)], na.rm = TRUE)
mean(example1$n_offspring_post_pruning[!is.na(example1$time_infection_absolute)], na.rm = TRUE)
mean(example1$n_offspring_post_pruning_2nd_chance[!is.na(example1$time_infection_absolute)], na.rm = TRUE)

## Example 2: No ring-vaccination, quarantine but symptom onset happens very late relative to generation time, so we expect it to be ineffective
### Should expect R to be unchanged i.e. approximately the value of mn_offspring
example2 <- basic_ring_vaccination_sim(mn_offspring = mn_offspring,
                                       disp_offspring = disp_offspring,
                                       t0 = t0,
                                       population = population,
                                       check_final_size = check_final_size,
                                       initial_immune = initial_immune,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = 0,
                                       generation_time = function(n) { rgamma(n, shape = 5, rate = 2) },
                                       infection_to_onset = function(n) { rgamma(n, shape = 50, rate = 2) },
                                       vaccine_start = 100000,
                                       vaccine_coverage = 0,
                                       vaccine_efficacy_infection = 0,
                                       vaccine_efficacy_transmission = 0,
                                       vaccine_logistical_delay = 10000,
                                       vaccine_protection_delay = 10000,
                                       prob_quarantine = 0.5,
                                       onset_to_quarantine = function(n) { rgamma(n, shape = 5, rate = 2) },
                                       quarantine_efficacy = 0.8,
                                       seed = seed)
mean(example2$n_offspring[!is.na(example2$time_infection_absolute)], na.rm = TRUE)
mean(example2$n_offspring_new[!is.na(example2$time_infection_absolute)], na.rm = TRUE)
mean(example2$n_offspring_quarantine[!is.na(example2$time_infection_absolute)], na.rm = TRUE)
mean(example2$n_offspring_post_pruning[!is.na(example2$time_infection_absolute)], na.rm = TRUE)
mean(example2$n_offspring_post_pruning_2nd_chance[!is.na(example2$time_infection_absolute)], na.rm = TRUE)

## Example 3: No ring-vaccination, quarantine and all infections are asymptomatic so we expect it to be ineffective
### Should expect R to be unchanged i.e. approximately the value of mn_offspring
example3 <- basic_ring_vaccination_sim(mn_offspring = mn_offspring,
                                       disp_offspring = disp_offspring,
                                       t0 = t0,
                                       population = population,
                                       check_final_size = check_final_size,
                                       initial_immune = initial_immune,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = 1,
                                       generation_time = function(n) { rgamma(n, shape = 50, rate = 2) },
                                       infection_to_onset = function(n) { rgamma(n, shape = 2, rate = 2) },
                                       vaccine_start = 100000,
                                       vaccine_coverage = 0,
                                       vaccine_efficacy_infection = 0,
                                       vaccine_efficacy_transmission = 0,
                                       vaccine_logistical_delay = 10000,
                                       vaccine_protection_delay = 10000,
                                       prob_quarantine = 0.5,
                                       onset_to_quarantine = function(n) { rgamma(n, shape = 5, rate = 2) },
                                       quarantine_efficacy = 0.8,
                                       seed = seed)
mean(example3$n_offspring[!is.na(example3$time_infection_absolute)], na.rm = TRUE)
mean(example3$n_offspring_new[!is.na(example3$time_infection_absolute)], na.rm = TRUE)
mean(example3$n_offspring_quarantine[!is.na(example3$time_infection_absolute)], na.rm = TRUE)
mean(example3$n_offspring_post_pruning[!is.na(example3$time_infection_absolute)], na.rm = TRUE)
mean(example3$n_offspring_post_pruning_2nd_chance[!is.na(example3$time_infection_absolute)], na.rm = TRUE)

## Example 4: No quarantine, ring-vaccination and symptom onset happens very early relative to generation time. Perfect coverage so individuals effectively get one chance at protection.
##            No infections asymptomatic.
### Should expect R to be ~ (mn_offspring * 1 * 0.5)
example4 <- basic_ring_vaccination_sim(mn_offspring = mn_offspring,
                                       disp_offspring = disp_offspring,
                                       t0 = t0,
                                       population = population,
                                       check_final_size = check_final_size,
                                       initial_immune = initial_immune,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = 0,
                                       generation_time = function(n) { rgamma(n, shape = 50, rate = 2) },
                                       infection_to_onset = function(n) { rgamma(n, shape = 2, rate = 2) },
                                       vaccine_start = 0,
                                       vaccine_coverage = 1,
                                       vaccine_efficacy_infection = 0.5,
                                       vaccine_efficacy_transmission = 0,
                                       vaccine_logistical_delay = 1,
                                       vaccine_protection_delay = 1,
                                       prob_quarantine = 0,
                                       onset_to_quarantine = function(n) { rgamma(n, shape = 5, rate = 2) },
                                       quarantine_efficacy = 0,
                                       seed = seed)
mean(example4$n_offspring[!is.na(example4$time_infection_absolute)], na.rm = TRUE)
mean(example4$n_offspring_new[!is.na(example4$time_infection_absolute)], na.rm = TRUE)
mean(example4$n_offspring_quarantine[!is.na(example4$time_infection_absolute)], na.rm = TRUE)
mean(example4$n_offspring_post_pruning[!is.na(example4$time_infection_absolute)], na.rm = TRUE)
mean(example4$n_offspring_post_pruning_2nd_chance[!is.na(example4$time_infection_absolute)], na.rm = TRUE)


## Example 5: No quarantine, ring-vaccination and symptom onset happens very early relative to generation time.
##            Perfect coverage so individuals effectively get one chance at protection <- I wrote this originally
##            and whilst it's technically true, some infections are asymptomatic and so individuals might get their "first chance"
##            at the "2nd chance of vaccination stage" (i.e. if index was asymptomatic hence no ring vacc, but secondary was symptomatic)
example5 <- basic_ring_vaccination_sim(mn_offspring = mn_offspring,
                                       disp_offspring = disp_offspring,
                                       t0 = t0,
                                       population = population,
                                       check_final_size = check_final_size,
                                       initial_immune = initial_immune,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = 0.5,
                                       generation_time = function(n) { rgamma(n, shape = 50, rate = 2) },
                                       infection_to_onset = function(n) { rgamma(n, shape = 1, rate = 2) },
                                       vaccine_start = 0,
                                       vaccine_coverage = 1,
                                       vaccine_efficacy_infection = 0.5,
                                       vaccine_efficacy_transmission = 0,
                                       vaccine_logistical_delay = 0.1,
                                       vaccine_protection_delay = 0.1,
                                       prob_quarantine = 0,
                                       onset_to_quarantine = function(n) { rgamma(n, shape = 5, rate = 2) },
                                       quarantine_efficacy = 0,
                                       seed = 100)
mean(example5$n_offspring[!is.na(example5$time_infection_absolute)], na.rm = TRUE)
mean(example5$n_offspring_new[!is.na(example5$time_infection_absolute)], na.rm = TRUE)
mean(example5$n_offspring_quarantine[!is.na(example5$time_infection_absolute)], na.rm = TRUE)
x <- c()
for (i in 1:nrow(example5)) {
  if (example5$secondary_offspring_generated[i] & example5$tertiary_offspring_generated[i]) {
    temp_total_offspring <- example5$n_offspring_post_pruning[i]
    indices <- which(example5$parent == example5$id[i])
    temp_total_offspring <- temp_total_offspring + sum(example5$n_offspring_post_pruning[indices])
    x <- c(x, temp_total_offspring)
  }
}
mean(sqrt(x)) # number of offspring per infection (mean)
mean(example5$n_offspring_post_pruning_2nd_chance[!is.na(example5$time_infection_absolute)], na.rm = TRUE)
mean(example5$n_offspring_post_pruning[!is.na(example5$time_infection_absolute)], na.rm = TRUE) # unclear to me currently why this is ~3.2 and lines 162-170 comes to 3 (which is what I expected)

# I need to understand why doing the below isn't reasonable in this case (and why in others it is)
mean(example5$n_offspring_post_pruning[!is.na(example5$time_infection_absolute) &
                                         example5$secondary_offspring_generated &
                                         example5$tertiary_offspring_generated], na.rm = TRUE) # unclear to me currently why this is ~3.2 and lines 162-170 comes to 3 (which is what I expected)

## Example 6: No quarantine, ring-vaccination and symptom onset happens very early relative to generation time.
##            Imperfect coverage so individuals effectively get two chances at protection.
example6 <- basic_ring_vaccination_sim(mn_offspring = mn_offspring,
                                       disp_offspring = disp_offspring,
                                       t0 = t0,
                                       population = population,
                                       check_final_size = check_final_size,
                                       initial_immune = initial_immune,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = 0.5,
                                       generation_time = function(n) { rgamma(n, shape = 50, rate = 2) },
                                       infection_to_onset = function(n) { rgamma(n, shape = 1, rate = 2) },
                                       vaccine_start = 0,
                                       vaccine_coverage = 0.5,
                                       vaccine_efficacy_infection = 0.5,
                                       vaccine_efficacy_transmission = 0,
                                       vaccine_logistical_delay = 0.1,
                                       vaccine_protection_delay = 0.1,
                                       prob_quarantine = 0,
                                       onset_to_quarantine = function(n) { rgamma(n, shape = 5, rate = 2) },
                                       quarantine_efficacy = 0,
                                       seed = 100)
mean(example6$n_offspring[!is.na(example6$time_infection_absolute)], na.rm = TRUE)
mean(example6$n_offspring_new[!is.na(example6$time_infection_absolute)], na.rm = TRUE)
mean(example6$n_offspring_quarantine[!is.na(example6$time_infection_absolute)], na.rm = TRUE)
x <- c()
for (i in 1:nrow(example6)) {
  if (example6$secondary_offspring_generated[i] & example6$tertiary_offspring_generated[i]) {
    temp_total_offspring <- example6$n_offspring_post_pruning[i]
    indices <- which(example6$parent == example6$id[i])
    temp_total_offspring <- temp_total_offspring + sum(example6$n_offspring_post_pruning[indices])
    x <- c(x, temp_total_offspring)
  }
}
mean(sqrt(x)) # number of offspring per infection (mean)
mean(example6$n_offspring_post_pruning_2nd_chance[!is.na(example6$time_infection_absolute)], na.rm = TRUE)
mean(example6$n_offspring_post_pruning[!is.na(example6$time_infection_absolute)], na.rm = TRUE) # unclear to me currently why this is ~3.2 and lines 162-170 comes to 3 (which is what I expected)

## Example 7: No quarantine, ring-vaccination and symptom onset happens very early relative to generation time.
##            Imperfect coverage so individuals effectively get two chances at protection.
##            Vaccine provides reduction in transmission in breakthrough infections.
example7 <- basic_ring_vaccination_sim(mn_offspring = 3,
                                       disp_offspring = disp_offspring,
                                       t0 = t0,
                                       population = population,
                                       check_final_size = 25000,
                                       initial_immune = initial_immune,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = 0,
                                       generation_time = function(n) { rgamma(n, shape = 50, rate = 2) },
                                       infection_to_onset = function(n) { rgamma(n, shape = 1, rate = 2) },
                                       vaccine_start = 0,
                                       vaccine_coverage = 1,
                                       vaccine_efficacy_infection = 0,
                                       vaccine_efficacy_transmission = 0.5,
                                       vaccine_logistical_delay = 0.1,
                                       vaccine_protection_delay = 0.1,
                                       prob_quarantine = 0,
                                       onset_to_quarantine = function(n) { rgamma(n, shape = 5, rate = 2) },
                                       quarantine_efficacy = 0,
                                       seed = 100)
mean(example7$n_offspring[!is.na(example7$time_infection_absolute)], na.rm = TRUE)
mean(example7$n_offspring_new[!is.na(example7$time_infection_absolute)], na.rm = TRUE)
mean(example7$n_offspring_quarantine[!is.na(example7$time_infection_absolute)], na.rm = TRUE)
mean(example7$n_offspring_post_pruning[!is.na(example7$time_infection_absolute)], na.rm = TRUE)
mean(example7$n_offspring_post_pruning_2nd_chance[!is.na(example7$time_infection_absolute)], na.rm = TRUE) # this is because when vaccine_start is 1, none of the secondary and tertiary cases of seeding cases are vaccinated in the first round
                                                                                                           # (and there's a lot of them) - they're vaccinated at the second attempt (I think).

x <- example7[!is.na(example7$time_infection_absolute) &
                !is.na(example7$n_offspring_post_pruning_2nd_chance), ]
mean(x$n_offspring_post_pruning)
mean(x$n_offspring_post_pruning_2nd_chance)

x <- c()
for (i in 1:nrow(example7)) {
  if (example7$secondary_offspring_generated[i] & example7$tertiary_offspring_generated[i]) {
    temp_total_offspring <- example7$n_offspring_post_pruning[i]
    indices <- which(example7$parent == example7$id[i])
    temp_total_offspring <- temp_total_offspring + sum(example7$n_offspring_post_pruning[indices])
    x <- c(x, temp_total_offspring)
  }
}
mean(sqrt(x)) # number of offspring per infection (mean)

mean(example7$n_offspring_post_pruning[!is.na(example7$time_infection_absolute) & example7$generation != 1], na.rm = TRUE)
mean(example7$n_offspring_post_pruning_2nd_chance[!is.na(example7$time_infection_absolute) & example7$generation != 1], na.rm = TRUE)

mean(example7$n_offspring_post_pruning[!is.na(example7$time_infection_absolute) &
                                         example7$secondary_offspring_generated &
                                         example7$tertiary_offspring_generated], na.rm = TRUE)
mean(example7$n_offspring_post_pruning_2nd_chance[!is.na(example7$time_infection_absolute) &
                                         example7$secondary_offspring_generated &
                                         example7$tertiary_offspring_generated], na.rm = TRUE) # unclear to me currently why this is ~3.2 and lines 162-170 comes to 3 (which is what I expected)

hist(example7$n_offspring - example7$n_offspring_new, breaks = 20)
hist(example7$n_offspring_new - example7$n_offspring_quarantine, breaks = 20)
hist(example7$n_offspring_quarantine - example7$n_offspring_post_pruning, breaks = 20)
hist(example7$n_offspring_post_pruning - example7$n_offspring_post_pruning_2nd_chance, breaks = 20)














## Scrap Code

# generation_time <- function(n) { rgamma(n, shape = 25, rate = 2) }
# prop_asymptomatic <- 0
# infection_to_onset <- function(n) { rgamma(n, shape = 2, rate = 2) }
# vaccine_start <- 0
# vaccine_coverage <- 0.5
# vaccine_efficacy_infection <- 1
# vaccine_efficacy_transmission <- 0
# vaccine_logistical_delay <- 0.01
# vaccine_protection_delay <- 0.01
# prob_quarantine <- 0
# onset_to_quarantine <- function(n) { rgamma(n, shape = 5, rate = 2) }
# quarantine_efficacy <- 0
#
# length(unique(tdf$id))
#
# x <- c()
# for (i in 1:nrow(tdf)) {
#   if (tdf$secondary_offspring_generated[i] & tdf$tertiary_offspring_generated[i]) {
#     temp_total_offspring <- tdf$n_offspring_post_pruning[i]
#     indices <- which(tdf$parent == tdf$id[i])
#     temp_total_offspring <- temp_total_offspring + sum(tdf$n_offspring_post_pruning[indices])
#     x <- c(x, temp_total_offspring)
#   }
# }
#
# mean(sqrt(x))
# hist(sqrt(x), breaks = 20)
#
# hist(tdf$n_offspring[!is.na(tdf$time_infection_absolute)], breaks = 20)
# hist(tdf$n_offspring_new[!is.na(tdf$time_infection_absolute)], breaks = 20)
# hist(tdf$n_offspring_quarantine[!is.na(tdf$time_infection_absolute)], breaks = 20)
# hist(tdf$n_offspring_post_pruning[!is.na(tdf$time_infection_absolute)], breaks = 20)
# hist(tdf$n_offspring_post_pruning_2nd_chance[!is.na(tdf$time_infection_absolute)], breaks = 20)
#
# mean(tdf$n_offspring[!is.na(tdf$time_infection_absolute)], na.rm = TRUE)
# mean(tdf$n_offspring_new[!is.na(tdf$time_infection_absolute)], na.rm = TRUE)
# mean(tdf$n_offspring_quarantine[!is.na(tdf$time_infection_absolute)], na.rm = TRUE)
# mean(tdf$n_offspring_post_pruning[!is.na(tdf$time_infection_absolute)], na.rm = TRUE)
# mean(tdf$n_offspring_post_pruning_2nd_chance[!is.na(tdf$time_infection_absolute)], na.rm = TRUE)


# sum(example5$vaccinated[!is.na(example5$time_infection_absolute)])
# table(example5$vaccinated, useNA = "ifany")
# table(example5$vaccinated_2nd_chance, useNA = "ifany")
#
# table(example5$asymptomatic, useNA = "ifany")
#
# mean(example5$n_offspring_post_pruning[!is.na(example5$time_infection_absolute) & example5$asymptomatic == 0], na.rm = TRUE)
# mean(example5$n_offspring_post_pruning[!is.na(example5$time_infection_absolute) & example5$asymptomatic == 1], na.rm = TRUE)
#
# mean(example5$n_offspring_post_pruning_2nd_chance[!is.na(example5$time_infection_absolute) & example5$asymptomatic == 0], na.rm = TRUE)
# mean(example5$n_offspring_post_pruning_2nd_chance[!is.na(example5$time_infection_absolute) & example5$asymptomatic == 1], na.rm = TRUE)
#
# sum(!is.na(example5$time_infection_absolute) & example5$asymptomatic == 0)
# sum(!is.na(example5$time_infection_absolute) & example5$asymptomatic == 1)
#
# 0.5 * 6
#
# hist(example5$n_offspring[!is.na(example5$time_infection_absolute)], na.rm = TRUE)
# hist(example5$n_offspring_new[!is.na(example5$time_infection_absolute)], na.rm = TRUE)
# hist(example5$n_offspring_quarantine[!is.na(example5$time_infection_absolute)], na.rm = TRUE)
# hist(example5$n_offspring_post_pruning[!is.na(example5$time_infection_absolute)], na.rm = TRUE)
# hist(example5$n_offspring_post_pruning_2nd_chance[!is.na(example5$time_infection_absolute)], na.rm = TRUE)

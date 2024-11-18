#########################################################################################################################################################################################################
## Overview: The basic_ring_vaccination_sim function simulates the spread of an infectious disease within a population using a stochasting branching process,
##           incorporating the effects of ring vaccination and quarantine measures. The function begins by setting up the initial conditions and pre-allocating
##           a data frame to store details of the infections that are generated details. It uses a negative binomial distribution to model the number of secondary
##           infections each individual produces, adjusted for the susceptible population size. The simulation progresses iteratively, generating secondary and tertiary infections
##           for each index infection, and applying both ring vaccination and quarantine interventions. The simulation continues until there are no new infections
##           to process or the maximum number of infections (check_final_size) is reached. The function outputs a data frame containing detailed information on each simulated infection,
##           see below (lines 44-74) for the meaning of each data frame column.
#########################################################################################################################################################################################################

#########################################################################################################################################################################################################
## Notes (just jotting down thoughts/uncertainties and/or things to address as they come to mind)
##
##  Note #1: We use two rings for ring vaccination (i.e. contacts of the index case AND contacts of these contacts can get vaccinated) but we DON'T re-vaccinate people who have already
##  been vaccinated and for whom vaccination fails to protect them. This can produce some weird behaviour - for example, imperfect coverage with a good vaccine is better than perfect coverage
##  with a worse vaccine. Because in the latter, you vaccinate everyone first time around and they don't get another chance at being vaccinated if it fails. In the former, unvaccinated individuals
##  get another chance at being vaccinated the second time around. Don't think this is necessarily a problem, but good to be aware of.
##
#########################################################################################################################################################################################################

#########################################################################################################################################################################################################
## To Consider Building In:
##
## - Make prob_contact_traced a function of the number of active infections (or something that means it's time varying and related to the size of the epidemic)
## - Explicit tracking of the number of vaccine doses required - will require consideration of contacts that DON'T otherwise become infections (could be done post-hoc instead)
## - Limits to the number of doses and its impact on probability of someone being successfully vaccinated
## - Ability to specify whether 1 ring or 2 rings of ring-vaccination
## - Ability to specify whether vaccine coverage is due to being missed by health system (in which case 2nd attempt might work) or due to hesitancy/refusal (in which case 2nd attemtp won't work)
##
#########################################################################################################################################################################################################

#########################################################################################################################################################################################################
## To Definitely Build In:
##
## - Prob quarantine conditional on symptoms, rather than being independent of it
##
#########################################################################################################################################################################################################

#' @export
basic_ring_vaccination_sim <- function(mn_offspring,                       # mean of the offspring distribution
                                       disp_offspring,                     # overdispersion of the offspring distribution
                                       generation_time,                    # generation time distribution
                                       t0 = 0,                             # simulation starting time
                                       population,                         # total population size
                                       check_final_size,                   # maximum number of infections to simulate
                                       initial_immune = 0,                 # initial number of individuals who are immune
                                       seeding_cases,                      # number of seeding cases to start the epidemic with
                                       prop_asymptomatic,                  # probability an infection is asymptomatic
                                       infection_to_onset,                 # infection to symptom onset time distribution
                                       vaccine_start,                      # time at which the vaccine becomes available
                                       vaccine_coverage,                   # probability that each eligible individual gets vaccinated
                                       vaccine_efficacy_infection,         # vaccine efficacy against infection
                                       vaccine_efficacy_transmission,      # reduction in transmissibility of breakthrough infections in vaccinated individuals
                                       vaccine_logistical_delay,           # delay between symptom onset and vaccination of contacts (and contacts of contacts) occurring
                                       vaccine_protection_delay,           # delay between vaccination and protection developing
                                       prob_quarantine,                    # probability that an individual isolates
                                       onset_to_quarantine,                # symptom onset to quarantine time distribution
                                       quarantine_efficacy,                # effectiveness of the quarantine at reducing onwards transmission
                                       seed                                # stochastic seed
                                       ) {

  ## Setting the seed
  set.seed(seed)

  ## Setting up the offspring distribution
  offspring_fun <- function(n, susc) {
    new_mn <- mn_offspring * susc/population
    size <- new_mn/(disp_offspring - 1)
    truncdist::rtrunc(n, spec = "nbinom", b = susc, mu = new_mn, size = size)
  }
  susc <- population - initial_immune - 1L ## number of susceptibles remaining

  # Pre-allocate a dataframe with the maximum number of infections that we will simulate
  tdf <- data.frame(
    id = integer(check_final_size),                                   # unique identifier for each infected individual
    parent = integer(check_final_size),                               # the parent infection of the infected individual
    generation = integer(check_final_size),                           # the generation (of infections) that the infected individual belongs to
    time_infection_absolute = NA_real_,                               # time of infection, relative to the start of the simulation
    time_infection_relative_parent = NA_real_,                        # time of infection, relative to the parent's time of infection
    time_onset_absolute = NA_real_,                                   # time of symptom onset, relative to the start of the simulation
    time_onset_relative_parent = NA_real_,                            # time of symptom onset, relative to the parent's time of infection
    vaccinated = integer(check_final_size),                           # whether or not the individual is vaccinated via ring-vaccination at their first opportunity (i.e. as a tertiary infection)
    vaccinated_as_secondary = NA_real_,                               # whether the individual was vaccinated via ring vaccination as a secondary infection
    vaccinated_as_tertiary = NA_real_,                                # whether the individual was vaccinated via ring vaccination as a tertiary infection
    vaccinated_2nd_chance = integer(check_final_size),                # whether or not the individual is vaccinated via ring-vaccination at their first opportunity (i.e. as a secondary infection)
    overwritten = NA_real_,                                           # whether or not individual's info was overwritten in that second chance ring vaccination
    time_vaccinated = numeric(check_final_size),                      # if the individual is vaccinated, the time (relative to the start of the simulation) they are vaccinated
    vaccinated_before_infection = integer(check_final_size),          # binary indicator for whether the individual would have received the vaccination before they are infected
    vaccinated_after_infection = integer(check_final_size),           # binary indicator for whether the individual would have received the vaccination after they are infected
    time_protected = numeric(check_final_size),                       # if the individual is vaccinated, the time (relative to the start of the simulation) protection due to vaccination develops
    protected_before_infection = integer(check_final_size),           # binary indicator for whether the individual would have developed protection before they are infected
    protected_after_infection = integer(check_final_size),            # binary indicator for whether the individual would have developed protection after they are infected
    asymptomatic = integer(check_final_size),                         # binary indicator for whether the infected individual is asymptomatic
    quarantined = integer(check_final_size),                          # binary indicator for whether the infected individual quarantines
    time_quarantined_relative_time_infection = NA_real_,              # if the individual quarantines, when do they quarantine (relative to the time of their infection)
    time_quarantined_relative_time_onset = NA_real_,                  # if the individual quarantines, when do they quarantine (relative to the time of their symptom onset)
    time_quarantined_absolute = NA_real_,                             # if the individual quarantines, when do they quarantine (relative to the start of the epidemic)
    n_offspring = integer(check_final_size),                          # number of secondary infections initially generated for this individual
    n_offspring_new = integer(check_final_size),                      # number of secondary infections taking into account (optional) reduced transmissibility of breakthrough infections in vaccinated individuals
    n_offspring_quarantine = integer(check_final_size),               # number of secondary infections after isolation/quarantining
    n_offspring_post_pruning = integer(check_final_size),             # number of secondary infections after ring-vaccination
    n_offspring_post_pruning_2nd_chance = integer(check_final_size),  # number of secondary infections after ring-vaccination second attempt
    secondary_offspring_generated = FALSE,                            # indicator tracking whether secondary infections have been generated for a particular infected individual
    tertiary_offspring_generated = FALSE,                             # indicator tracking whether tertiary infections have been generated for a particular infected individual
    stringsAsFactors = FALSE)

  # Initialize the dataframe with the seeding cases
  tdf[1:seeding_cases, ] <- data.frame(
    id = seq_len(seeding_cases),
    parent = NA_integer_,
    generation = 1L,
    time_infection_absolute = t0 + seq(from = 0, to = 0.01, length.out = seeding_cases),
    time_infection_relative_parent = 0.0,
    time_onset_absolute = NA,
    time_onset_relative_parent = NA,
    vaccinated = 0,
    vaccinated_as_secondary = NA_real_,
    vaccinated_as_tertiary = NA_real_,
    vaccinated_2nd_chance = 0,
    overwritten = 0,
    time_vaccinated = NA,
    vaccinated_before_infection = NA,
    vaccinated_after_infection = NA,
    time_protected = NA,
    protected_before_infection = NA,
    protected_after_infection = NA,
    asymptomatic = integer(seeding_cases),
    quarantined = integer(seeding_cases),
    time_quarantined_relative_time_infection = NA,
    time_quarantined_relative_time_onset = NA,
    time_quarantined_absolute = NA,
    n_offspring = NA_integer_,
    n_offspring_new = NA_integer_,
    n_offspring_quarantine = NA_integer_,
    n_offspring_post_pruning = NA_integer_,
    n_offspring_post_pruning_2nd_chance = NA_integer_,
    secondary_offspring_generated = FALSE,
    tertiary_offspring_generated = FALSE)

  ##########################################################################################################################################
  # Generating secondary and tertiary infections from an index case, then pruning the resulting transmission tree via ring-vaccination
  ##########################################################################################################################################
  while ((any(is.na(tdf$n_offspring)) | any(tdf$n_offspring_post_pruning != 0 & tdf$secondary_offspring_generated == FALSE)) & nrow(tdf) <= check_final_size) {

    ## Get the total number of infections in the dataframe currently (so we can figure out how to label the new infections)
    current_max_row <- max(which(!is.na(tdf$time_infection_absolute)))
    current_max_id <- tdf$id[current_max_row]

    ## Getting the timings of the earliest/oldest infection we haven't yet generated tertiary infections for
    ### NOTE: IF WE WANT TO TOGGLE 1 VS 2 RINGS, WE'LL HAVE TO CHANGE THIS SECTION TO HAVE A DIFFERENT APPROACH AS tertiary_offspring_generated WILL ALWAYS
    ### BE 0 WHEN WE ONLY HAVE 1 RING (AS WE'LL ONLY GENERATE SECONDARY INFECTIONS EACH TIME).
    index_time_infection <- min(tdf$time_infection_absolute[tdf$tertiary_offspring_generated == 0 &  ## NOTE: DO WE WANT SECONDARY OFFSPRING GENERATED IN HERE AS A CONDITION AS WELL? I THINK WE DON'T NEED TO.
                                                              !is.na(tdf$time_infection_absolute)])                  # timing of the earliest unsimulated infection
    index_idx <- which(tdf$time_infection_absolute == index_time_infection & !tdf$tertiary_offspring_generated)[1]   # get the row of the earliest unsimulated infection
    index_id <- tdf$id[index_idx]                                                                                    # id of the earliest unsimulated infection
    index_t <- tdf$time_infection_absolute[index_idx]                                                                # infection time of the earliest unsimulated infection
    index_gen <- tdf$generation[index_idx]                                                                           # generation of the earliest unsimulated infection
    index_vaccinated <- tdf$vaccinated[index_idx]                                                                    # whether or not the index case (the "parent") is vaccinated
    index_time_vaccinated <- tdf$time_vaccinated[index_idx]                                                          # when the index case (the "parent") was vaccinated
    index_time_protected <- tdf$time_protected[index_idx]                                                            # when the index case (the "parent") was protected
    index_onset_time <- infection_to_onset(n = 1)                                                                    # generate the time from infection to symptom onset for the index case
    index_asymptomatic <- tdf$asymptomatic[index_idx]                                                                # whether or not the index case (the "parent") is asymptomatic (influences whether contacts get ring vaccinated or not)
    index_secondary_offspring_generated <- tdf$secondary_offspring_generated[index_idx]                              # whether or not secondary offspring have already been generated for this infection
    ## NOTE: needs to change to make probability of quarantining conditional on having symptoms
    index_quarantine <- rbinom(n = 1, size = 1, prob = prob_quarantine)                                              # whether or not the index case isolates
    index_quarantine_time <- ifelse(index_quarantine == 1, onset_to_quarantine(n = 1), NA)                           # if the infection isolates, how soon after symptom onset they do so

    ## Adding the onset times and quarantine times to this infection's information in the dataframe
    tdf$time_onset_relative_parent[index_idx] <- index_onset_time
    index_onset_time_absolute <- index_onset_time + index_time_infection
    tdf$time_onset_absolute[index_idx] <- index_onset_time_absolute
    tdf$quarantined[index_idx] <- index_quarantine     ## NOTE: needs to change to make probability of quarantining conditional on having symptoms
    tdf$time_quarantined_relative_time_onset[index_idx] <- index_quarantine_time
    tdf$time_quarantined_relative_time_infection[index_idx] <- index_onset_time + index_quarantine_time
    tdf$time_quarantined_absolute[index_idx] <- index_time_infection + index_onset_time + index_quarantine_time

    ## Calculating time to vaccinate the contacts of this index infection
    time_to_contact_vaccination <- index_onset_time + vaccine_logistical_delay                         ## Time between index case infected and contacts being ring vaccinated
    time_to_contact_vaccination_protection <- time_to_contact_vaccination + vaccine_protection_delay   ## Time between index case infected and contacts being protected by the vaccination

    ##########################################################################################################################
    # Generating secondary infections from the index infection
    # - Note: If we've already generated the secondary offspring for this index infection (when they were a secondary
    #         infection to a preceding index infection), we only generate tertiary offspring for this index infection
    #         and skip the step directly below.
    ##########################################################################################################################

    ############################################################################################
    # Generate secondary offspring if we haven't yet generated them for this index infection
    ############################################################################################
    if (index_secondary_offspring_generated == FALSE) {

      ##########################################################################################################
      # Generating secondary infections and infection times, taking vaccination status of index into account
      ##########################################################################################################
      index_n_offspring <- offspring_fun(1, susc)      # raw number of secondary infections
      tdf$n_offspring[index_idx] <- index_n_offspring
      tdf$secondary_offspring_generated[index_idx] <- TRUE

      ## If index infection was vaccinated, account for vaccination's effect on transmission in breakthrough infections
      if (index_vaccinated == 1) {
        if (index_time_protected < index_time_infection) {
          index_n_offspring <- sum(rbinom(n = index_n_offspring, size = 1, prob = 1 - vaccine_efficacy_transmission))
        }
      }
      tdf$n_offspring_new[index_idx] <- index_n_offspring # number of secondary infections after accounting for vaccination's effect on transmission in breakthrough infections
      secondary_infection_times <- generation_time(index_n_offspring)

      ###########################################################################################################################################################
      ## Removing any infections that are averted due to quarantining (assumed to occur index_quarantine_time after symptom onset, which is index_onset_time)
      ###########################################################################################################################################################
      ## If index infection quarantines, reduce secondary infections - note that quarantine only occurs if infection has symptoms
      if (index_quarantine == 1 & index_asymptomatic == 0) { ## NOTE: needs to change to make probability of quarantining conditional on having symptoms
        index_quarantine_possible_avert <- ifelse((index_onset_time + index_quarantine_time) < secondary_infection_times, 1, 0)           # if an infection occurs later than quarantining, it can be averted by quarantine
        index_quarantine_averted <- rbinom(n = index_n_offspring, size = 1, prob = index_quarantine_possible_avert * quarantine_efficacy) # is the infection actually averted by the quarantining (which can be imperfect)
        index_n_offspring <- index_n_offspring - sum(index_quarantine_averted)                                                            # removing infections averted by quarantine from index_n_offspring
        secondary_infection_times <- secondary_infection_times[which(index_quarantine_averted == 0)]                                      # removing infections averted by quarantine from secondary_infection_times
        tdf$n_offspring_quarantine[index_idx] <- index_n_offspring                                                                        # number of secondary infections after accounting for vaccination's effect on transmission in breakthrough infections AND quarantine
      ## If no quarantining, no reduction in secondary infections
      } else {
        tdf$n_offspring_quarantine[index_idx] <- index_n_offspring
      }

      ###################################################################################################################################################
      # If this index infection generates secondary infections, calculate their infection times and whether they're prevented by ring vaccination etc
      ###################################################################################################################################################
      if (index_n_offspring > 0) {

        ###################################################################################################################################################################################
        # If the vaccine hasn't yet been deployed or infection is asymptomatic, no secondary infections are prevented by ring vaccination - add secondary infections to the dataframe
        ###################################################################################################################################################################################
        if ((index_onset_time_absolute < vaccine_start) | index_asymptomatic == 1) {

          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "id"] <- c(current_max_id + seq_len(index_n_offspring))
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "parent"] <- index_id
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "generation"] <- index_gen + 1L
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "time_infection_absolute"] <- secondary_infection_times + index_t
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "time_infection_relative_parent"] <- secondary_infection_times
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "time_onset_absolute"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "time_onset_relative_parent"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "vaccinated"] <- 0
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "vaccinated_as_secondary"] <- 0
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "vaccinated_as_tertiary"] <- 0
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "vaccinated_2nd_chance"] <- 0
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "overwritten"] <- 0
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "time_vaccinated"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "vaccinated_before_infection"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "vaccinated_after_infection"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "time_protected"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "protected_before_infection"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "protected_after_infection"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "asymptomatic"] <- rbinom(n = index_n_offspring, size = 1, prob = prop_asymptomatic) # whether these secondary infections are asymptomatic
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "n_offspring"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "n_offspring_new"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "n_offspring_quarantine"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "n_offspring_post_pruning"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "n_offspring_post_pruning_2nd_chance"] <- NA
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "secondary_offspring_generated"] <- FALSE
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "tertiary_offspring_generated"] <- FALSE

          index_pruned_n_offspring <- index_n_offspring # number of secondary infections after accounting for vaccination's effect on transmission in breakthrough infections AND quarantine AND ring vaccination (which doesn't occur in this case)
          tdf$n_offspring_post_pruning[index_idx] <- index_pruned_n_offspring

        ###################################################################################################################################################################################
        # If ring-vaccination can occur, calculate timings of secondary infections relative to ring vaccination, assess which infections are prevented, prune the transmission tree
        # and add the remaining infections (i.e. those NOT averted by ring-vaccination) to the dataframe
        ###################################################################################################################################################################################
        } else {

          ## Creating vectors to store whether or not secondary infections are successfully ring-vaccinated, protected and/or prevented
          secondary_vaccinated_successfully <- vector(mode = "integer", length = index_n_offspring)  ## vector of whether secondary infections get vaccinated in time
          secondary_protected_successfully <- vector(mode = "integer", length = index_n_offspring)   ## vector of whether secondary infections get protected in time
          secondary_infection_retained <- vector(mode = "integer", length = index_n_offspring)       ## vector of whether secondary infections are retained (i.e. not prevented by ring vaccination)
          secondary_infection_retained[1:length(secondary_infection_retained)] <- 1                  ## default to infections being retained; and then flow through below to see if they get removed

          ## Looping over secondary infections and evaluating whether they're vaccinated, protected and/or prevented
          for (i in 1:index_n_offspring) {

            # If infection would otherwise occur AFTER (potential) ring-vaccination
            if (time_to_contact_vaccination <= secondary_infection_times[i]) {
              secondary_vaccinated_successfully[i] <- rbinom(n = 1, size = 1, prob = vaccine_coverage)                                    # are they vaccinated?
              secondary_protected_successfully[i] <- ifelse(time_to_contact_vaccination_protection <= secondary_infection_times[i], 1, 0) # does protection develop in time?
              successful_protection <- rbinom(n = 1, size = 1, prob = vaccine_efficacy_infection * secondary_protected_successfully[i])   # does protection successfully stop infection?
              secondary_infection_retained[i] <- ifelse(secondary_vaccinated_successfully[i] == 1 & successful_protection == 1, 0, 1)     # is the infection retained (i.e. ring-vaccination fails to avert)

            # If infection would otherwise occur BEFORE ring-vaccination
            } else {
              secondary_vaccinated_successfully[i] <- 0 ## note that we're eliding together "unvaccinated" and "vaccinated after infection occurs" here
              secondary_infection_retained[i] <- 1
            }
          }

          ## Pruning the secondary infections after applying ring-vaccination
          index_pruned_n_offspring <- sum(secondary_infection_retained)          # number of secondary infections after accounting for ring vaccination, vaccination's effect on transmission in breakthrough infections AND quarantine
          tdf$n_offspring_post_pruning[index_idx] <- index_pruned_n_offspring
          retained_index <- which(secondary_infection_retained == 1)                                                                             # which secondary infections were NOT averted by ring-vaccination and thus are retained for inclusion in the dataframe
          secondary_pruned_infection_times <- secondary_infection_times[retained_index]                                                          # infection times of retained infections
          secondary_vaccinated <- ifelse(secondary_vaccinated_successfully[retained_index] == 0, 0, 1)                                           # of the retained infections, which are vaccinated
          secondary_time_vaccinated <- ifelse(secondary_vaccinated_successfully[retained_index] == 1, time_to_contact_vaccination, NA)           # of the retained infections, when are they vaccinated (relative to infection time of index)
          secondary_vaccinated_before_infection <- ifelse(secondary_vaccinated_successfully[retained_index] == 0, NA, 1)                         # (currently we combine all individuals not vaccinated and vaccinated after infection into "unvaccinated", so all vaccinated individuals necessarily got vaccinated before infection)
          secondary_time_protected <- ifelse(secondary_vaccinated_successfully[retained_index] == 1, time_to_contact_vaccination_protection, NA) # of the retained infections, when are they protected (relative to infection time of index)
          secondary_protected_before_infection <- ifelse(is.na(secondary_time_protected), NA, ifelse(secondary_time_protected <= secondary_pruned_infection_times, 1, 0))  # retained infections were protected by vaccine before infection (i.e. vaccine protection successfully developed but failed to protect)
          secondary_protected_after_infection <- ifelse(is.na(secondary_time_protected), NA, ifelse(secondary_time_protected > secondary_pruned_infection_times, 1, 0))    # retained infections were NOT protected by vaccine before infection (i.e. vaccine protection did not develop in time to protect)

          ## Adding secondary infections that aren't prevented by ring-vaccination to the overall dataframe
          if (index_pruned_n_offspring != 0) {

            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "id"] <- c(current_max_id + seq_len(index_pruned_n_offspring))
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "parent"] <- index_id
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "generation"] <- index_gen + 1L
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "time_infection_absolute"] <- secondary_infection_times[retained_index] + index_t
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "time_infection_relative_parent"] <- secondary_infection_times[retained_index]
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "time_onset_absolute"] <- NA
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "time_onset_relative_parent"] <- NA
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "vaccinated"] <- secondary_vaccinated
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "vaccinated_as_secondary"] <- secondary_vaccinated
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "vaccinated_as_tertiary"] <- 0
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "vaccinated_2nd_chance"] <- 0
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "overwritten"] <- 0
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "time_vaccinated"] <- secondary_time_vaccinated + index_t
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "vaccinated_before_infection"] <- secondary_vaccinated_before_infection
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "vaccinated_after_infection"] <- ifelse(is.na(secondary_vaccinated_before_infection), NA, 0)
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "time_protected"] <- secondary_time_protected + index_t
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "protected_before_infection"] <- secondary_protected_before_infection
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "protected_after_infection"] <- secondary_protected_after_infection
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "asymptomatic"] <- rbinom(n = index_pruned_n_offspring, size = 1, prob = prop_asymptomatic) # whether these secondary infections are asymptomatic
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "n_offspring"] <- NA
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "n_offspring_new"] <- NA
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "n_offspring_quarantine"] <- NA
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "n_offspring_post_pruning"] <- NA
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "n_offspring_post_pruning_2nd_chance"] <- NA
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "secondary_offspring_generated"] <- FALSE
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "tertiary_offspring_generated"] <- FALSE

          }
        }

      ###################################################################################################################################################
      # If this index infection does not end up generating secondary infections, skip all those steps
      ###################################################################################################################################################
      } else {
        index_pruned_n_offspring <- 0
        tdf$n_offspring_post_pruning[index_idx] <- 0
        tdf$n_offspring_post_pruning_2nd_chance[index_idx] <- 0
      }

    ##################################################################################################################################################################################
    # If we have already generated secondary offspring for the index infection, calculate whether any of these infections are removed by the second chance at being ring-vaccinated
    # - Note: If we've already generated secondary offspring for this index (when the now index was itself a secondary infection, and the now secondary offspring were tertiary
    #         infections), they've already had a chance at being ring-vaccinated. That first attempt failed to prevent them being infected (hence their inclusion
    #         and presence in the dataframe). HOWEVER, if they weren't vaccinated during the first chance at being ring-vaccinated (when they were the tertiary of an index
    #         i.e. a contact of a contact), they now have another opportunity to be ring-vaccinated (as a contact of their index infection, which was the secondary infection of
    #         the previous index infection). We only give them this second chance if they weren't vaccinated in the first round - those vaccinated but where
    #         vaccination failed to prevent infection are not re-vaccinated.
    ##################################################################################################################################################################################
    } else if (index_secondary_offspring_generated == TRUE) {

      index_pruned_n_offspring_pre_second_chance <- tdf$n_offspring_post_pruning[index_idx] # number of infections generated by this index, taking into account first chance at ring-vaccination that has already happened

      ## Only simulate the second attempt at ring-vaccination if there are remaining infections to potentially avert
      if (index_pruned_n_offspring_pre_second_chance > 0) {

        ###################################################################################################################################################
        # If this index infection generates secondary infections and these weren't averted during the first chance at ring-vaccination,
        # calculate whether they're prevented by the second chance at ring vaccination
        ## - Note that these infections will already have had an opportunity to be ring-vaccinated - when they were a tertiary infection associated
        ##   with an index infection of a generation older than our current index infection (they're now secondary infections of an index infection
        ##   of a newer generation) - if they were NOT vaccinated during that last round of ring vaccination, we calculate here whether they are
        ##   successfully ring-vaccinated this time around.
        ###################################################################################################################################################

        ##############################################################################
        # Only do this is if the vaccine is available and the index is symptomatic
        ##############################################################################
        if ((index_onset_time_absolute > vaccine_start) & index_asymptomatic == 0) {

          ## Creating vectors to store whether or not infections are vaccinated, protected and/or prevented
          secondary_vaccinated_successfully_2nd_chance <- vector(mode = "integer", length = index_pruned_n_offspring_pre_second_chance)  ## vector of whether secondary infections get vaccinated in time on the 2nd attempt at ring vaccination
          secondary_protected_successfully_2nd_chance <- vector(mode = "integer", length = index_pruned_n_offspring_pre_second_chance)   ## vector of whether secondary infections get protected in time on the 2nd attempt at ring vaccination
          secondary_infection_retained_2nd_chance <- vector(mode = "integer", length = index_pruned_n_offspring_pre_second_chance)       ## vector of whether secondary infections are retained (i.e. not prevented by ring vaccination) on the 2nd attempt at ring vaccination
          secondary_infection_retained_2nd_chance[1:length(secondary_infection_retained_2nd_chance)] <- 1                                ## default to infections being retained; and then flow through below to see if they get removed on the 2nd attempt at ring vaccination

          ## Getting the secondary infection IDs, times and vaccination status (from their first chance at ring vaccination when they were tertiary infections)
          secondary_infection_ids <- tdf$id[tdf$parent == index_id & !is.na(tdf$parent)]
          secondary_infection_times <- tdf$time_infection_relative_parent[tdf$parent == index_id & !is.na(tdf$parent)]
          secondary_vaccination_status <- tdf$vaccinated[tdf$parent == index_id & !is.na(tdf$parent)] # whether those individuals were vaccinated at the first opportunity they had for ring-vaccination (i.e. when they were a tertiary infection)

          if (index_pruned_n_offspring_pre_second_chance != length(secondary_infection_ids)) {
            stop("something's gone wrong with second round of pruning")
          }

          ## Looping over secondary infections and evaluating whether they're vaccinated, protected and/or prevented in this second chance at ring-vaccination
          for (i in 1:index_pruned_n_offspring_pre_second_chance) {

            # If infection would otherwise occur AFTER (potential) vaccination AND they were not vaccinated during the first attempt at ring-vaccination
            if (time_to_contact_vaccination <= secondary_infection_times[i] & secondary_vaccination_status[i] != 1) { # we don't re-vaccinate folks who were vaccinated as tertiary infections but where vaccination was unsuccessful in averting the infection
              secondary_vaccinated_successfully_2nd_chance[i] <- rbinom(n = 1, size = 1, prob = vaccine_coverage)                                                        # are they vaccinated?
              secondary_protected_successfully_2nd_chance[i] <- ifelse(time_to_contact_vaccination_protection <= secondary_infection_times[i], 1, 0)                     # does protection develop in time?
              successful_protection_2nd_chance <- rbinom(n = 1, size = 1, prob = vaccine_efficacy_infection * secondary_protected_successfully_2nd_chance[i])            # does protection successfully stop infection?
              secondary_infection_retained_2nd_chance[i] <- ifelse(secondary_vaccinated_successfully_2nd_chance[i] == 1 & successful_protection_2nd_chance == 1, 0, 1)   # is the infection retained (i.e. ring-vaccination second attempt fails to avert)

            # If infection would otherwise occur BEFORE vaccination OR they were already vaccinated in the previous round (when they were tertiary infections)
            } else {
              secondary_vaccinated_successfully_2nd_chance[i] <- 0 ## note that we're eliding together "unvaccinated" and "vaccinated after infection occurs"
              secondary_infection_retained_2nd_chance[i] <- 1
            }
          }

          ## Pruning the secondary infections after applying the second attempt ring-vaccination
          index_pruned_n_offspring <- sum(secondary_infection_retained_2nd_chance)
          tdf$n_offspring_post_pruning_2nd_chance[tdf$id == index_id & !is.na(tdf$parent)] <- index_pruned_n_offspring

          ## Removing the infections averted in this second ring vaccination attempt
          removed_index_2nd_chance <- which(secondary_infection_retained_2nd_chance == 0) # which secondary infections are averted by ring-vaccination second attempt
          removed_id_2nd_chance <- secondary_infection_ids[removed_index_2nd_chance]      # id of secondary infections averted by ring-vaccination second attempt
          if (sum(removed_id_2nd_chance) != 0) {
            tdf <- tdf[!(tdf$id %in% removed_id_2nd_chance), ]                            # removing these secondary infections from the main dataframe
          }

          ## Calculating the updated information for those infections vaccinated but NOT averted in this second attempt
          retained_index_2nd_chance <- which(secondary_infection_retained_2nd_chance == 1 & secondary_vaccination_status != 1) # Subsetting to only not averted infections vaccinated in THIS round of ring-vaccination - ignoring those from the previous round of ring-vaccination as they've already been vaccinated and aren't considered here
          secondary_pruned_infection_times_2nd_chance <- secondary_infection_times[retained_index_2nd_chance]                                                                      # infection times of retained infections
          secondary_vaccinated_2nd_chance <- ifelse(secondary_vaccinated_successfully_2nd_chance[retained_index_2nd_chance] == 0, 0, 1)                                            # of the retained infections, whether they are vaccinated in this 2nd attempt
          secondary_time_vaccinated_2nd_chance <- ifelse(secondary_vaccinated_successfully_2nd_chance[retained_index_2nd_chance] == 1, time_to_contact_vaccination, NA)            # of the retained infections, when they are vaccinated in this 2nd attempt
          secondary_vaccinated_before_infection_2nd_chance <- ifelse(secondary_vaccinated_successfully_2nd_chance[retained_index_2nd_chance] == 0, NA, 1)                          # of the retained infections, whether they are vaccinated before infection occurs
          secondary_time_protected_2nd_chance <- ifelse(secondary_vaccinated_successfully_2nd_chance[retained_index_2nd_chance] == 1, time_to_contact_vaccination_protection, NA)  # of the retained infections, when are they protected (relative to infection time of index)
          secondary_protected_before_infection_2nd_chance <- ifelse(is.na(secondary_time_protected_2nd_chance), NA, ifelse(secondary_time_protected_2nd_chance <= secondary_pruned_infection_times_2nd_chance, 1, 0))
          secondary_protected_after_infection_2nd_chance <- ifelse(is.na(secondary_time_protected_2nd_chance), NA, ifelse(secondary_time_protected_2nd_chance > secondary_pruned_infection_times_2nd_chance, 1, 0))

          ## Modifying the information in the infections successfully ring vaccinated (but not averted) in this second attempt to take new vaccination info into account
          retained_id_2nd_chance <- secondary_infection_ids[retained_index_2nd_chance]   # id of secondary infections newly ring-vaccination in the second attempt but which were not averted
          index_for_overwrite <- which(tdf$id %in% retained_id_2nd_chance)
          if (!(identical(integer(0), index_for_overwrite))) { # only overwrite if there's actually stuff to overwrite
            tdf[index_for_overwrite, "vaccinated"] <- secondary_vaccinated_2nd_chance
            tdf[index_for_overwrite, "vaccinated_as_secondary"] <- secondary_vaccinated_2nd_chance
            tdf[index_for_overwrite, "vaccinated_2nd_chance"] <- secondary_vaccinated_2nd_chance
            tdf[index_for_overwrite, "time_vaccinated"] <- secondary_time_vaccinated_2nd_chance + index_t
            tdf[index_for_overwrite, "vaccinated_before_infection"] <- secondary_vaccinated_before_infection_2nd_chance
            tdf[index_for_overwrite, "vaccinated_after_infection"] <- ifelse(is.na(secondary_vaccinated_before_infection_2nd_chance), NA, 0)
            tdf[index_for_overwrite, "time_protected"] <- secondary_time_protected_2nd_chance + index_t
            tdf[index_for_overwrite, "protected_before_infection"] <- secondary_protected_before_infection_2nd_chance
            tdf[index_for_overwrite, "protected_after_infection"] <- secondary_protected_after_infection_2nd_chance
            tdf[index_for_overwrite, "overwritten"] <- 1
          }

        ########################################################################################################
        ## If vaccine isn't available or the infection is asymptomatic, then no extra infections are averted
        ########################################################################################################
        } else {
          index_pruned_n_offspring <- index_pruned_n_offspring_pre_second_chance
          tdf$n_offspring_post_pruning_2nd_chance[tdf$id == index_id & !is.na(tdf$parent)] <- index_pruned_n_offspring_pre_second_chance
        }

      ###################################################################################################################################################
      # If this index infection does not generate secondary infections, skip all those steps
      ###################################################################################################################################################
      } else {
        index_pruned_n_offspring <- 0
        tdf$n_offspring_post_pruning_2nd_chance[tdf$id == index_id & !is.na(tdf$parent)] <- 0
      }

    }

    ##########################################################################################################################
    # Generating tertiary infections for each of the secondary infections
    # - Note: Tertiary infections are the offspring of secondary infections i.e. they are the offspring arising from the
    #         offspring of the index infection. They are therefore contacts of the contacts of the index infection
    #         and are therefore eligible for ring-vaccination in response to detection of the index infection. We
    #         simulate them here, establish whether or not they're successfully vaccinated, and prune the transmission tree
    #         as appropriate
    ##########################################################################################################################

    ############################################################################################
    # Generate tertiary offspring if there are secondary offspring to generate them
    ############################################################################################

    ## NOTE: I think adding a condition here like "if (second_ring == TRUE)" and have it being skipped otherwise would enable us to toggle between 1 ring and 2 rings?
    ##       Would require some small changes else (e.g. lines at the very beginning that determines which index case we pick to start simulating from) but overall don't think
    ##       it would be that complicated.
    if (index_pruned_n_offspring != 0) {

      ## Looping through each secondary infection and generating tertiary offspring from them
      for (i in 1:index_pruned_n_offspring) {

        current_max_row <- max(which(!is.na(tdf$time_infection_absolute)))  # total number of infections in the dataframe currently (so we can figure out how to label the new infections)
        current_max_id <- tdf$id[current_max_row]

        ## For the index infection being considered, getting the timings of the earliest/oldest secondary infection that we haven't yet generated tertiary infections for
        secondary_time_infection <- tdf$time_infection_absolute[which(tdf$parent == index_id)][i]                                   # timing of the secondary infections of the index infection being considered
        secondary_idx <- which(tdf$time_infection_absolute == secondary_time_infection & !tdf$secondary_offspring_generated)[1]     # get the id of the earliest secondary infection we haven't generated tertiary infections for
                                                                                                                                    # (note: tertiary is relative to index, for the secondary infection they are their secondary offspring)
        secondary_id <- tdf$id[secondary_idx]                                                                                       # id of the earliest unsimulated secondary infection
        secondary_t <- tdf$time_infection_absolute[secondary_idx]                                                                   # infection time of the earliest unsimulated secondary infection
        secondary_gen <- tdf$generation[secondary_idx]                                                                              # generation of the earliest unsimulated secondary infection
        secondary_vaccinated <- tdf$vaccinated[secondary_idx]                                                                       # whether this secondary infection got vaccinated during ring-vaccination attempts
        secondary_time_vaccinated <- tdf$time_vaccinated[secondary_idx]                                                             # when this secondary infection being considered got vaccinated
        secondary_time_protected <- tdf$time_protected[secondary_idx]                                                               # when this secondary infection being considered got protected by vaccination
        secondary_onset_time <- infection_to_onset(n = 1)                                                                           # generate the time from infection to symptom onset for this secondary infection
        secondary_asymptomatic <- tdf$asymptomatic[secondary_idx]                                                                   # whether or not this secondary infection is asymptomatic (influences whether their contacts/offspring can get ring vaccinated)
        secondary_time_infection_relative <- tdf$time_infection_relative[secondary_idx]                                             # time of infection relative to the parent (the index infection)
        secondary_quarantine <- rbinom(n = 1, size = 1, prob = prob_quarantine)                                                     # whether or not this secondary infection quarantines
        secondary_time_quarantine <- ifelse(secondary_quarantine == 1, onset_to_quarantine(n = 1), NA)                              # if they quarantine, when do they quarantine (how long after symptom onset)
        ## NOTE: needs to change to make probability of quarantining conditional on having symptoms

        ## Adding the onset times to this secondary infection's information in the dataframe
        tdf$time_onset_relative_parent[secondary_idx] <- secondary_onset_time
        secondary_onset_time_absolute <- secondary_onset_time + secondary_time_infection
        tdf$time_onset_absolute[secondary_idx] <- secondary_onset_time_absolute
        tdf$quarantined[secondary_idx] <- secondary_quarantine     ## NOTE: needs to change to make probability of quarantining conditional on having symptoms
        tdf$time_quarantined_relative_time_onset[secondary_idx] <- secondary_time_quarantine
        tdf$time_quarantined_relative_time_infection[secondary_idx] <- secondary_onset_time + secondary_time_quarantine
        tdf$time_quarantined_absolute[secondary_idx] <- secondary_time_infection + secondary_onset_time + secondary_time_quarantine

        ###################################################################################################################
        # Generating infections and infection times, taking vaccination status of index into account
        # - Note: Relative to the index infection being considered, these are tertiary offspring (because they are
        #         the offspring of the index's offspring). But relative to the secondary infections (i.e. the index's
        #         offspring), they are themselves secondary offspring.
        ###################################################################################################################
        secondary_n_offspring <- offspring_fun(1, susc)            # number of offspring this particular secondary infection produces
        tdf$n_offspring[secondary_idx] <- secondary_n_offspring
        if (secondary_vaccinated == 1) {
          if (secondary_time_protected < secondary_time_infection) {
            secondary_n_offspring <- sum(rbinom(n = secondary_n_offspring, size = 1, prob = 1 - vaccine_efficacy_transmission))
          }
        }
        tdf$n_offspring_new[secondary_idx] <- secondary_n_offspring
        tdf$secondary_offspring_generated[secondary_idx] <- TRUE
        tertiary_infection_times <- generation_time(secondary_n_offspring)

        ################################################################################################################################################################
        ## Removing any infections that are averted due to quarantining (assumed to occur secondary_time_quarantine after symptom onset, which is secondary_onset_time)
        ################################################################################################################################################################
        if (secondary_quarantine == 1 & secondary_asymptomatic == 0) {     ## NOTE: needs to change to make probability of quarantining conditional on having symptoms
          secondary_quarantine_possible_avert <- ifelse((secondary_onset_time + secondary_time_quarantine) < tertiary_infection_times, 1, 0)            # if infection occurs later than quarantining, it can be averted by quarantine
          secondary_quarantine_averted <- rbinom(n = secondary_n_offspring, size = 1, prob = secondary_quarantine_possible_avert * quarantine_efficacy) # is the infection actually averted by the quarantining (which can be imperfect)
          secondary_n_offspring <- secondary_n_offspring - sum(secondary_quarantine_averted)                                                            # removing infections averted by quarantine from index_n_offspring
          tertiary_infection_times <- tertiary_infection_times[which(secondary_quarantine_averted == 0)]                                                # removing infections averted by quarantine from secondary_infection_times
          tdf$n_offspring_quarantine[secondary_idx] <- secondary_n_offspring                                                                            # number of secondary infections after accounting for vaccination's effect on transmission in breakthrough infections AND quarantine
        ## If no quarantining, no reduction in secondary infections
        } else {
          tdf$n_offspring_quarantine[secondary_idx] <- secondary_n_offspring
        }

        ###################################################################################################################################################
        # If this secondary infection generates tertiary infections, calculate their infection times and whether they're prevented by ring vaccination etc
        ###################################################################################################################################################
        if (secondary_n_offspring > 0) {

          ###################################################################################################################################################################################
          # If the vaccine hasn't yet been deployed or infection is asymptomatic, no tertiary infections are prevented by ring vaccination - add tertiary infections to the dataframe
          ###################################################################################################################################################################################
          if ((index_onset_time_absolute < vaccine_start) | index_asymptomatic == 1) {

            tdf$n_offspring_post_pruning[secondary_idx] <- secondary_n_offspring

            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "id"] <- c(current_max_id + seq_len(secondary_n_offspring))
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "parent"] <- secondary_id
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "generation"] <- secondary_gen + 1L
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "time_infection_absolute"] <- tertiary_infection_times + secondary_t
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "time_infection_relative_parent"] <- tertiary_infection_times
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "time_onset_absolute"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "time_onset_relative_parent"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "vaccinated"] <- 0
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "vaccinated_as_secondary"] <- 0
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "vaccinated_as_tertiary"] <- 0
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "vaccinated_2nd_chance"] <- 0
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "overwritten"] <- 0
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "time_vaccinated"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "vaccinated_before_infection"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "vaccinated_after_infection"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "time_protected"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "protected_before_infection"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "protected_after_infection"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "asymptomatic"] <- rbinom(n = secondary_n_offspring, size = 1, prob = prop_asymptomatic)
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "n_offspring"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "n_offspring_new"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "n_offspring_quarantine"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "n_offspring_post_pruning"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "n_offspring_post_pruning_2nd_chance"] <- NA
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "secondary_offspring_generated"] <- FALSE
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "tertiary_offspring_generated"] <- FALSE

          ###################################################################################################################################################################################
          # If ring-vaccination can occur, calculate timings of tertiary infections relative to ring vaccination, assess which infections are prevented, prune the transmission tree
          # and add the remaining infections (i.e. those NOT averted by ring-vaccination) to the dataframe
          ###################################################################################################################################################################################
          } else {

            ## Creating vectors to store whether or not infections are vaccinated, protected and/or prevented
            tertiary_vaccinated_successfully <- vector(mode = "integer", length = secondary_n_offspring)          ## vector of whether tertiary infections get vaccinated in time
            tertiary_protected_successfully <- vector(mode = "integer", length = secondary_n_offspring)           ## vector of whether tertiary infections get protected in time
            tertiary_infection_retained <- vector(mode = "integer", length = secondary_n_offspring)               ## vector of whether tertiary infections are retained (i.e. not prevented by ring vaccination)
            tertiary_infection_retained[1:length(tertiary_infection_retained)] <- 1                               ## default to infections being retained; and then flow through below to see if they get removed

            ## Looping over secondary infections and evaluating whether they're vaccinated, protected and/or prevented
            for (i in 1:secondary_n_offspring) {

              # If infection would otherwise occur AFTER (potential) vaccination (note these are tertiary infections so we have to account for timing of secondary infections given that ring vaccination
              # is activated in response to the index infection)
              if (time_to_contact_vaccination <= (secondary_time_infection_relative + tertiary_infection_times[i])) {
                tertiary_vaccinated_successfully[i] <- rbinom(n = 1, size = 1, prob = vaccine_coverage)                                                                          # are they vaccinated?
                tertiary_protected_successfully[i] <- ifelse(time_to_contact_vaccination_protection <= (secondary_time_infection_relative + tertiary_infection_times[i]), 1, 0)  # does protection develop in time?
                successful_protection <- rbinom(n = 1, size = 1, prob = vaccine_efficacy_infection * tertiary_protected_successfully[i])                                         # does protection successfully stop infection?
                tertiary_infection_retained[i] <- ifelse(tertiary_vaccinated_successfully[i] == 1 & successful_protection == 1, 0, 1)                                            # is the infection retained (i.e. ring-vaccination fails to avert)

              # If infection would otherwise occur BEFORE vaccination
              } else {
                tertiary_vaccinated_successfully[i] <- 0 ## note that we're eliding together "unvaccinated" and "vaccinated after infection occurs" here
                tertiary_infection_retained[i] <- 1
              }
            }

            ## Pruning the secondary infections after applying ring-vaccination
            secondary_pruned_n_offspring <- sum(tertiary_infection_retained)             # number of tertiary infections remaining after ring vaccination
            tdf$n_offspring_post_pruning[secondary_idx] <- secondary_pruned_n_offspring
            retained_index <- which(tertiary_infection_retained == 1)                    # which tertiary infections are not averted by ring vaccination (i.e. they should be retained to be included in the infections dataframe)
            tertiary_pruned_infection_times <- tertiary_infection_times[retained_index]  # infection times of the tertiary infections not averted by ring vaccination
            tertiary_vaccinated <- ifelse(tertiary_vaccinated_successfully[retained_index] == 0, 0, 1)                                            # of the retained infections, which are vaccinated
            tertiary_time_vaccinated <- ifelse(tertiary_vaccinated_successfully[retained_index] == 1, time_to_contact_vaccination, NA)            # of the retained infections, when are they vaccinated (relative to infection time of index)
            tertiary_vaccinated_before_infection <- ifelse(tertiary_vaccinated_successfully[retained_index] == 0, NA, 1)                          # (currently we combine all individuals not vaccinated and vaccinated after infection into "unvaccinated", so all vaccinated individuals necessarily got vaccinated before infection)
            tertiary_time_protected <- ifelse(tertiary_vaccinated_successfully[retained_index] == 1, time_to_contact_vaccination_protection, NA)  # of the retained infections, when are they protected (relative to infection time of index)
            tertiary_protected_before_infection <- ifelse(is.na(tertiary_time_protected), NA, ifelse(tertiary_time_protected <= tertiary_pruned_infection_times, 1, 0))
            tertiary_protected_after_infection <- ifelse(is.na(tertiary_time_protected), NA, ifelse(tertiary_time_protected > tertiary_pruned_infection_times, 1, 0))

            ## Adding tertiary infections that aren't prevented by ring-vaccination to the overall dataframe
            if (secondary_pruned_n_offspring != 0) {

              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "id"] <- c(current_max_id + seq_len(secondary_pruned_n_offspring))
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "parent"] <- secondary_id
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "generation"] <- secondary_gen + 1L
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "time_infection_absolute"] <- tertiary_infection_times[retained_index] + secondary_t
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "time_infection_relative_parent"] <- tertiary_infection_times[retained_index]
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "time_onset_absolute"] <- NA
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "time_onset_relative_parent"] <- NA
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "vaccinated"] <- tertiary_vaccinated
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "vaccinated_as_secondary"] <- 0
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "vaccinated_as_tertiary"] <- tertiary_vaccinated
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "vaccinated_2nd_chance"] <- 0
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "overwritten"] <- 0
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "time_vaccinated"] <- tertiary_time_vaccinated + index_t
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "vaccinated_before_infection"] <- tertiary_vaccinated_before_infection
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "vaccinated_after_infection"] <- ifelse(is.na(tertiary_vaccinated_before_infection), NA, 0)
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "time_protected"] <- tertiary_time_protected + index_t
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "protected_before_infection"] <- tertiary_protected_before_infection
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "protected_after_infection"] <- tertiary_protected_after_infection
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "asymptomatic"] <- rbinom(n = secondary_pruned_n_offspring, size = 1, prob = prop_asymptomatic)
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "n_offspring"] <- NA
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "n_offspring_new"] <- NA
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "n_offspring_quarantine"] <- NA
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "n_offspring_post_pruning"] <- NA
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "n_offspring_post_pruning_2nd_chance"] <- NA
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "secondary_offspring_generated"] <- FALSE
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "tertiary_offspring_generated"] <- FALSE
            }
          }
        ###################################################################################################################################################
        # If this index infection does not end up generating tertiary infections after accounting for ring vaccination
        ###################################################################################################################################################
        } else {
          tdf$n_offspring_post_pruning[secondary_idx] <- 0
        }
      }
    }
    ###################################################################################################################################################
    # Having generated tertiary offspring, update the index infection's entry to reflect this and adjust the susceptible population appropriately
    ###################################################################################################################################################
    tdf$tertiary_offspring_generated[index_idx] <- TRUE
    susc <- population - sum(tdf$n_offspring_post_pruning, na.rm = TRUE) # Note: this is crude - might need to do something more nuanced
  }
  tdf <- tdf[order(tdf$time_infection_absolute, tdf$id), ]
  return(tdf)
}

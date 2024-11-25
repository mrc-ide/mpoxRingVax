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
## To Consider Building In:
##
## - Make prob_contact_traced a function of the number of active infections (or something that means it's time varying and related to the size of the epidemic)
## - Explicit tracking of the number of vaccine doses required - will require consideration of contacts that DON'T otherwise become infections (could be done post-hoc instead)
## - Limits to the number of doses and its impact on probability of someone being successfully vaccinated
## - Ability to specify whether 1 ring or 2 rings of ring-vaccination
## - Ability to specify whether vaccine coverage is due to being missed by health system (in which case 2nd attempt might work) or due to hesitancy/refusal (in which case 2nd attemtp won't work)
##
#########################################################################################################################################################################################################
#' @export
mpox_ring_vaccination_sim <- function(## Sexual Transmission Parameters
                                      mn_offspring_sexual_SW,                 # mean of the sexual offspring distribution for SW
                                      disp_offspring_sexual_SW,               # overdispersion of the sexual offspring distribution for SW
                                      mn_offspring_sexual_PBS,                # mean of the sexual offspring distribution for PBS
                                      disp_offspring_sexual_PBS,              # mean of the sexual offspring distribution for PBS
                                      mn_offspring_sexual_genPop,             # mean of the sexual offspring distribution for genPop
                                      disp_offspring_sexual_genPop,           # mean of the sexual offspring distribution for genPop
                                      sexual_transmission_occupation_matrix,  # rows are the proportions of cases in each occupation a particular occupation gives rises to

                                      ## Household Transmission Parameters
                                      mn_offspring_hh,                        # mean of the household offspring distribution
                                      disp_offspring_hh,                      # overdispersion of the household offspring distribution

                                      ## Community Transmission Parameters
                                      mn_offspring_community,                 # mean of the community offspring distribution
                                      disp_offspring_community,               # overdispersion of the community offspring distribution
                                      community_transmission_age_matrix,      # rows are the proportions of cases in each age-group a particular age-group gives rises to (note this needs to take into account smallpox vaccination)

                                      ## Natural History Parameters
                                      prop_asymptomatic,                      # probability an infection is asymptomatic
                                      generation_time,                        # generation time distribution
                                      infection_to_onset,                     # infection to symptom onset time distribution

                                      ## Vaccine-Related Parameters
                                      vaccine_start,                          # time at which the vaccine becomes available
                                      vaccine_coverage,                       # probability that each eligible individual gets vaccinated
                                      vaccine_efficacy_infection,             # vaccine efficacy against infection
                                      vaccine_efficacy_transmission,          # reduction in transmissibility of breakthrough infections in vaccinated individuals
                                      vaccine_logistical_delay,               # delay between symptom onset and vaccination of contacts (and contacts of contacts) occurring
                                      vaccine_protection_delay,               # delay between vaccination and protection developing

                                      ## Quarantine Related Parameters
                                      prob_quarantine,                        # probability that an individual isolates
                                      onset_to_quarantine,                    # symptom onset to quarantine time distribution
                                      quarantine_efficacy_household,          # effectiveness of the quarantine at reducing onwards transmission for household infections
                                      quarantine_efficacy_else,               # effectiveness of the quarantine at reducing onwards transmission for infections from other sources

                                      ## Miscellaneous Parameters
                                      synthetic_household_df,                 # dataframe of synthetic drc household, age and occupation population
                                      t0 = 0,                                 # simulation starting time
                                      population,                             # total population size
                                      check_final_size,                       # maximum number of infections to simulate
                                      initial_immune = 0,                     # initial number of individuals who are immune
                                      seeding_cases,                          # number of seeding cases to start the epidemic with
                                      seed                                    # stochastic seed
                                      ) {

  ## Setting the seed
  set.seed(seed)

  ## Setting up the number of susceptibles
  susc <- population - initial_immune

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
    transmission_route = NA_character_,                               # how the individual was infected - sexual, household or community
    occupation = NA_character_,                                       # characteristic of the individual - sex worker (SW), person who buys sex (PBS) or general population (genPop)
    age = NA_character_,                                              # age group of the individual (0-5, 5-18, 18+) - SW and PBS can only be 18+
    hh_id = integer(check_final_size),                                # id of the household an individual resides in
    hh_member_index = integer(check_final_size),                      # id number for each individual within a household
    hh_size = integer(check_final_size),                              # size of the household
    hh_ages = I(vector("list", check_final_size)),                    # ages of all the individuals in the household
    hh_occupations = I(vector("list", check_final_size)),             # occupation of all the individuals in the household
    hh_infections = NA_real_,                                         # total number of infections that have occurred in the household (note: this is dynamically updated)
    hh_infected_index = I(vector("list", check_final_size)),          # id numbers of all individuals in the household who have been infected (note: this is dynamically updated)
    stringsAsFactors = FALSE)

  # Sampling the households for the initial seeding cases - we seed this epidemic in sex workers (i.e. all initial seeding cases are assumed to be SWs)
  sampled_rows <- synthetic_household_df[get("contains_SW") == 1][sample(.N, seeding_cases)]  # sampling from our synthetic DRC age/occ/hh dataframe households that contains a SW
  seeding_case_indices <- unlist(lapply(sampled_rows$hh_occupations, function(occupations) {  # Find the household member id of the individual who is a SW
    which(occupations == "SW")   # only ever 1 SW or in a household, so this works
  }))
  seeding_case_occupations <- rep("SW", seeding_cases)
  seeding_case_ages <- rep("18+", seeding_cases)        ## all SWs are assumed to be 18+

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
    tertiary_offspring_generated = FALSE,
    transmission_route = NA_character_,
    occupation = seeding_case_occupations,
    age = seeding_case_ages,
    hh_id = 1:seeding_cases,
    hh_member_index = seeding_case_indices,
    hh_size = sampled_rows$hh_size,
    hh_ages = I(sampled_rows$hh_ages_grouped),
    hh_occupations = I(sampled_rows$hh_occupations),
    hh_infections = 1,
    hh_infected_index = I(as.list(seeding_case_indices)))

  ##########################################################################################################################################
  # Generating secondary and tertiary infections from an index case, then pruning the resulting transmission tree via ring-vaccination
  ##########################################################################################################################################
  # We use the words INDEX, SECONDARY and TERTIARY to indicate how infections are related to each other. SECONDARY offspring refers to
  # direct offspring of a particular INDEX infection. TERTIARY offspring are offspring of those SECONDARY offspring, i.e. they are
  # indirect descendants of the INDEX infection, separated by a generation. We use these terms relationally i.e. the TERTIARY offspring
  # of a particular INDEX infection will ALSO be the SECONDARY offspring of that INDEX infection's SECONDARY offspring.
  #
  # During each loop of the "while" loop below, we select an infection for whom we have not yet generated the full complement of offspring
  # (i.e. both SECONDARY and TERTIARY offspring). This is what index_secondary_offspring_generated and index_tertiary_offspring_generated
  # track.
  #
  # If SECONDARY offspring have NOT yet been generated for this particular INDEX infection:
  #
  #  Then by definition TERTIARY offspring have not been generated either. We therefore begin by using the offspring function
  #  to generate SECONDARY offspring. What follows then is a series of decision points and potential pruning steps where potential offspring
  #  are averted. These steps depend on the timing of the generated SECONDARY infections, and properties of the INDEX infection:
  #
  #  - If the INDEX infection had previously been vaccinated and we are modelling an impact of vaccination in reducing transmissibility
  #    in breakthrough infections of vaccinated individuals, we prune the offspring dataframe to reflect this reduced transmissibility.
  #    This pruning does NOT depend on the timings of the SECONDARY offspring generated by the INDEX infection.
  #
  #  - If the INDEX infection is symptomatic (i.e. knows they're infected) and is quarantining, then we further prune the offspring
  #    dataframe to reflect the infections averted by quarantining. We model different impacts of quarantining on household infections
  #    vs those arising from other sources. This pruning DOES depend on the timings of the SECONDARY offspring generated by the
  #    INDEX infection. Only infections that would otherwise be generated AFTER the infection quarantines have a chance to be averted.
  #
  #  - If the INDEX infection is symptomatic and the vaccine is available, then we implement a round of ring vaccination. During this
  #    round of ring vaccination, and prune the transmission tree to reflect infections averted by ring vaccination. This pruning
  #    DOES depend on the timings of the SECONDARY offspring generated by the INDEX infection. Only infections who would otherwise
  #    be generated after they are both vaccinated and protected by vaccination have a chance to be averted.
  #
  #  - If after all of this, there are still SECONDARY offspring remaining (i.e. they haven't all been averted), we then generate TERTIARY
  #    offspring for each of the SECONDARY offspring. We flow through the same series of decision points as above (prior vaccination,
  #    quarantining and ring-vaccination) and again prune the transmission tree as relevant.
  #
  # If SECONDARY offspring HAVE ALREADY been generated for this particular INDEX infection:
  #
  #  Remember, we use INDEX, SECONDARY and TERTIARY relationally, i.e. the TERTIARY offspring of a particular INDEX infection will
  #  ALSO be the SECONDARY offspring of that INDEX infection's SECONDARY offspring. What this means is that when we select a new INDEX
  #  infection for whom we have not yet generated both SECONDARY and TERTIARY offspring for, we have usually already generated their
  #  SECONDARY offspring (when they were a SECONDARY offspring themself for a prior INDEX infection). If that is the case, we skip the
  #  steps above (which were implemented before) and instead do the following:
  #
  #  - If the new INDEX infection is symptomatic and the vaccine is available, we implement a round of ring vaccination. Note that
  #    this represents a second chance for these SECONDARY offspring at being ring vaccinated. All infections in this simulation
  #    framework get 2 potential chances at being ring-vaccinated. The first chance is when they are TERTIARY offspring of an INDEX
  #    infection (i.e. in the 2nd "ring" of the ring vaccination campaign"); and the second chance is when they are SECONDARY offspring
  #    of an INDEX infection (i.e. in the 1st "ring" of another ring vaccination campaign - which is NOW) - this new INDEX infection will
  #    have been the SECONDARY offspring of the INDEX infection in a previous step. We then prune the transmission tree further based on
  #    infections averted during this 2nd chance at ring vaccination.
  #
  #  - If there are still SECONDARY offspring remaining after this (i.e. they haven't all been averted by this 2nd chance at ring
  #    vaccination), we then generate TERTIARY offspring for each of the SECONDARY offspring. We flow through the same series of decision
  #    points as above (prior vaccination, quarantining and ring-vaccination) and again prune the transmission tree as relevant.
  #
  # We keep repeating this process i.e. i) select a new infection to be the INDEX, ii) generate SECONDARY AND/OR TERTIARY offspring (depending
  # on what has already been generated), iii) prune the resulting transmission tree connecting INDEX, SECONDARY and TERTIARY offspring based
  # on control measures; until either the outbreak has gone extinct (there are no more infections to generate offspring for) OR
  # we hit a certain number of simulated infections (this is what check_final_size controls).
  #
  ##########################################################################################################################################

  ## While we haven't hit the simulation cap size (check_final_size) and any infections exist where we have not yet generated the requisite offspring, repeat
  ## the loop described in the text above
  while ((any(is.na(tdf$n_offspring)) | any(tdf$n_offspring_post_pruning != 0 & tdf$secondary_offspring_generated == FALSE)) & nrow(tdf) <= check_final_size) {

    ## Get the total number of infections in the dataframe currently (so we can figure out how to label the new infections)
    current_max_row <- max(which(!is.na(tdf$time_infection_absolute)))
    current_max_id <- tdf$id[current_max_row]
    current_max_hh_id <- max(tdf$hh_id)

    ## Getting the timings of the earliest/oldest infection we haven't yet generated tertiary infections for - this is the "INDEX INFECTION"
    index_time_infection <- min(tdf$time_infection_absolute[tdf$tertiary_offspring_generated == 0 &
                                                              !is.na(tdf$time_infection_absolute)])                  # timing of the earliest unsimulated infection ("index" infection)
    index_idx <- which(tdf$time_infection_absolute == index_time_infection & !tdf$tertiary_offspring_generated)[1]   # get the row of the earliest unsimulated infection ("index" infection)
    index_id <- tdf$id[index_idx]                                                                                    # id of the earliest unsimulated infection ("index" infection)
    index_t <- tdf$time_infection_absolute[index_idx]                                                                # infection time of the earliest unsimulated infection ("index" infection)
    index_gen <- tdf$generation[index_idx]                                                                           # generation of the earliest unsimulated infection ("index" infection)
    index_vaccinated <- tdf$vaccinated[index_idx]                                                                    # whether or not the index infection (the "parent") is vaccinated
    index_time_vaccinated <- tdf$time_vaccinated[index_idx]                                                          # when the index infection (the "parent") was vaccinated
    index_time_protected <- tdf$time_protected[index_idx]                                                            # when the index infection (the "parent") was protected
    index_onset_time <- infection_to_onset(n = 1)                                                                    # generate the time from infection to symptom onset for the index infection (this is the time relative to index's parent infection)
    index_asymptomatic <- tdf$asymptomatic[index_idx]                                                                # whether or not the index infection (the "parent") is asymptomatic (influences whether contacts get ring vaccinated or not)
    index_secondary_offspring_generated <- tdf$secondary_offspring_generated[index_idx]                              # whether or not secondary offspring have already been generated for this infection
    index_quarantine <- rbinom(n = 1, size = 1, prob = prob_quarantine)                                              # whether or not the index infection isolates
    index_quarantine_time <- ifelse(index_quarantine == 1, onset_to_quarantine(n = 1), NA)                           # if the infection isolates, how soon after symptom onset they do so
    index_occupation <- tdf$occupation[index_idx]                                                                    # occupation of the index infection (SW, PBS or genPop)
    index_age_group <- tdf$age[index_idx]                                                                            # age-group of the index infection (0-5, 5-18 or 18+)
    index_hh_id <- tdf$hh_id[index_idx]                                                                              # household ID of the index infection
    index_hh_member_index <- tdf$hh_member_index[index_idx]                                                          # house member ID for each household (used to subset the ages and occupations vectors for each household)
    index_hh_size <- unlist(tdf$hh_size[index_idx])                                                                  # household size
    index_hh_ages <- unlist(tdf$hh_ages[index_idx])                                                                  # ages of all of the household members
    index_hh_occupations <- unlist(tdf$hh_occupations[index_idx])                                                    # occupations of all the household members
    index_hh_infections <- tdf$hh_infections[index_idx]                                                              # cumulative number of infections there have been in this particular household (Note: need to make sure this is updated for the index case when we've simulated from them)
    index_hh_infected_index <- unlist(tdf$hh_infected_index[index_idx])                                              # house member IDs of all infected household members

    ## Adding the onset times and quarantine times to this infection's information in the dataframe
    tdf$time_onset_relative_parent[index_idx] <- index_onset_time                                                    # adding index onset time (relative to index's parent)to the storage dataframe
    index_onset_time_absolute <- index_onset_time + index_time_infection                                             # converting index_onset_time (symptom onset time relative to index's parent) into absolute calendar time
    tdf$time_onset_absolute[index_idx] <- index_onset_time_absolute                                                  # adding index onset time (absolute calendar time) to the storage dataframe
    tdf$quarantined[index_idx] <- index_quarantine                                                                   # adding quarantine indicator to storage dataframe
    tdf$time_quarantined_relative_time_onset[index_idx] <- index_quarantine_time                                     # adding quarantine time relative to index's symptom onset to the storage dataframe
    tdf$time_quarantined_relative_time_infection[index_idx] <- index_onset_time + index_quarantine_time              # adding quarantine time relative to index's infection to the storage dataframe
    tdf$time_quarantined_absolute[index_idx] <- index_time_infection + index_onset_time + index_quarantine_time      # adding quarantine time in absolute calendar time to the storage dataframe

    ## Calculating time to vaccinate the contacts of this index infection
    time_to_contact_vaccination <- index_onset_time + vaccine_logistical_delay                                       # time between index case infected and contacts being ring vaccinated
    time_to_contact_vaccination_protection <- time_to_contact_vaccination + vaccine_protection_delay                 # time between index case infected and contacts being protected by the vaccination


    ##########################################################################################################################
    # Generating secondary infections from the index infection
    # - Note: If we've already generated the secondary offspring for this index infection (when they were a secondary
    #         infection to a preceding index infection), we skip the step directly below (which is about generating
    #         secondary offspring) and instead see whether any of the previously generated secondary offspring
    #         are prevented by the second chance at ring vaccination offered by having two rings of ring-vaccination.
    ##########################################################################################################################

    ############################################################################################
    # Generate secondary offspring if we haven't yet generated them for this index infection
    ############################################################################################
    if (index_secondary_offspring_generated == FALSE) {

      ##########################################################################################################
      # Generating secondary infections for this index infection
      ##########################################################################################################
      index_offspring_function_draw <- offspring_fun(synthetic_household_df = synthetic_household_df,
                                                     index_age_group = index_age_group,
                                                     index_hh_id = index_hh_id,
                                                     index_occupation = index_occupation,
                                                     max_hh_id = current_max_hh_id,

                                                     mn_offspring_sexual_SW = mn_offspring_sexual_SW,
                                                     disp_offspring_sexual_SW = disp_offspring_sexual_SW,
                                                     mn_offspring_sexual_PBS = mn_offspring_sexual_PBS,
                                                     disp_offspring_sexual_PBS = disp_offspring_sexual_PBS,
                                                     mn_offspring_sexual_genPop = mn_offspring_sexual_genPop,
                                                     disp_offspring_sexual_genPop = disp_offspring_sexual_genPop,
                                                     sexual_transmission_occupation_matrix = sexual_transmission_occupation_matrix,

                                                     mn_offspring_hh = mn_offspring_hh,
                                                     disp_offspring_hh = disp_offspring_hh,
                                                     index_hh_member_index = index_hh_member_index,
                                                     index_hh_infections_index = index_hh_infected_index,
                                                     index_hh_size = index_hh_size,
                                                     index_hh_prior_infections = index_hh_infections,
                                                     index_hh_ages = index_hh_ages,
                                                     index_hh_occupations = index_hh_occupations,

                                                     mn_offspring_community = mn_offspring_community,
                                                     disp_offspring_community = disp_offspring_community,
                                                     number_susceptible_community = susc,
                                                     population_community = population,
                                                     community_transmission_age_matrix = community_transmission_age_matrix)

      index_n_offspring <- index_offspring_function_draw$total_offspring
      tdf$n_offspring[index_idx] <- index_n_offspring
      tdf$secondary_offspring_generated[index_idx] <- TRUE

      ####################################################################################################################################
      # Removing secondary infections after taking vaccination's effect on transmission in breakthrough infections into account
      ####################################################################################################################################
      ## Only flow through this step if index infection was vaccinated and still infected (i.e. a breakthrough infection)
      if (index_vaccinated == 1) {

        ## Only do this step if the index's vaccination protection developed BEFORE they were infected (otherwise they wouldn't have any effect of vaccination)
        if (index_time_protected < index_time_infection) {

          ## Binomial draw to decide which infections are averted by the reduced transmissibility of the breakthrough infection
          index_offspring_vaccine_retained <- rbinom(n = index_n_offspring, size = 1, prob = 1 - vaccine_efficacy_transmission)

          ## Subsetting index_offspring_function_draw to get the infections that don't occur because of the reduced transmissibility of breakthrough infections
          index_offspring_vaccine_averted_index <- which(index_offspring_vaccine_retained == 0)
          index_offspring_vaccine_averted_transmission_route <- index_offspring_function_draw$offspring_characteristics$transmission_route[index_offspring_vaccine_averted_index] # getting the transmisssion route of each of the averted infections
          index_offspring_vaccine_averted_hh_member_index <- index_offspring_function_draw$offspring_characteristics$hh_member_index[index_offspring_vaccine_averted_index]       # getting the hh member ids of each of the averted infections

          ## Subsetting index_offspring_function_draw to reflect the averted infections because of reduced transmissibility of the breakthrough infection
          index_offspring_retained_index <- which(index_offspring_vaccine_retained == 1)

          ## Updating the offspring function draw index_offspring_function_draw to reflect loss of averted infections
          index_offspring_function_draw$offspring_characteristics <- index_offspring_function_draw$offspring_characteristics[index_offspring_retained_index, ]    # subsetting the offspring dataframe to only retain the infections that weren't averted
          index_offspring_function_draw$total_offspring <- length(index_offspring_retained_index)                                                                 # updating total number of offspring to reflect loss of averted infections
          index_offspring_function_draw$num_offspring_sexual <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "sexual")       # updating number of sexual infections to reflect loss of averted infections
          index_offspring_function_draw$num_offspring_hh <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "household")        # updating number of household infections to reflect loss of averted infections
          index_offspring_function_draw$num_offspring_community <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "community") # updating number of community infections to reflect loss averted infections

          ## If any of the averted infections are household infections, we ALSO need to update the information of any infections who share their household
          ## - Specifically, we need to update the cumulative household infections tracker and the index tracking which household members have been infected
          ##   and remove the averted infections from both of these
          ## - Note that because of the way the offspring function is set up, the "household" infections must all belong to the same household,
          ##   which is the same household as the index infection
          if ("household" %in% index_offspring_vaccine_averted_transmission_route) { # are any of the averted infections ones who got infected in a household?

            ## Removing the averted household infections from the cumulative total
            index_offspring_num_household_infections_averted <- sum(index_offspring_vaccine_averted_transmission_route == "household")
            index_offspring_function_draw$new_hh_cumulative_infections <- index_offspring_function_draw$new_hh_cumulative_infections - index_offspring_num_household_infections_averted
            index_offspring_function_draw$offspring_characteristics$hh_infections[index_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- index_offspring_function_draw$new_hh_cumulative_infections

            ## Modifying the list of all ids of infected household members to account for the averted infections
            index_offspring_vaccine_averted_household_member_id <- index_offspring_vaccine_averted_hh_member_index[which(index_offspring_vaccine_averted_transmission_route == "household")]
            index_offspring_vaccine_averted_household_member_index_for_removal <- which(unlist(index_offspring_function_draw$new_hh_infected_index) %in% index_offspring_vaccine_averted_household_member_id)
            index_offspring_function_draw$new_hh_infected_index <- list(unlist(index_offspring_function_draw$new_hh_infected_index)[-index_offspring_vaccine_averted_household_member_index_for_removal])
            index_offspring_function_draw$offspring_characteristics$hh_infected_index[index_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- I(index_offspring_function_draw$new_hh_infected_index)

          }

          ## New total number of offspring produced as a result of  reduced transmissibility of the breakthrough infection
          index_n_offspring <- length(index_offspring_retained_index)

        }
      }
      tdf$n_offspring_new[index_idx] <- index_n_offspring
      secondary_infection_times <- generation_time(index_n_offspring)

      ### GOT TO HERE

      ###########################################################################################################################################################
      ## Removing any infections that are averted due to quarantining (assumed to occur index_quarantine_time after symptom onset, which is index_onset_time)
      ###########################################################################################################################################################
      ## If index infection quarantines, reduce secondary infections - note that quarantine only occurs if infection has symptoms
      if (index_quarantine == 1 & index_asymptomatic == 0) {

        ## Check if infections are quarantining and then do binomial draw to decide which infections are averted by it
        index_quarantine_possible_avert <- ifelse((index_onset_time + index_quarantine_time) < secondary_infection_times, 1, 0)           # if an infection occurs later than quarantining, it can be averted by quarantine
        index_quarantine_efficacy <- ifelse(index_offspring_function_draw$offspring_characteristics$transmission_route == "household",    # depending on where transmission occurs, use different quarantining efficacy
                                            quarantine_efficacy_household,
                                            quarantine_efficacy_else)
        index_offspring_quarantine_retained <- rbinom(n = index_n_offspring, size = 1, prob = 1 - (index_quarantine_possible_avert * index_quarantine_efficacy)) # is the infection averted by the quarantining (which can be imperfect)

        ## Subsetting index_offspring_function_draw to get the infections that don't occur because of the quarantining
        index_offspring_quarantine_averted_index <- which(index_offspring_quarantine_retained == 0)
        index_offspring_quarantine_averted_transmission_route <- index_offspring_function_draw$offspring_characteristics$transmission_route[index_offspring_quarantine_averted_index]
        index_offspring_quarantine_averted_hh_member_index <- index_offspring_function_draw$offspring_characteristics$hh_member_index[index_offspring_quarantine_averted_index]

        ## Subsetting index_offspring_function_draw to reflect the averted infections because of quarantining
        index_offspring_retained_index <- which(index_offspring_quarantine_retained == 1)
        index_offspring_function_draw$offspring_characteristics <- index_offspring_function_draw$offspring_characteristics[index_offspring_retained_index, ]
        index_offspring_function_draw$total_offspring <- length(index_offspring_retained_index)
        index_offspring_function_draw$num_offspring_sexual <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "sexual")
        index_offspring_function_draw$num_offspring_hh <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "household")
        index_offspring_function_draw$num_offspring_community <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "community")

        ## If any of the averted infections are household ones, update the info of any remaining household infections there
        if ("household" %in% index_offspring_quarantine_averted_transmission_route) {

          ## Removing the averted household infections from the cumulative total
          index_offspring_num_household_infections_averted <- sum(index_offspring_quarantine_averted_transmission_route == "household")
          index_offspring_function_draw$new_hh_cumulative_infections <- index_offspring_function_draw$new_hh_cumulative_infections - index_offspring_num_household_infections_averted
          index_offspring_function_draw$offspring_characteristics$hh_infections[index_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- index_offspring_function_draw$new_hh_cumulative_infections

          ## Modifying the list of all ids of infected household members to account for the averted infections
          index_offspring_quarantine_averted_household_member_id <- index_offspring_quarantine_averted_hh_member_index[which(index_offspring_quarantine_averted_transmission_route == "household")]
          index_offspring_quarantine_averted_household_member_index_for_removal <- which(unlist(index_offspring_function_draw$new_hh_infected_index) %in% index_offspring_quarantine_averted_household_member_id)
          index_offspring_function_draw$new_hh_infected_index <- list(unlist(index_offspring_function_draw$new_hh_infected_index)[-index_offspring_quarantine_averted_household_member_index_for_removal])
          index_offspring_function_draw$offspring_characteristics$hh_infected_index[index_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- I(index_offspring_function_draw$new_hh_infected_index)

        }

        # Updating index_n_offspring and secondary_infection_times in light of new removals due to quarantining
        index_n_offspring <- sum(index_offspring_quarantine_retained)                                                                     # accounting for infections averted by quarantine from index_n_offspring
        secondary_infection_times <- secondary_infection_times[index_offspring_retained_index]                                            # removing infections averted by quarantine from secondary_infection_times
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

          ## Add information of new infections to the main storage dataframe
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
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "transmission_route"] <- index_offspring_function_draw$offspring_characteristics$transmission_route
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "occupation"] <- index_offspring_function_draw$offspring_characteristics$occupation
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "age"] <- index_offspring_function_draw$offspring_characteristics$age
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "hh_id"] <- index_offspring_function_draw$offspring_characteristics$hh_id
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "hh_member_index"] <- index_offspring_function_draw$offspring_characteristics$hh_member_index
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "hh_size"] <- index_offspring_function_draw$offspring_characteristics$hh_size
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "hh_ages"] <- list(index_offspring_function_draw$offspring_characteristics$hh_ages)
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "hh_occupations"] <- list(index_offspring_function_draw$offspring_characteristics$hh_occupations)
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "hh_infections"] <- index_offspring_function_draw$offspring_characteristics$hh_infections
          tdf[(current_max_row+1):(current_max_row+index_n_offspring), "hh_infected_index"] <- list(index_offspring_function_draw$offspring_characteristics$hh_infected_index)

          # No ring vaccination so index_pruned_n_offspring is same as index_n_offspring
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

          ## Subsetting index_offspring_function_draw to get the infections that don't occur because of the reduced transmissibility of breakthrough infections
          index_offspring_ring_vaccine_averted_index <- which(secondary_infection_retained == 0)
          index_offspring_ring_vaccine_averted_transmission_route <- index_offspring_function_draw$offspring_characteristics$transmission_route[index_offspring_ring_vaccine_averted_index]
          index_offspring_ring_vaccine_averted_hh_member_index <- index_offspring_function_draw$offspring_characteristics$hh_member_index[index_offspring_ring_vaccine_averted_index]

          ## Updating index_offspring_function_draw to reflect the infections averted because of ring vaccination
          index_offspring_function_draw$offspring_characteristics <- index_offspring_function_draw$offspring_characteristics[retained_index, ]
          index_offspring_function_draw$total_offspring <- length(retained_index)
          index_offspring_function_draw$num_offspring_sexual <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "sexual")
          index_offspring_function_draw$num_offspring_hh <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "household")
          index_offspring_function_draw$num_offspring_community <- sum(index_offspring_function_draw$offspring_characteristics$transmission_route == "community")

          ## If any of the averted infections are household ones, update the info of any remaining household infections there
          if ("household" %in% index_offspring_ring_vaccine_averted_transmission_route) {

            ## Removing the averted household infections from the cumulative total
            index_offspring_num_household_infections_averted <- sum(index_offspring_ring_vaccine_averted_transmission_route == "household")
            index_offspring_function_draw$new_hh_cumulative_infections <- index_offspring_function_draw$new_hh_cumulative_infections - index_offspring_num_household_infections_averted
            index_offspring_function_draw$offspring_characteristics$hh_infections[index_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- index_offspring_function_draw$new_hh_cumulative_infections

            ## Modifying the list of all ids of infected household members to account for the averted infections
            index_offspring_ring_vaccine_averted_household_member_id <- index_offspring_ring_vaccine_averted_hh_member_index[which(index_offspring_ring_vaccine_averted_transmission_route == "household")]
            index_offspring_ring_vaccine_averted_household_member_index_for_removal <- which(unlist(index_offspring_function_draw$new_hh_infected_index) %in% index_offspring_ring_vaccine_averted_household_member_id)
            index_offspring_function_draw$new_hh_infected_index <- list(unlist(index_offspring_function_draw$new_hh_infected_index)[-index_offspring_ring_vaccine_averted_household_member_index_for_removal])
            index_offspring_function_draw$offspring_characteristics$hh_infected_index[index_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- I(index_offspring_function_draw$new_hh_infected_index)
            ## CFWNOTE: when there's only 1 hh_infected_index left, this is being coerced to numeric rather than list. Need to check that's not introducing any issues later on.

          }

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
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "transmission_route"] <- index_offspring_function_draw$offspring_characteristics$transmission_route
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "occupation"] <- index_offspring_function_draw$offspring_characteristics$occupation
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "age"] <- index_offspring_function_draw$offspring_characteristics$age
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "hh_id"] <- index_offspring_function_draw$offspring_characteristics$hh_id
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "hh_member_index"] <- index_offspring_function_draw$offspring_characteristics$hh_member_index
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "hh_size"] <- index_offspring_function_draw$offspring_characteristics$hh_size
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "hh_ages"] <- list(index_offspring_function_draw$offspring_characteristics$hh_ages)
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "hh_occupations"] <- list(index_offspring_function_draw$offspring_characteristics$hh_occupations)
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "hh_infections"] <- index_offspring_function_draw$offspring_characteristics$hh_infections
            tdf[(current_max_row+1):(current_max_row+index_pruned_n_offspring), "hh_infected_index"] <- list(index_offspring_function_draw$offspring_characteristics$hh_infected_index)

          }
        }

        ## Update all the previous infections in the household regarding the additional household infections (if there are any)
        if (any(index_offspring_function_draw$offspring_characteristics$transmission_route == "household")) {

          ## Getting all infections previously generated in this household
          hh_row_index <- which(tdf$hh_id == index_hh_id)

          ## Updating the tally of cumulative number of infections in the household for the index case and all other prior infections in the household
          tdf$hh_infections[hh_row_index] <- index_offspring_function_draw$new_hh_cumulative_infections

          ## Updating the list of the indices of all infected individuals in the household for the index case
          tdf$hh_infected_index[hh_row_index] <- index_offspring_function_draw$new_hh_infected_index
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

          ## Removing the infections averted in this second ring vaccination attempt and modifying hh-related contents of the other (retained) household members that are in the same household as averted infections
          removed_index_2nd_chance <- which(secondary_infection_retained_2nd_chance == 0) # which secondary infections are averted by ring-vaccination second attempt and need to be removed
          removed_id_2nd_chance <- secondary_infection_ids[removed_index_2nd_chance]      # id of secondary infections averted by ring-vaccination second attempt and which need to be removed
          removed_hh_id_2nd_chance <- tdf$hh_id[tdf$id %in% removed_id_2nd_chance]        # hh id of secondary infections averted by ring-vaccination second attempt and which need to be removed
          removed_hh_id_2nd_chance_unique <- unique(removed_hh_id_2nd_chance)

          # if(index_id == 10) {
          #   stop("error throwing time")
          # }

          if (sum(removed_id_2nd_chance) != 0) { # only do this if infections are being removed by 2nd chance ring vax

            ## Looping through each of the households with infections removed, finding other infections that share their household, and modifying these infections' info
            ## to reflect removal of this infection
            for (i in 1:length(removed_hh_id_2nd_chance_unique)) { # CFWNOTE: Need to check this is behaving as it should, especially w.r.t to the list modification at the end of this segment

              ## Getting the id of all hh members of that particular household
              removed_temp_hh_id <- removed_hh_id_2nd_chance_unique[i]
              removed_temp_id <- removed_id_2nd_chance[removed_hh_id_2nd_chance == removed_temp_hh_id]   # id of the hh member(s) we're removing
              removed_temp_id_hh_members <- tdf$id[tdf$hh_id == removed_temp_hh_id]                      # ids of the other hh members with that hh id in the dataframe
              removed_temp_id_hh_members_excluding_removed <- removed_temp_id_hh_members[!(removed_temp_id_hh_members %in% removed_temp_id)]
              ### CFWNOTE - wherever else I've done the minus thing, replace with the exclamation point thing

              ## If hh members are removed from a household with other infections already in it, modify the remaining hh members' info to reflect that
              if (length(removed_temp_id) > 0 & length(removed_temp_id_hh_members_excluding_removed) > 0) {

                if (length(unique(tdf$hh_infections[tdf$id %in% removed_temp_id_hh_members_excluding_removed])) != 1) {
                  stop("something wrong with household size and cumulative infections differing across individuals in the same household")
                }
                # print(index_id)

                ## Modifying the cumulative number of household infections recorded for other household infections of this removed infectio
                tdf$hh_infections[tdf$id %in% removed_temp_id_hh_members_excluding_removed] <- tdf$hh_infections[tdf$id %in% removed_temp_id_hh_members_excluding_removed] - length(removed_temp_id) # removing this particular averted infection from the tally of total infected in the household

                ## Modifying the list of the ids of all infected household members of this removed infection
                current_hh_infected_index <- unlist(tdf$hh_infected_index[tdf$id %in% removed_temp_id_hh_members_excluding_removed][1]) # [1] is there because there will potentially be multiple hh members in tdf, all of whom will have the same hh_infected_index so only first is needed here (esp. as all are being overwritten below)
                removed_temp_hh_infected_index <- tdf$hh_member_index[tdf$id %in% removed_temp_id]
                removed_temp_hh_infected_index_for_removal <- which(current_hh_infected_index %in% removed_temp_hh_infected_index)
                updated_hh_infected_index <- current_hh_infected_index[-removed_temp_hh_infected_index_for_removal]
                tdf$hh_infected_index[tdf$id %in% removed_temp_id_hh_members_excluding_removed] <- list(updated_hh_infected_index) # removing this particular averted infection from the list of all infected members in the household and updating the remaining household members' entries

              }

            }

            # Removing the secondary infections averted by their 2nd chance of ring vaccination from the main dataframe
            tdf <- tdf[!(tdf$id %in% removed_id_2nd_chance), ]
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
        current_max_hh_id <- max(tdf$hh_id)

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
        secondary_occupation <- tdf$occupation[secondary_idx]                                                                    # occupation of the secondary infection (SW, PBS or genPop)
        secondary_age_group <- tdf$age[secondary_idx]                                                                            # age-group of the secondary infection (0-5, 5-18 or 18+)
        secondary_hh_id <- tdf$hh_id[secondary_idx]                                                                              # household ID of the secondary infection
        secondary_hh_member_index <- tdf$hh_member_index[secondary_idx]                                                          # house member ID for each household (used to subset the ages and occupations vectors for each household)
        secondary_hh_size <- unlist(tdf$hh_size[secondary_idx])                                                                  # household size
        secondary_hh_ages <- unlist(tdf$hh_ages[secondary_idx])                                                                  # ages of all of the household members
        secondary_hh_occupations <- unlist(tdf$hh_occupations[secondary_idx])                                                    # occupations of all the household members
        secondary_hh_infections <- tdf$hh_infections[secondary_idx]                                                              # cumulative number of infections there have been in this particular household (Note: need to make sure this is updated for the index case when we've simulated from them)
        secondary_hh_infected_index <- unlist(tdf$hh_infected_index[secondary_idx])                                              # house member IDs of all infected household members

        ## Adding the onset times to this secondary infection's information in the dataframe
        tdf$time_onset_relative_parent[secondary_idx] <- secondary_onset_time
        secondary_onset_time_absolute <- secondary_onset_time + secondary_time_infection
        tdf$time_onset_absolute[secondary_idx] <- secondary_onset_time_absolute
        tdf$quarantined[secondary_idx] <- secondary_quarantine
        tdf$time_quarantined_relative_time_onset[secondary_idx] <- secondary_time_quarantine
        tdf$time_quarantined_relative_time_infection[secondary_idx] <- secondary_onset_time + secondary_time_quarantine
        tdf$time_quarantined_absolute[secondary_idx] <- secondary_time_infection + secondary_onset_time + secondary_time_quarantine

        ###################################################################################################################
        # Generating infections and infection times, taking vaccination status of index into account
        # - Note: Relative to the index infection being considered, these are tertiary offspring (because they are
        #         the offspring of the index's offspring). But relative to the secondary infections (i.e. the index's
        #         offspring), they are themselves secondary offspring.
        ###################################################################################################################
        secondary_offspring_function_draw <- offspring_fun(synthetic_household_df = synthetic_household_df,
                                                           index_age_group = secondary_age_group,
                                                           index_hh_id = secondary_hh_id,
                                                           index_occupation = secondary_occupation,
                                                           max_hh_id = current_max_hh_id,

                                                           mn_offspring_sexual_SW = mn_offspring_sexual_SW,
                                                           disp_offspring_sexual_SW = disp_offspring_sexual_SW,
                                                           mn_offspring_sexual_PBS = mn_offspring_sexual_PBS,
                                                           disp_offspring_sexual_PBS = disp_offspring_sexual_PBS,
                                                           mn_offspring_sexual_genPop = mn_offspring_sexual_genPop,
                                                           disp_offspring_sexual_genPop = disp_offspring_sexual_genPop,
                                                           sexual_transmission_occupation_matrix = sexual_transmission_occupation_matrix,

                                                           mn_offspring_hh = mn_offspring_hh,
                                                           disp_offspring_hh = disp_offspring_hh,
                                                           index_hh_member_index = secondary_hh_member_index,
                                                           index_hh_infections_index = secondary_hh_infected_index,
                                                           index_hh_size = secondary_hh_size,
                                                           index_hh_prior_infections = secondary_hh_infections,
                                                           index_hh_ages = secondary_hh_ages,
                                                           index_hh_occupations = secondary_hh_occupations,

                                                           mn_offspring_community = mn_offspring_community,
                                                           disp_offspring_community = disp_offspring_community,
                                                           number_susceptible_community = susc,
                                                           population_community = population,
                                                           community_transmission_age_matrix = community_transmission_age_matrix)
        print("tertiary gen done")

        secondary_n_offspring <- secondary_offspring_function_draw$total_offspring
        tdf$n_offspring[secondary_idx] <- secondary_n_offspring
        tdf$secondary_offspring_generated[secondary_idx] <- TRUE

        ####################################################################################################################################
        # Removing tertiary infections after taking vaccination's effect on transmission in breakthrough infections into account
        #   CFWNote: Need to check, when vaccine_efficacy_transmission is 1, does that introduce NAs that get passed through to next section
        ####################################################################################################################################
        if (secondary_vaccinated == 1) { # Only do this if secondary infection was vaccinated

          if (secondary_time_protected < secondary_time_infection) {

            ## Binomial draw to decide which infections are averted by the reduced transmissibility of the breakthrough infection
            secondary_offspring_retained_index <- rbinom(n = secondary_n_offspring, size = 1, prob = 1 - vaccine_efficacy_transmission)

            ## Subsetting secondary_offspring_function_draw to get the infections that don't occur because of the reduced transmissibility of breakthrough infections
            secondary_offspring_vaccine_averted_index <- which(secondary_offspring_retained_index == 0)
            secondary_offspring_vaccine_averted_transmission_route <- secondary_offspring_function_draw$offspring_characteristics$transmission_route[secondary_offspring_vaccine_averted_index]
            secondary_offspring_vaccine_averted_hh_member_index <- secondary_offspring_function_draw$offspring_characteristics$hh_member_index[secondary_offspring_vaccine_averted_index]

            ## Subsetting secondary_offspring_function_draw to reflect the averted infections because of reduced transmissibility of the breakthrough infection
            secondary_offspring_retained_index <- which(secondary_offspring_retained_index == 1)
            secondary_offspring_function_draw$offspring_characteristics <- secondary_offspring_function_draw$offspring_characteristics[secondary_offspring_retained_index, ]
            secondary_offspring_function_draw$total_offspring <- length(secondary_offspring_retained_index)
            secondary_offspring_function_draw$num_offspring_sexual <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "sexual")
            secondary_offspring_function_draw$num_offspring_hh <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household")
            secondary_offspring_function_draw$num_offspring_community <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "community")

            ## If any of the averted infections are household ones, update the info of any remaining household infections there
            if ("household" %in% secondary_offspring_vaccine_averted_transmission_route) {

              ## Removing the averted household infections from the cumulative total
              secondary_offspring_num_household_infections_averted <- sum(secondary_offspring_vaccine_averted_transmission_route == "household")
              secondary_offspring_function_draw$new_hh_cumulative_infections <- secondary_offspring_function_draw$new_hh_cumulative_infections - secondary_offspring_num_household_infections_averted
              secondary_offspring_function_draw$offspring_characteristics$hh_infections[secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- secondary_offspring_function_draw$new_hh_cumulative_infections

              ## Modifying the list of all ids of infected household members to account for the averted infections
              secondary_offspring_vaccine_averted_household_member_id <- secondary_offspring_vaccine_averted_hh_member_index[which(secondary_offspring_vaccine_averted_transmission_route == "household")]
              secondary_offspring_vaccine_averted_household_member_index_for_removal <- which(unlist(secondary_offspring_function_draw$new_hh_infected_index) %in% secondary_offspring_vaccine_averted_household_member_id)
              secondary_offspring_function_draw$new_hh_infected_index <- list(unlist(secondary_offspring_function_draw$new_hh_infected_index)[-secondary_offspring_vaccine_averted_household_member_index_for_removal])
              secondary_offspring_function_draw$offspring_characteristics$hh_infected_index[secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- I(secondary_offspring_function_draw$new_hh_infected_index)

            }

            ## New total number of offspring produced as a result of  reduced transmissibility of the breakthrough infection
            secondary_n_offspring <- length(secondary_offspring_retained_index)

          }
        }
        tdf$n_offspring_new[secondary_idx] <- secondary_n_offspring
        tertiary_infection_times <- generation_time(secondary_n_offspring)

        ################################################################################################################################################################
        ## Removing any infections that are averted due to quarantining (assumed to occur secondary_time_quarantine after symptom onset, which is secondary_onset_time)
        ################################################################################################################################################################
        ## If secondary infection quarantines, reduce tertiary infections - note that quarantine only occurs if infection has symptoms
        if (secondary_quarantine == 1 & secondary_asymptomatic == 0) {

          ## Check if infections are quarantining and then do binomial draw to decide which infections are averted by it
          secondary_quarantine_possible_avert <- ifelse((secondary_onset_time + secondary_time_quarantine) < tertiary_infection_times, 1, 0)  # if an infection occurs later than quarantining, it can be averted by quarantine
          secondary_quarantine_efficacy <- ifelse(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household",      # depending on where transmission occurs, use different quarantining efficacy
                                                  quarantine_efficacy_household,
                                                  quarantine_efficacy_else)
          secondary_offspring_quarantine_retained <- rbinom(n = secondary_n_offspring, size = 1, prob = 1 - (secondary_quarantine_possible_avert * secondary_quarantine_efficacy)) # is the infection averted by the quarantining (which can be imperfect)

          ## Subsetting secondary_offspring_function_draw to get the infections that don't occur because of the quarantining and their properties
          secondary_offspring_quarantine_averted_index <- which(secondary_offspring_quarantine_retained == 0)
          secondary_offspring_quarantine_averted_transmission_route <- secondary_offspring_function_draw$offspring_characteristics$transmission_route[secondary_offspring_quarantine_averted_index]
          secondary_offspring_quarantine_averted_hh_member_index <- secondary_offspring_function_draw$offspring_characteristics$hh_member_index[secondary_offspring_quarantine_averted_index]

          ## Subsetting secondary_offspring_function_draw to reflect the averted infections because of quarantining
          secondary_offspring_retained_index <- which(secondary_offspring_quarantine_retained == 1)
          secondary_offspring_function_draw$offspring_characteristics <- secondary_offspring_function_draw$offspring_characteristics[secondary_offspring_retained_index, ]
          secondary_offspring_function_draw$total_offspring <- length(secondary_offspring_retained_index)
          secondary_offspring_function_draw$num_offspring_sexual <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "sexual")
          secondary_offspring_function_draw$num_offspring_hh <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household")
          secondary_offspring_function_draw$num_offspring_community <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "community")

          ## If any of the averted infections are household ones, update the info of any remaining household infections there
          if ("household" %in% secondary_offspring_quarantine_averted_transmission_route) {

            ## Removing the averted household infections from the cumulative total
            secondary_offspring_num_household_infections_averted <- sum(secondary_offspring_quarantine_averted_transmission_route == "household")
            secondary_offspring_function_draw$new_hh_cumulative_infections <- secondary_offspring_function_draw$new_hh_cumulative_infections - secondary_offspring_num_household_infections_averted
            secondary_offspring_function_draw$offspring_characteristics$hh_infections[secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- secondary_offspring_function_draw$new_hh_cumulative_infections

            ## Modifying the list of all ids of infected household members to account for the averted infections
            secondary_offspring_quarantine_averted_household_member_id <- secondary_offspring_quarantine_averted_hh_member_index[which(secondary_offspring_quarantine_averted_transmission_route == "household")]
            secondary_offspring_quarantine_averted_household_member_index_for_removal <- which(unlist(secondary_offspring_function_draw$new_hh_infected_index) %in% secondary_offspring_quarantine_averted_household_member_id)
            secondary_offspring_function_draw$new_hh_infected_index <- list(unlist(secondary_offspring_function_draw$new_hh_infected_index)[-secondary_offspring_quarantine_averted_household_member_index_for_removal])
            secondary_offspring_function_draw$offspring_characteristics$hh_infected_index[secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- I(secondary_offspring_function_draw$new_hh_infected_index)

          }

          # Updating secondary_n_offspring and tertiary_infection_times in light of new removals due to quarantining
          secondary_n_offspring <- sum(secondary_offspring_quarantine_retained)                                        # accounting for infections averted by quarantine from index_n_offspring
          tertiary_infection_times <- tertiary_infection_times[secondary_offspring_retained_index]                     # removing infections averted by quarantine from secondary_infection_times
          tdf$n_offspring_quarantine[index_idx] <- secondary_n_offspring                                               # number of secondary infections after accounting for vaccination's effect on transmission in breakthrough infections AND quarantine

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

            ## Add information of new infections to the main storage dataframe
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
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "transmission_route"] <- secondary_offspring_function_draw$offspring_characteristics$transmission_route
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "occupation"] <- secondary_offspring_function_draw$offspring_characteristics$occupation
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "age"] <- secondary_offspring_function_draw$offspring_characteristics$age
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "hh_id"] <- secondary_offspring_function_draw$offspring_characteristics$hh_id
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "hh_member_index"] <- secondary_offspring_function_draw$offspring_characteristics$hh_member_index
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "hh_size"] <- secondary_offspring_function_draw$offspring_characteristics$hh_size
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "hh_ages"] <- list(secondary_offspring_function_draw$offspring_characteristics$hh_ages)
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "hh_occupations"] <- list(secondary_offspring_function_draw$offspring_characteristics$hh_occupations)
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "hh_infections"] <- secondary_offspring_function_draw$offspring_characteristics$hh_infections
            tdf[(current_max_row+1):(current_max_row+secondary_n_offspring), "hh_infected_index"] <- list(secondary_offspring_function_draw$offspring_characteristics$hh_infected_index)

            # No ring vaccination so n_offspring_post_pruning is same as secondary_n_offspring
            tdf$n_offspring_post_pruning[secondary_idx] <- secondary_n_offspring

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

            ## Subsetting secondary_offspring_function_draw to get the infections that don't occur because of the reduced transmissibility of breakthrough infections
            secondary_offspring_ring_vaccine_averted_index <- which(tertiary_infection_retained == 0)
            secondary_offspring_ring_vaccine_averted_transmission_route <- secondary_offspring_function_draw$offspring_characteristics$transmission_route[secondary_offspring_ring_vaccine_averted_index]
            secondary_offspring_ring_vaccine_averted_hh_member_index <- secondary_offspring_function_draw$offspring_characteristics$hh_member_index[secondary_offspring_ring_vaccine_averted_index]

            ## Updating secondary_offspring_function_draw to reflect the infections averted because of ring vaccination
            secondary_offspring_function_draw$offspring_characteristics <- secondary_offspring_function_draw$offspring_characteristics[retained_index, ]
            secondary_offspring_function_draw$total_offspring <- length(retained_index)
            secondary_offspring_function_draw$num_offspring_sexual <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "sexual")
            secondary_offspring_function_draw$num_offspring_hh <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household")
            secondary_offspring_function_draw$num_offspring_community <- sum(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "community")

            ## If any of the averted infections are household ones, update the info of any remaining household infections there
            if ("household" %in% secondary_offspring_ring_vaccine_averted_transmission_route) {

              ## Removing the averted household infections from the cumulative total
              secondary_offspring_num_household_infections_averted <- sum(secondary_offspring_ring_vaccine_averted_transmission_route == "household")
              secondary_offspring_function_draw$new_hh_cumulative_infections <- secondary_offspring_function_draw$new_hh_cumulative_infections - secondary_offspring_num_household_infections_averted
              secondary_offspring_function_draw$offspring_characteristics$hh_infections[secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- secondary_offspring_function_draw$new_hh_cumulative_infections

              ## Modifying the list of all ids of infected household members to account for the averted infections
              secondary_offspring_ring_vaccine_averted_household_member_id <- secondary_offspring_ring_vaccine_averted_hh_member_index[which(secondary_offspring_ring_vaccine_averted_transmission_route == "household")]
              secondary_offspring_ring_vaccine_averted_household_member_index_for_removal <- which(unlist(secondary_offspring_function_draw$new_hh_infected_index) %in% secondary_offspring_ring_vaccine_averted_household_member_id)
              secondary_offspring_function_draw$new_hh_infected_index <- list(unlist(secondary_offspring_function_draw$new_hh_infected_index)[-secondary_offspring_ring_vaccine_averted_household_member_index_for_removal])
              secondary_offspring_function_draw$offspring_characteristics$hh_infected_index[secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household"] <- I(secondary_offspring_function_draw$new_hh_infected_index)
              ## CFWNOTE: when there's only 1 hh_infected_index left, this is being coerced to numeric rather than list. Need to check that's not introducing any issues later on.

            }

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
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "transmission_route"] <- secondary_offspring_function_draw$offspring_characteristics$transmission_route
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "occupation"] <- secondary_offspring_function_draw$offspring_characteristics$occupation
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "age"] <- secondary_offspring_function_draw$offspring_characteristics$age
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "hh_id"] <- secondary_offspring_function_draw$offspring_characteristics$hh_id
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "hh_member_index"] <- secondary_offspring_function_draw$offspring_characteristics$hh_member_index
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "hh_size"] <- secondary_offspring_function_draw$offspring_characteristics$hh_size
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "hh_ages"] <- list(secondary_offspring_function_draw$offspring_characteristics$hh_ages)
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "hh_occupations"] <- list(secondary_offspring_function_draw$offspring_characteristics$hh_occupations)
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "hh_infections"] <- secondary_offspring_function_draw$offspring_characteristics$hh_infections
              tdf[(current_max_row+1):(current_max_row+secondary_pruned_n_offspring), "hh_infected_index"] <- list(secondary_offspring_function_draw$offspring_characteristics$hh_infected_index)
            }
          }

          ## Update all the previous infections in the household regarding the additional household infections (if there are any left after all that pruning)
          if (any(secondary_offspring_function_draw$offspring_characteristics$transmission_route == "household")) {

            ## Getting all infections previously generated in this household
            hh_row_index <- which(tdf$hh_id == secondary_hh_id)

            ## Updating the tally of cumulative number of infections in the household for the index case and all other prior infections in the household
            tdf$hh_infections[hh_row_index] <- secondary_offspring_function_draw$new_hh_cumulative_infections

            ## Updating the list of the indices of all infected individuals in the household for the index case
            tdf$hh_infected_index[hh_row_index] <- secondary_offspring_function_draw$new_hh_infected_index
          }

        ###################################################################################################################################################
        # If this secondary infection does not end up generating tertiary infections after accounting for ring vaccination
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

#' @export
offspring_fun <- function(synthetic_household_df,                    # dataframe of synthetic drc household/age/occ population
                          index_age_group,                           # age of the index infection (0-5, 5-18, 18+)
                          index_hh_id,                               # household id of the index infection
                          index_occupation,                          # occupation of the index infection
                          max_hh_id,                                 # maximum current household id

                          mn_offspring_sexual_SW,                    # mean of the sexual offspring distribution for SW
                          disp_offspring_sexual_SW,                  # overdispersion of the sexual offspring distribution for SW
                          mn_offspring_sexual_PBS,                   # mean of the sexual offspring distribution for PBS
                          disp_offspring_sexual_PBS,                 # mean of the sexual offspring distribution for PBS
                          mn_offspring_sexual_genPop,                # mean of the sexual offspring distribution for genPop
                          disp_offspring_sexual_genPop,              # mean of the sexual offspring distribution for genPop
                          sexual_transmission_occupation_matrix,     # rows are the proportions of cases in each occupation a particular occupation gives rises to

                          mn_offspring_hh,                           # mean of the household offspring distribution
                          disp_offspring_hh,                         # overdispersion of the household offspring distribution
                          index_hh_member_index,                     # index for the index infection in this household
                          index_hh_infections_index,                 # index for all the infections (including the index infection) in this household
                          index_hh_size,                             # size of index infection's household (including the index infection)
                          index_hh_prior_infections,                 # number of prior infections in the index's household (including the index infection)
                          index_hh_ages,                             # ages of index infection's household (including the index infection)
                          index_hh_occupations,                      # occupations of the index infection's household (including the index infection)

                          mn_offspring_community,                    # mean of the community offspring distribution
                          disp_offspring_community,                  # overdispersion of the community offspring distribution
                          number_susceptible_community,              # number of susceptible individuals remaining in the community
                          population_community,                      # overall population size
                          community_transmission_age_matrix,         # rows are the proportions of cases in each age-group a particular age-group gives rises to (note this needs to take into account smallpox vaccination)
) {

  # Helper vectors that will assist with subsetting later on
  age_group_to_column <- c("0-5" = "contains_0_5", "5-18" = "contains_5_18", "18+" = "contains_18plus")
  occupation_to_column <- c("genPop" = "contains_genPop", "PBS" = "contains_PBS", "SW" = "contains_SW")
  age_group_to_index <- c("0-5" = 1, "5-18" = 2, "18+" = 3)
  occupation_to_index <- c("genPop" = 1, "PBS" = 2, "SW" = 3)

  ###################################################################################################################################
  ## SEXUAL INFECTIONS: Generating infections from sexual transmission and their characteristics
  ## - Note that we assume:
  ##     1) All SW and PBS are >18
  ##     2) No sexual transmission to/from <18
  ##     3) Single mpox introduction per household (thus all infections are in newly instantiated and separate households)
  ###################################################################################################################################

  ## Generating sexual transmission offspring

  # Sexual transmission only in those aged >18+, in a manner dependent on their occupation
  if (index_age_group %in% c("18+")) {
    if (index_occupation == "SW") {
      num_offspring_sexual <- rnbinom(n = 1, mu = mn_offspring_sexual_SW, size = disp_offspring_sexual_SW)
    } else if (index_occupation == "PBS") {
      num_offspring_sexual <- rnbinom(n = 1, mu = mn_offspring_sexual_PBS, size = disp_offspring_sexual_PBS)
    } else if (index_occupation == "genPop") {
      num_offspring_sexual <- rnbinom(n = 1, mu = mn_offspring_sexual_genPop, size = disp_offspring_sexual_genPop)
    } else {
      stop("something has gone wrong with the sexual offspring distribution generation of infections")
    }

  ## If individuals are <18, no sexual transmission occurs
  } else {
    num_offspring_sexual <- 0
  }

  ## If sexual offspring are generated, create them and imbue them with all the required characteristics
  if (num_offspring_sexual > 0) {

    offspring_sexual_list <- vector("list", num_offspring_sexual)         # list to temporarily store outputs
    occupation_matrix_index <- occupation_to_index[index_occupation]
    offspring_sexual_occupations <- sample(c("SW", "PBS", "genPop"),      # sample occupations of sexual offspring
                                           size = num_offspring_sexual,
                                           prob = sexual_transmission_occupation_matrix[occupation_matrix_index, ],
                                           replace = TRUE)
    offspring_sexual_new_hh_id <- max_hh_id + 1:num_offspring_sexual # enumerate the household id of each of the sexual offspring (generated in a new household)

    ## Looping through sexual offspring and sampling details of their occupation and household that they belong to, and their occupation
    for (i in 1:num_offspring_sexual) {

      ## Sampling a household for sexual offspring to be produced into, based on the occupation group of the index infection
      offspring_sexual_occupation_group <- offspring_sexual_occupations[i]         # occupation of the offspring
      col_name <- occupation_to_column[offspring_sexual_occupation_group]          # relevant column in the synthetic hh df for sampling from
      sampled_row <- synthetic_household_df[get(col_name) == 1][sample(.N, 1)]     # sampling from synthetic hh df a hh that contains a member of the offspring's occupation

      ## Picking the particular household member infected
      possible_hh_member_indices <- which(unlist(sampled_row$hh_occupations) %in% offspring_sexual_occupation_group)
      if (length(possible_hh_member_indices) == 1) {
        hh_member_index <- possible_hh_member_indices
      } else {
        hh_member_index <- sample(possible_hh_member_indices, 1) # sampling which individual with that occupation in the hh is the infection
      }

      ## Creating dataframe with information on this particular offspring
      offspring_sexual_list[[i]] <- data.table::data.table(transmission_route = "sexual",                                     # transmission route
                                                           occupation = unlist(sampled_row$hh_occupations)[hh_member_index],  # occupation of the offspring
                                                           age = "18+",                                                       # age of the offspring
                                                           hh_id = offspring_sexual_new_hh_id[i],                             # hh id of the new hh that offspring is in
                                                           hh_member_index = hh_member_index,                                 # index of this infection within hh members
                                                           hh_size = sampled_row$hh_size,                                     # hh size
                                                           hh_ages = sampled_row$hh_ages,                                     # ages of individuals in the household
                                                           hh_occupations = sampled_row$hh_occupations,                       # occupations of individuals in the household
                                                           hh_infections = 1,                                                 # cumulative number of infections in that household so far
                                                           hh_infected_index = hh_member_index)                               # index of all infections in this hh that have been infected
    }
    offspring_characteristics_df_sexual <- data.table::rbindlist(offspring_sexual_list, fill = TRUE)

  ## If no sexual offspring generated, initialise an empty sexual offspring df
  } else {
    offspring_characteristics_df_sexual <- data.frame(transmission_route = character(0),
                                                      occupation = character(0),
                                                      age = character(0),
                                                      hh_id = numeric(0),
                                                      hh_member_index = numeric(0),
                                                      hh_size = numeric(0),
                                                      hh_ages = numeric(0),
                                                      hh_occupations = numeric(0),
                                                      hh_infections = numeric(0),
                                                      hh_infected_index = numeric(0))
  }

  #########################################################################################################################################

  ##########################################################################################################################################
  ## HOUSEHOLD INFECTIONS: Generating infections from household transmission and their characteristics
  ## - Note that we assume:
  ##     1) Frequency-dependent transmission e.g. https://elifesciences.org/articles/70767 but might want to change this
  ##########################################################################################################################################
  if (index_hh_prior_infections != length(index_hh_infections_index)) {
    stop("something has gone wrong with keeping track of prior infections")
  }
  number_susceptible_hh <- index_hh_size - length(index_hh_infections_index)   # number of susceptibles remaining in the household
  new_mn_hh <- mn_offspring_hh * number_susceptible_hh/index_hh_size           # mean of household offspring distribution taking susceptible depletion into account
  num_offspring_hh <- min(c(rnbinom(n = 1, mu = new_mn_hh, size = disp_offspring_hh), number_susceptible_hh)) # number of household offspring, capped at number of individuals who can still be infected in the household
  ## note that we're not modifying overdispersion for susceptible depletion - we probably need to, but unclear if the way in the main function is correct - NEED TO CHECK!!!

  ## If household offspring are generated, create them and imbue them with all the required characteristics
  if (num_offspring_hh > 0) {

    ## Sampling which individuals in the household get infected
    offspring_hh_infections_index <- sample(seq.int(1, index_hh_size)[-index_hh_infections_index], num_offspring_hh, replace = FALSE)
    offspring_hh_infections_index <- offspring_hh_infections_index[order(offspring_hh_infections_index)]

    cumulative_hh_infected_index <- c(index_hh_infections_index, offspring_hh_infections_index)
    cumulative_hh_infected_index <- cumulative_hh_infected_index[order(cumulative_hh_infected_index)]

    ## Creating dataframe with information on these new offspring in the household
    offspring_characteristics_df_hh <- data.table::data.table(transmission_route = "household",                                                      # transmission route
                                                              occupation = index_hh_occupations[offspring_hh_infections_index],                      # occupation of new hh infections
                                                              age = index_hh_ages[offspring_hh_infections_index],                                    # age of new hh infections
                                                              hh_id = index_hh_id,                                                                   # hh id of new infections
                                                              hh_member_index = offspring_hh_infections_index,                                       # index of each hh member (i.e. which element in index_hh_ages etc they are)
                                                              hh_size = index_hh_size,                                                               # size of the household
                                                              hh_ages = list(index_hh_ages),                                                         # ages of individuals in the household
                                                              hh_occupations = list(index_hh_occupations),                                           # occupations of individuals in the household
                                                              hh_infections = index_hh_prior_infections + num_offspring_hh,                          # total number of infected individuals in the household
                                                              hh_infected_index = list(cumulative_hh_infected_index))                                # indices of all the infected individuals in the household

  ## If no household offspring generated, initialise an empty household offspring df
  } else {
    offspring_characteristics_df_hh <- data.frame(transmission_route = character(0),
                                                  occupation = character(0),
                                                  age = character(0),
                                                  hh_id = numeric(0),
                                                  hh_member_index = numeric(0),
                                                  hh_size = numeric(0),
                                                  hh_ages = numeric(0),
                                                  hh_occupations = numeric(0),
                                                  hh_infections = numeric(0),
                                                  hh_infected_index = numeric(0))
  }

  ##########################################################################################################################################

  ##########################################################################################################################################
  ## COMMUNITY INFECTIONS: Generating infections from community transmission and their characteristics
  ## - Note that we assume:
  ##     1) Single mpox introduction per household (thus all infections are in newly instantiated and separate households)
  ##########################################################################################################################################
  new_mn_community <- mn_offspring_community * number_susceptible_community/population_community   # mean of community offspring distribution taking susceptible depletion into account
  num_offspring_community <- rnbinom(n = 1, mu = new_mn_community, size = disp_offspring_community)    # number of community offspring generated by index case
  ## note that we're not modifying overdispersion for susceptible depletion - we probably need to, but unclear if the way in the main function is correct - NEED TO CHECK!!!

  ## If community offspring are generated, create them and imbue them with all the required characteristics
  if (num_offspring_community > 0) {

    offspring_community_list <- vector("list", num_offspring_community)  # list to temporarily store outputs
    age_matrix_index <- age_group_to_index[index_age_group]              # getting the relevant row of the community_transmission_age_matrix needed for sampling
    offspring_community_ages <- sample(c("0-5", "5-18", "18+"),          # sample ages of community offspring
                                       size = num_offspring_community,
                                       prob = community_transmission_age_matrix[age_matrix_index, ],
                                       replace = TRUE)
    offspring_community_new_hh_id <- max(offspring_sexual_new_hh_id) + 1:num_offspring_community # enumerate the household id of each of the community offspring (each into a new household)

    ## Looping through community offspring and sampling details of their occupation and household that they belong to, and their occupation
    for (i in 1:num_offspring_community) {

      ## Sampling a household for community offspring to be produced into, based on the age group of the index infection
      offspring_community_age_group <- offspring_community_ages[i]
      col_name <- age_group_to_column[offspring_community_age_group]
      sampled_row <- synthetic_household_df[get(col_name) == 1][sample(.N, 1)]

      ## Picking the particular household member infected
      possible_hh_member_indices <- which(unlist(sampled_row$hh_ages_grouped) %in% offspring_community_age_group)
      if (length(possible_hh_member_indices) == 1) {
        hh_member_index <- possible_hh_member_indices
      } else {
        hh_member_index <- sample(possible_hh_member_indices, 1) # sampling which individual with that occupation in the hh is the infection
      }

      ## Creating dataframe with information on this particular offspring
      offspring_community_list[[i]] <- data.table::data.table(transmission_route = "community",
                                                              occupation = unlist(sampled_row$hh_occupations)[hh_member_index],
                                                              age = offspring_community_age_group,
                                                              hh_id = offspring_community_new_hh_id[i],
                                                              hh_member_index = hh_member_index,
                                                              hh_size = sampled_row$hh_size,
                                                              hh_ages = sampled_row$hh_ages,
                                                              hh_occupations = sampled_row$hh_occupations,
                                                              hh_infections = 1,
                                                              hh_infected_index = hh_member_index)
    }
    offspring_characteristics_df_community <- data.table::rbindlist(offspring_community_list, fill = TRUE)

  ## If no community offspring generated, initialise an empty community offspring df
  } else {
    offspring_characteristics_df_community <- data.frame(transmission_route = character(0),
                                                         occupation = character(0),
                                                         age = character(0),
                                                         hh_id = numeric(0),
                                                         hh_member_index = numeric(0),
                                                         hh_size = numeric(0),
                                                         hh_ages = numeric(0),
                                                         hh_occupations = numeric(0),
                                                         hh_infections = numeric(0),
                                                         hh_infected_index = numeric(0))
  }

  ##########################################################################################################################################

  return(list(total_offspring = num_offspring_sexual + num_offspring_hh + num_offspring_community,
              num_offspring_sexual = num_offspring_sexual,
              num_offspring_hh = num_offspring_hh,
              new_hh_infected_index = list(cumulative_hh_infected_index),  ## make sure to use this to update the number of infections left in the household for the index infection!!!!
              num_offspring_community = num_offspring_community,
              offspring_characteristics = rbind(offspring_characteristics_df_sexual, offspring_characteristics_df_hh, offspring_characteristics_df_community)))
}

# load("data/south_kivu_final.rda")
# synthetic_household_df <- south_kivu_final
# household_df <- south_kivu_final[8946, ]
# index_age_group <- "18+"
# index_hh_id <- 8946
# index_occupation <- "SW"
# max_hh_id <- 20
# mn_offspring_sexual_SW <- 5
# disp_offspring_sexual_SW <- 100
# mn_offspring_sexual_PBS <- 5
# disp_offspring_sexual_PBS <- 10
# mn_offspring_sexual_genPop <- 5
# disp_offspring_sexual_genPop <- 10
# sexual_transmission_occupation_matrix <- matrix(data = c(0.00, 0.10, 0.90,
#                                                          0.20, 0.10, 0.70,
#                                                          0.00, 0.90, 0.10), nrow = 3, ncol = 3, byrow = TRUE)
# mn_offspring_hh <- 7
# disp_offspring_hh <- 100
# index_hh_member_index <- which(unlist(household_df$hh_occupations) == "SW")
# index_hh_infections_index <- c(1, 4, 14)
# index_hh_size <- household_df$hh_size
# index_hh_prior_infections <- length(index_hh_infections_index)
# index_hh_ages <- unlist(household_df$hh_ages_grouped)
# index_hh_occupations <- unlist(household_df$hh_occupations)
# mn_offspring_community <- 5
# disp_offspring_community <- 100
# number_susceptible_community <- 10^5
# population_community <- 10^5
# community_transmission_age_matrix <- matrix(data = c(0.34, 0.33, 0.33,
#                                                      0.34, 0.33, 0.33,
#                                                      0.34, 0.33, 0.33), nrow = 3, ncol = 3, byrow = TRUE)

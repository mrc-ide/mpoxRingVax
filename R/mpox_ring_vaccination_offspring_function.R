#' @export
offspring_fun <- function(index_age_group,                           # age of the index infection (0-5, 5-18, 18+)
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
                          index_hh_size,                             # size of index infection's household (excluding the index infection)
                          index_hh_prior_infections,                 # number of prior infections in the index's household (excluding the index infection)
                          index_hh_ages,                             # ages of index infection's household (excluding the index infection and others already infected)
                          index_hh_occupations,                      # occupations of the index infection's household (excluding the index infection and others already infected)

                          mn_offspring_community,                    # mean of the community offspring distribution
                          disp_offspring_community,                  # overdispersion of the community offspring distribution
                          number_susceptible_community,              # number of susceptible individuals remaining in the community
                          population_community,                      # overall population size
                          community_transmission_age_matrix,         # rows are the proportions of cases in each age-group a particular age-group gives rises to (note this needs to take into account smallpox vaccination)
                          community_transmission_occupation_matrix   # rows are the proportions of cases in each age-group a particular age-group gives rises to
) {

  #########################################################################################################################
  ## Generating infections from sexual transmission and their characteristics
  ## - Note that we assume:
  ##     1) All SW and PBS are >18
  ##     2) No sexual transmission to/from <18
  ##     3) Single mpox introduction per household (thus all infections are in newly instantiated and separate households)
  #########################################################################################################################
  if (index_age_group %in% c("0-5", "5-18")) {
    num_offspring_sexual <- 0
    offspring_sexual_new_hh_id <- max_hh_id
    offspring_characteristics_df_sexual <- data.frame(transmission_route = character(0),
                                                      occupation = character(0),
                                                      age = character(0),
                                                      hh_id = numeric(0))
  } else {
    if (index_occupation == "SW") {
      num_offspring_sexual <- rnbinom(n = 1, mu = mn_offspring_sexual_SW, size = disp_offspring_sexual_SW)
    } else if (index_occupation == "PBS") {
      num_offspring_sexual <- rnbinom(n = 1, mu = mn_offspring_sexual_PBS, size = disp_offspring_sexual_PBS)
    } else if (index_occupation == "genPop") {
      num_offspring_sexual <- rnbinom(n = 1, mu = mn_offspring_sexual_genPop, size = disp_offspring_sexual_genPop)
    } else {
      stop("something has gone wrong with the sexual offspring distribution generation of infections")
    }

    if (num_offspring_sexual > 0) {
      offspring_sexual_occupation <- sample(x = c("SW", "PBS", "genPop"), size = num_offspring_sexual, prob = sexual_transmission_occupation_matrix[index_occupation, ], replace = TRUE)
      offspring_sexual_ages <- rep("18+", num_offspring_sexual)
      offspring_sexual_new_hh_id <- max_hh_id + 1:num_offspring_sexual
      offspring_characteristics_df_sexual <- data.frame(transmission_route = "sexual",
                                                        occupation = offspring_sexual_occupation,
                                                        age = offspring_sexual_ages,
                                                        hh_id = offspring_sexual_new_hh_id)
    } else {
      offspring_characteristics_df_sexual <- data.frame(transmission_route = character(0),
                                                        occupation = character(0),
                                                        age = character(0),
                                                        hh_id = numeric(0))
    }
  }

  ##########################################################################################################################################
  ## Generating infections from household transmission and their characteristics
  ## - Note that we aassume:
  ##     1) Frequency-dependent transmission a la https://elifesciences.org/articles/70767 but might want to change this
  ##########################################################################################################################################
  number_susceptible_hh <- index_hh_size - index_hh_prior_infections
  new_mn_hh <- mn_offspring_hh * number_susceptible_hh/index_hh_size
  new_disp_hh <- new_mn_hh/(disp_offspring_hh - 1)
  num_offspring_hh <- min(c(rnbinom(n = 1, mu = new_mn_hh, size = new_disp_hh), number_susceptible_hh))
  if (num_offspring_hh > 0) {
    infected_hh_index <- sample(x = 1:number_susceptible_hh, size = num_offspring_hh, replace = FALSE)
    offspring_hh_ages <- index_hh_ages[infected_hh_index]
    offspring_hh_occupation <- index_hh_occupations[infected_hh_index]
    offspring_hh_id <- rep(index_hh_id, num_offspring_hh)
    offspring_characteristics_df_hh <- data.frame(transmission_route = "household",
                                                  occupation = offspring_hh_occupation,
                                                  age = offspring_hh_ages,
                                                  hh_id = offspring_hh_id)
  } else {
    offspring_characteristics_df_hh <- data.frame(transmission_route = character(0),
                                                  occupation = character(0),
                                                  age = character(0),
                                                  hh_id = numeric(0))
  }

  ##########################################################################################################################################
  ## Generating infections from community transmission and their characteristics
  ## - Note that we assume:
  ##     1) Single mpox introduction per household (thus all infections are in newly instantiated and separate households)
  ##########################################################################################################################################
  new_mn_community <- mn_offspring_community * number_susceptible_community/population_community
  new_disp_community <- new_mn_community/(disp_offspring_community - 1)
  num_offspring_community <- rnbinom(n = 1, mu = new_mn_community, size = new_disp_community)
  if (num_offspring_community > 0) {
    offspring_community_occupation <- sample(c("SW", "PBS", "genPop"), size = num_offspring_community, prob = community_transmission_occupation_matrix[index_occupation, ], replace = TRUE)
    offspring_community_ages <- sample(c("0-5", "5-18", "18+"), size = num_offspring_community, prob = community_transmission_age_matrix[index_age_group, ], replace = TRUE)
    offspring_community_new_hh_id <- max(offspring_sexual_new_hh_id) + 1:num_offspring_community
    offspring_characteristics_df_community <- data.frame(transmission_route = "community",
                                                         occupation = offspring_community_occupation,
                                                         age = offspring_community_ages,
                                                         hh_id = offspring_community_new_hh_id)
  } else {
    offspring_characteristics_df_community <- data.frame(transmission_route = character(0),
                                                         occupation = character(0),
                                                         age = character(0),
                                                         hh_id = numeric(0))
  }


  return(list(total_offspring = num_offspring_sexual + num_offspring_hh + num_offspring_community,
              num_offspring_sexual = num_offspring_sexual,
              num_offspring_hh = num_offspring_hh,
              num_offspring_community = num_offspring_community,
              offspring_characteristics = rbind(offspring_characteristics_df_sexual, offspring_characteristics_df_hh, offspring_characteristics_df_community)))
}

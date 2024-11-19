#' @export

load("data/south_kivu_final.rda")

household <- south_kivu_final[south_kivu_final$hh_id == 506, ]

index_age_group <- 3
index_hh_id <- 10
index_occupation <- 1
max_hh_id <- 12

mn_offspring_sexual_SW <- 8
disp_offspring_sexual_SW <- 1.1
mn_offspring_sexual_PBS <- 4
disp_offspring_sexual_PBS <- 1.1
mn_offspring_sexual_genPop <- 1
disp_offspring_sexual_genPop <- 1.1
sexual_transmission_occupation_matrix <- matrix(data = c(0.00, 0.90, 0.10,
                                                         0.34, 0.33, 0.33,
                                                         0.00, 0.10, 0.90), nrow = 3, ncol = 3, byrow = TRUE)
mn_offspring_hh <- 3
disp_offspring_hh <- 1.1
index_hh_size <- household$hh_size
index_hh_prior_infections <- 1
index_hh_occupations <- household$hh_occupations

index_hh_ages <- unlist(household$hh_ages_grouped)

mn_offspring_community <- 2
disp_offspring_community <- 1.1
number_susceptible_community <- 10^5
population_community <- 10^5
community_transmission_age_matrix <- matrix(data = c(0.34, 0.33, 0.33,
                                                     0.34, 0.33, 0.33,
                                                     0.34, 0.33, 0.33), nrow = 3, ncol = 3, byrow = TRUE)

## Use community_transmission_age_matrix to figure out the ages of the individuals infected, and then select a household that as at least one of those people

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

                          drc_hh_age_occ_df                          # dataframe of synthetic drc household/age/occ population
) {

  ## add in "household member number" to track index of occupations and age that the person is in the vector in the list
  ## change it so that it doesn't exclude the index infection

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
      offspring_sexual_occupation <- sample(x = c("SW", "PBS", "genPop"), size = num_offspring_sexual,
                                            prob = sexual_transmission_occupation_matrix[index_occupation, ], replace = TRUE)
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
  ##     1) Frequency-dependent transmission e.g. https://elifesciences.org/articles/70767 but might want to change this
  ##########################################################################################################################################
  number_susceptible_hh <- index_hh_size - index_hh_prior_infections ## NOTE THAT INDEX_HH_SIZE MUST EXCLUDE THE INDEX CASE
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
  ## COMMUNITY INFECTIONS: Generating infections from community transmission and their characteristics
  ## - Note that we assume:
  ##     1) Single mpox introduction per household (thus all infections are in newly instantiated and separate households)
  ##########################################################################################################################################
  new_mn_community <- mn_offspring_community * number_susceptible_community/population_community        # mean of community offspring distribution taking susceptible depletion into account
  new_disp_community <- new_mn_community/(disp_offspring_community - 1)
  num_offspring_community <- rnbinom(n = 1, mu = new_mn_community, size = new_disp_community)
  age_group_to_column <- c("0-5" = "contains_0_5", "5-18" = "contains_5_18", "18+" = "contains_18plus")

  ## Populating the df if there are community offspring
  if (num_offspring_community > 0) {

    offspring_list <- vector("list", num_offspring_community)
    offspring_community_ages <- sample(c("0-5", "5-18", "18+"), size = num_offspring_community, prob = community_transmission_age_matrix[index_age_group, ], replace = TRUE)
    offspring_community_new_hh_id <- max(offspring_sexual_new_hh_id) + 1:num_offspring_community

    for (i in 1:num_offspring_community) {
      offspring_community_age_group <- offspring_community_ages[i]
      col_name <- age_group_to_column[offspring_community_age_group]
      sampled_row <- south_kivu_final[get(col_name) == 1][sample(.N, 1)]
      hh_member_index <- sample(which(unlist(sampled_row$hh_ages_grouped) %in% offspring_community_age_group), 1)
      offspring_list[[i]] <- data.table::data.table(transmission_route = "community",
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

    offspring_characteristics_df_community <- data.table::rbindlist(offspring_list, fill = TRUE)

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
              num_offspring_community = num_offspring_community,
              offspring_characteristics = rbind(offspring_characteristics_df_sexual, offspring_characteristics_df_hh, offspring_characteristics_df_community)))
}

# Loading most recent DRC DHS survey
library(rdhs); library(haven); library(ggplot2); library(dplyr)
set_rdhs_config(email = "charles.whittaker16@imperial.ac.uk",
                project = "mpox Outbreak Response - Understanding Household Size Distributions in DRC")
survs <- dhs_surveys(countryIds = c("CD"), surveyYear = 2013)
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileType = "HR", fileFormat="flat")
datasets$path <- unlist(get_datasets(datasets$FileName))

# Extracting out household size distribution
x <- readRDS(datasets$path)
labels <- colnames(x)
hh_variable_labels <- search_variable_labels(datasets$FileName, "*hold")[, 1:2]
weight_variable_labels <- search_variable_labels(datasets$FileName, "weight")[, 1:2]
y <- x[, c("hv001", "hv004", "hv024",  "hv002", "hv005", "hv009")]
hist(y$hv009, breaks = 20)
z <- y[sample(x = 1:nrow(y), size = nrow(y), replace = TRUE, prob = y$hv005 / 10^6), ]
hist(z$hv009, breaks = 30)

# Getting the ages of individuals in each household
hh_age_labels <- c(paste0("hv105_0", 1:9), paste0("hv105_1", 0:9), paste0("hv105_2", 0:4))
for (i in 1:nrow(y)) {
  temp_hh_ages <- x[i, hh_age_labels]
  temp_hh_ages <- temp_hh_ages[which(!is.na(temp_hh_ages))]
  y$hh_size_from_ages[i] <- length(temp_hh_ages)
  y$hh_ages[i]  <- list(unname(temp_hh_ages))
  if (i %% 1000 == 0) {
    print(i)
  }
}
sum(y$hv009 != y$hh_size_from_ages)
y$province <- haven::as_factor(y$hv024)
ggplot(y, aes(x = hv009, fill = province)) +
  geom_histogram(bins = 20, col = "black") +
  facet_wrap(. ~ province, scales = "free_y")

# Selecting South Kivu specifically and then populating with SW/PBS/genPop estimates
south_kivu <- y %>%
  filter(province == "sud-kivu")
south_kivu <- rbind(south_kivu, south_kivu, south_kivu, south_kivu, south_kivu,
                    south_kivu, south_kivu, south_kivu, south_kivu, south_kivu)
prop_sud_kivu_under20 <- sum(unlist(south_kivu$hh_ages) >= 20) / length(unlist(south_kivu$hh_ages))
prop_sud_kivu_SW_over20 <- 0.007 * 0.5 * 1 / (1 - prop_sud_kivu_under20)
prop_sud_kivu_PBS_over20 <- 0.11 * 0.5 * 1 / (1 - prop_sud_kivu_under20)

# Creating a column for whether it contains over 20s (hence eligible to have SW/PBS)
for (i in 1:nrow(south_kivu)) {
  temp_ages <- unlist(south_kivu$hh_ages[i][[1]])
  over_20 <- temp_ages >= 20
  if (sum(over_20) > 0) {
    south_kivu$eligible[i] <- "Yes"
    south_kivu$num_eligible[i] <- sum(over_20)
  } else {
    south_kivu$eligible[i] <- "No"
    south_kivu$num_eligible[i] <- 0
  }
  south_kivu$assigned[i] <- "No"
}

# Getting ages into groupings, and creating hh occupation column to modify later (initially all set to genPop)
for (i in 1:nrow(south_kivu)) {
  temp_ages <- unlist(south_kivu$hh_ages[i][[1]])
  south_kivu$hh_occupations[i] <- list(rep("genPop", length(temp_ages)))
  hh_ages_grouped <- ifelse(temp_ages <= 5, "0-5", ifelse(temp_ages <= 18, "5-18", "18+"))
  south_kivu$hh_ages_grouped[i] <- list(factor(hh_ages_grouped, levels = c("0-5", "5-18", "18+")))
}

# Populating the synthetic households with SW and PBS
## Note - we assume households can only have 1 of (SW/PBS) at most
ages <- unlist(south_kivu$hh_ages)
num_people_over20 <- sum(ages >= 20)
num_SW <- round(num_people_over20 * prop_sud_kivu_SW_over20)
num_PBS <- round(num_people_over20 * prop_sud_kivu_PBS_over20)

indices_eligible <- which(south_kivu$eligible == "Yes")
indices_assigned <- sample(x = indices_eligible, size = num_SW + num_PBS, replace = FALSE, prob = south_kivu$num_eligible[indices_eligible])
indices_assigned_SW <- indices_assigned[1:num_SW]
indices_assignd_PBS <- indices_assigned[(num_SW + 1):(num_SW + num_PBS)]

## Assigning the SW to their households
for (i in 1:length(indices_assigned_SW)) {
  index <- indices_assigned_SW[i]
  temp_ages <- unlist(south_kivu$hh_ages[index][[1]])
  temp_occ <- rep("genPop", length(temp_ages))
  over_20 <- which(temp_ages >= 20)
  indiv_occ_index <- over_20[sample(1:length(over_20), 1)]
  temp_occ[indiv_occ_index] <- "SW"
  south_kivu$assigned[index] <- "Yes"
  south_kivu$hh_occupations[index] <- list(unname(temp_occ))
}

## Assigning the PBS to their households
for (i in 1:length(indices_assignd_PBS)) {
  index <- indices_assignd_PBS[i]
  temp_ages <- unlist(south_kivu$hh_ages[index][[1]])
  temp_occ <- rep("genPop", length(temp_ages))
  over_20 <- which(temp_ages >= 20)
  indiv_occ_index <- over_20[sample(1:length(over_20), 1)]
  temp_occ[indiv_occ_index] <- "PBS"
  south_kivu$assigned[index] <- "Yes"
  south_kivu$hh_occupations[index] <- list(unname(temp_occ))
}

occupations <- unlist(south_kivu$hh_occupations)
ages <- unlist(south_kivu$hh_ages)
age_groups <- unlist(south_kivu$hh_ages_grouped)

table(occupations)
table(ages)
table(age_groups)
table(age_groups, occupations)

colnames(south_kivu)
south_kivu$hh_id <- 1:nrow(south_kivu)
south_kivu_final <- south_kivu %>%
  select(hh_id, hh_size_from_ages, hh_ages_grouped, hh_occupations) %>%
  rename(hh_size = hh_size_from_ages)

south_kivu_final <- data.table::as.data.table(south_kivu_final)
usethis::use_data_raw("south_kivu_final")

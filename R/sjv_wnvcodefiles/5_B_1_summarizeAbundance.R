
# Packages ----------------------------------------------------------------

library(sf)
library(tidyverse)
library(tmap)
library(lubridate)
library(scales)
library(PooledInfRate)


tmap_mode("view")

# Function to read in column data types correctly
source("R/wnvcodefiles/1_0_dataRead.R")


dir.create("Data/SJV/wnv/5_summarizeClusters/")

# Read in Data ------------------------------------------------------------


# id : Cluster key
clustKey <- read_csv("Data/SJV/wnv/3_clusterExport/idClusterKey.csv")

# Abundance data
abund <- readAbund("Data/SJV/wnv/2_spatialExploration/abundID.csv") %>%
  inner_join(clustKey, by = "id") %>%
  mutate(
    collection_date = mdy(collection_date) %>% floor_date(unit = "week")
  )

# Pool data
pools <- readPools("Data/SJV/wnv/2_spatialExploration/poolsID.csv") %>%
  inner_join(clustKey, by = "id")  %>%
  mutate(
    collection_date = mdy(collection_date) %>% floor_date(unit = "week")
  )

# WNV testing data
test <- readTest("Data/SJV/wnv/2_spatialExploration/testID.csv") %>%
  inner_join(clustKey, by = "id")   %>%
  mutate(
    collection_date = mdy(collection_date) #%>% floor_date(unit = "week")
  )





# select Columns of Interest ----------------------------------------------

# This creates a vector of all columns across the three datasets. It is then
# subsetted so that extraneous columns are removed.

abundNames <- colnames(abund)
poolsNames <- colnames(pools)
testNames <- colnames(test)


combineNames <- c(abundNames, poolsNames, testNames) %>%
  unique()

# Removes extraneous columns (calculated_county, zip, etc.)
# combineNames <- combineNames[28:140]



# Abundance and Testing Cleaning ------------------------------------------



# Format abundance data for easier reading
abundFilter <- abund %>%
  # Removes location info columns
  dplyr::select(any_of(combineNames)) %>%
  # Filters to species of interest
  # dplyr::select(contains(c("culex_pipiens", "culex_tarsalis", "culex_quinquefasciatus")))  %>%
  dplyr::select(-matches("culex_((erythrothorax)|(stigmatosoma)|(territans)|(thriambus)).*"))  %>%
  # Trap problem is a column in the raw data where numbers might not be accurate
  filter(trap_problem == "N") %>%
  # Format of easier reading
  dplyr::select(id, collection_date, everything()) %>%
  dplyr::select(-ends_with("gravid")) %>%
  filter(trap_type != "GRVD")





# Summarize abundance data ------------------------------------------------



# Group by cluster / week and plot different species abundance
## wnv vectors

chooseSpecies <- abundFilter %>%
  pivot_longer(contains(c("culex", "aedes_aegypti")), names_to = "speciesLife") %>%
  mutate(
    species = case_when(
      str_detect(speciesLife, "restuans") ~ "pipiensMix",
      str_starts(speciesLife, "culex_pipiens") ~ "pipiens",
      str_starts(speciesLife, "culex_tarsalis") ~ "tarsalis",
      str_starts(speciesLife, "culex_quinquefasciatus") ~ "quinquefasciatus",
      str_starts(speciesLife, "aedes_aegypti") ~ "aegypti",
      T ~ "other")
  ) %>%
  mutate(
    lifeCycle = case_when(str_detect(speciesLife, "females") ~ "female",
                          str_detect(speciesLife, "males") ~ "male",
                          str_detect(speciesLife, "(eggs)|(larvae)|(pupae)") ~ "juvenile",
                          T ~ "unknown")
  ) %>%
  mutate(
    value = replace_na(value, 0)
  )


unique(pull(chooseSpecies, lifeCycle))

# see if getting rid of different species would help

cs2 <- chooseSpecies %>%
  mutate(
    woy = week(collection_date),
    year = year(collection_date) %>% as.factor(),
    trap_nights = if_else(trap_nights == 0, 1, trap_nights),
    traps = num_trap * trap_nights
  ) %>%
  group_by(clust, woy, year, species, lifeCycle, collection_date) %>%
  summarize(
    points = n(),
    trapNights = sum(traps),
    totalMosquitoes = sum(value)
  ) %>%
  # Final adjustment
  mutate(
    # Adjust mosquito numbers for number of traps / trap nights
    mosPerTrapNight = (totalMosquitoes / trapNights) %>% round(., 2)
  ) %>%
  unite("wy", year, woy, remove = FALSE)

cs2 %>%
  filter(totalMosquitoes > 0)


#Write out

write_csv(cs2, "Data/SJV/wnv/5_summarizeClusters/abundSpeciesSummary.csv")



glimpse(cs2)

# Filter out non females


cs2 %>%
  filter(lifeCycle == "female")



# Check high mosquito numbers ---------------------------------------------


abundClust <- abundFilter %>%
  ungroup() %>%
  # Sum mosquito numbers across columns
  mutate(
    # Set 0 trap night value to 1
    trap_nights = if_else(trap_nights == 0, 1, trap_nights),
    mos = rowSums(dplyr::select(., starts_with("culex")), na.rm = TRUE),
    traps = num_trap * trap_nights
  ) %>%
  # Get number of trapping nights, mosquitoes by cluster
  group_by(clust, collection_date) %>%
  summarise(
    points = n(),
    trapNights = sum(traps),
    totalMosquitoes = sum(mos)
  ) %>%
  # Final adjustment
  mutate(
    # Adjust mosquito numbers for number of traps / trap nights
    mosPerTrapNight = (totalMosquitoes / trapNights) %>% round(., 2),

    # Get some date metrics
    woy = week(collection_date),
    year = year(collection_date) %>% as.factor()
  ) %>%
  unite("wy", year, woy, remove = FALSE)


abundClust %>%
  filter(
    totalMosquitoes > 4000
  )



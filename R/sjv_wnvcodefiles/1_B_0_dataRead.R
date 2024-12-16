# packages
library(dplyr)

# read in old data for reference to match with the script

old_abundance <- read.csv("Data/0_input/MosquitoAbundance2010.2020.csv")
old_pools <- read.csv("Data/0_input/MosquitoPools2010.2020.csv")
old_testing <- read.csv("Data/0_input/MosquitoPoolsTestResults2010.2020.csv")

# read in new data
KMVCD_abundance <- read.csv("Data/0_input/KMVCD_Mosquito_Collections_2010_2023.csv") %>%
  mutate(collection_date = mdy(collection_date) %>% as.character())

write.csv(KMVCD_abundance, "Data/KMVCD/wnv/Raw/Mosquito_Abundance_2010_2023.csv")

KMVCD_pools <- read.csv("Data/0_input/KMVCD_Mosquito_Pools_2010_2023.csv") %>%
  filter(test_target == "WNV")

write.csv(KMVCD_pools, "Data/KMVCD/wnv/Raw/Mosquito_Pools_2010_2023.csv")

# create KMVCD_testing
KMVCD_testing <- KMVCD_abundance %>%
  right_join(KMVCD_pools, by = intersect(names(KMVCD_abundance), names(KMVCD_pools))) %>%
  dplyr::select(any_of(names(old_testing)))

write.csv(KMVCD_testing, "Data/KMVCD/wnv/Raw/Mosquito_Testing_2010_2023.csv")

KMVCD_testing <- read.csv("Data/KMVCD/wnv/Raw/Mosquito_Testing_2010_2023.csv")




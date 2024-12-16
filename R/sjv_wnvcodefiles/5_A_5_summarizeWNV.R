##############################################################################
###############################################################################

# FLOORING DATE AND CLUSTERING BEFORE RUNNING PIR

##############################################################################
###############################################################################






# Packages ----------------------------------------------------------------
# xfun::install_github("https://github.com/CDCgov/PooledInfRate")


library(sf)
library(tidyverse)
library(tmap)
library(lubridate)
library(scales)
library(PooledInfRate)


tmap_mode("view")

# Function to read in column data types correctly
source("R/sjv_wnvcodefiles/1_0_dataRead.R")


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

# # WNV testing data
test <- readTest("Data/SJV/wnv/2_spatialExploration/testID.csv") %>%
  inner_join(clustKey, by = "id")   %>%
  mutate(
    collection_date = mdy(collection_date) %>% floor_date(unit = "week")
  )





# Select Columns of Interest ----------------------------------------------

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

# Format testing data for easier reading
testCols <- test %>%
  dplyr::select(any_of(combineNames)) %>%
  dplyr::select(id, collection_date, everything())





# Summarize WNV Testing Data ----------------------------------------------


# Format Data according to PIR package, grouping for each ID date, and species


testSpecies <- testCols %>%
  arrange(clust, desc(collection_date)) %>%

  # Remove non female pools
  filter(sex_condition != "Unknown Sex") %>%

  # For every date, location, species summarize
  group_by(clust, collection_date, species, test_status) %>%
  summarise(
    traps = n(),
    count = sum(num_count %>% as.integer())
  ) %>%
  summarize(
    # Number of Mosquitoes across pools
    PoolSize = sum(count),
    # Number of positive pools
    Pos = sum(traps[test_status == "Confirmed"]) %>%
      as.numeric(),
    # Total Pools
    NumPools = sum(traps) %>%
      as.numeric()
  ) %>%
  ungroup()




# See how much of the data can be used in PIR -----------------------------


# When there are no negative Pools in a group, this is an issue.
testSpecies2Raw <- testSpecies %>%
  mutate(
    woy = week(collection_date),
    year = year(collection_date),
    month = month(collection_date)
  ) #%>%
# Filtering Pos < NumPools removes first set of problems
# filter(
#   Pos >= NumPools
# ) %>%
# distinct(clust, woy, year)




Pos <- testSpecies %>%
  filter(Pos > 0)


PosLeft <- Pos %>%
  filter(Pos < NumPools)


# Split table into list of 1 row tables for every observation
# This is to test whether or not the PIR function will work
testSpecies2 <- testSpecies2Raw %>%
  ungroup() %>%
  mutate(rowNum = row_number()) %>%
  group_split(rowNum)



# PIR function
mapPIR2 <- function(table) {


  newTable <- table


  estimation <-  pooledBin(Pos ~ m(PoolSize) + n(NumPools), data = newTable)


  toJoin <- estimation %>%
    as.data.frame() %>%
    as_tibble()

  out <- table %>%
    bind_cols(toJoin)



  out
}


# Returns NULL if PIR function works on row, and rowNumber if it throws error
PIRtest <- function(data) {

  out <- tryCatch(
    expr = {
      # Your code...
      # goes here...
      # ...
      mapPIR2(data)
      NULL
    },
    error = function(e){
      # (Optional)
      # Do this if an error is caught...
      return(as.numeric(data[[1, "rowNum"]]))
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>'
      # message(paste("Processed URL:", url))
      # message("Some other message at the end")
    }
  )

  return(out)


}


# Create vector of all errors
errors <- map(testSpecies2, PIRtest)

# Remove NULL values
errors2 <- errors %>%
  compact() %>%
  simplify()


# Remove error rowNumbers from dataset
testSpecies3 <- testSpecies2 %>%
  bind_rows() %>%
  filter(!(rowNum %in% errors2)) %>%
  group_split(rowNum)

# Run PIR on cleaned dataset
mappedPIR <- map(testSpecies3, mapPIR2) %>%
  bind_rows()



mappedPIR %>%
  filter(is.na(P))


# Subset data that PIR won't work for
testSpecies4 <- testSpecies2 %>%
  bind_rows() %>%
  filter(rowNum %in% errors2)

# Bind Rows
MIR <- mappedPIR %>%
  bind_rows(testSpecies4) %>%
  mutate(
    MIR = (Pos / PoolSize) * 1000,
    posTestFreq = Pos / NumPools
  )



# Write out
write_csv(MIR, "Data/SJV/wnv/5_summarizeClusters/testPIRMIRClusterTEST.csv")









# Appendix ----------------------------------------------------------------
#
# testRaw <- read_csv("Data/0_input/MosquitoPoolsTestResults2010.2020.csv")
#
# glimpse(testRaw)
#
#
# testRaw[,c(17, 9, 15, 18, 32)]
#
#
# test %>%
#   filter(collection_date == "2012-07-23") %>%
#   filter(id == 1260) %>%
#   view()
#
# testRaw %>%
#   filter(collection_date == "7/23/2012") %>%
#   view()
# #
# #
# # x <- c(1)
# #
# # p <- c(4)
# #
# # pooledBin(x, p)
# # pooledBin(1, 4)
# #
# #
# # testSpecies2[[3667]] %>% glimpse()
# #
# # testSpecies2[[3667]] %>% mapPIR2()
# #
# # tryTestMap <- map(testSpecies2, mapPIR2)
# #
#
#
# print(testSpecies2[[2385]])
#
# mapPIR2(testSpecies2[[2385]])
#
#
# print(removeRow, n = 24)
# mapPIR2(removeRow)
#
# mapPIR2(removeRow[20,])
#
#
# print(testSpecies2[[2385]])
#
# mapPIR2(testSpecies2[[2385]])
#
#
# print(removeRow, n = 24)
# mapPIR2(removeRow)
#
# mapPIR2(removeRow[20,])
#
#
# Positives <- testSpecies %>%
#   group_by(clust, collection_date) %>%
#   filter(
#     Pos > 0
#   ) %>%
#   select(clust, collection_date)
#
#
# # posExpand <- testSpecies %>%
# #   inner_join(Positives, by = c("clust", "collection_date"))
# #
# # toDrop <- posExpand %>%
# #   group_by(clust, collection_date) %>%
# #   mutate(
# #     flag = if_else(Pos == 0, 1, 0)
# #   ) %>%
# #   summarize(
# #     flag = sum(flag)
# #   ) %>%
# #   filter(flag < 1)
#
# #
# # pirPre <- testSpecies
# #   group_split(clust, year, species)
# #
# # #
# #
# # library(PooledInfRate)
# #
# # testCols %>%
# #   distinct(agency_pool_num)
# #
# # testCols %>%
# #   distinct(id)
# #
# #
# # testCols %>%
# #   group_by(id, collection_date, species) %>%
# #   mutate(
# #     n = n()
# #   ) %>%
# #   filter(n > 1)
# #
# #
# #
#
#
#
# PIRtest(testSpecies2[[631]])
#
#
# errorTibble %>%
#   add_row(testSpecies2[[9132]] %>% as.data.frame())
#
# testSpecies2[[54]]
#
# removeRow <- testSpecies2[[2385]] %>%
#   mutate(
#
#     # # Reduce number of pools / positive pools
#     # Pos = ceiling(Pos / 4),
#     # NumPools = ceiling(NumPools / 2)
#
#     # # Check Pos / Negative
#     # Pos = if_else(Pos > 0, 1, 0),
#     # NumPools = if_else(NumPools > 9, 9, NumPools)
#
#     PoolSize = if_else(PoolSize > 1048, 1048, PoolSize)
#
#   )

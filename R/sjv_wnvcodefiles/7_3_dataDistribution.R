
# Packages ----------------------------------------------------------------




library(sf)
library(tidyverse)
library(tmap)
library(lubridate)
library(scales)
library(hrbrthemes)


tmap_mode("view")

dir.create("Data/7_dataDistribution/")

# Read in Data ------------------------------------------------------------


abund <- read_csv("Data/5_summarizeClusters/abundSpeciesSummary.csv")
abundInt <- read_csv("Data/6_interpolateData/abundClusterInterpolation.csv")
abundExt <- read_csv("Data/6_interpolateData/abundClusterInterpolationInferNA.csv")

test <- read_csv("Data/5_summarizeClusters/testPIRMIRClusterTEST.csv")
testInt <- read_csv("Data/6_interpolateData/testClusterInterpolation.csv")
testExt <- read_csv("Data/6_interpolateData/testClusterInterpolationInferNA.csv")


dataList <- list(abund, abundInt, abundExt, test, testInt, testExt)



map(dataList, glimpse)


# Filter Data to 2019 + 2020 ----------------------------------------------




# Filter to 2019 and 2020 (ebird Data years)
yearFilter <- function(tibble) {
  tibble %>%
    filter((year == 2019) | (year == 2020))
}


# Filter Abundance and Testing to 2019-2020
dataYearsFilter <- map(dataList, yearFilter)





# Join Abundance and Testing Data -----------------------------------------



# Clean 2 datasets for easier use
abundCombo <- dataList[[1]] %>%
  full_join(dataList[[2]], by = c('clust', 'wy', 'woy', 'year', 'species', 'lifeCycle', 'points', 'trapNights', 'totalMosquitoes', 'mosPerTrapNight', 'collection_date')) %>%
  full_join(dataList[[3]], by = c('clust', 'wy', 'woy', 'year', 'species', 'lifeCycle', 'points', 'trapNights', 'totalMosquitoes', 'mosPerTrapNight', 'collection_date')) %>%
  dplyr::select(clust, woy, year, collection_date) %>%
  mutate(data = "abund")

testCombo <- dataList[[4]] %>%
  full_join(dataList[[5]], by = c('clust', 'collection_date', 'species', 'PoolSize', 'Pos', 'NumPools', 'woy', 'year', 'month', 'rowNum', 'P', 'Lower', 'Upper', 'MIR', 'posTestFreq')) %>%
  full_join(dataList[[6]], by = c('clust', 'collection_date', 'species', 'PoolSize', 'Pos', 'NumPools', 'woy', 'year', 'month', 'rowNum', 'P', 'Lower', 'Upper', 'MIR', 'posTestFreq')) %>%
  dplyr::select(clust, woy, year, collection_date) %>%
  mutate(data = "test")


dataColTiny <- list(abund = abundCombo, test = testCombo)

# bind_rows Join between the data


# This is for plots of both datasets to show difference in the amount of data
# between them, as well as when data was most frequently collected
abundTestCombo <- abundCombo %>%
  bind_rows(testCombo) %>%
  # Coerce "data" column as factor for plot visualization
  mutate(data = as.factor(data),
         yr = as.character(year),
         dt = as.character(data)) %>%
  # Create year-date signifier, just in case
  unite("yearData", dt, yr)



write_csv(abundTestCombo, "Data/7_dataDistribution/abundanceTestingClean.csv")



# Plot Temporal Distribution ----------------------------------------------



# Abundance Data histograms + Testing Data Histograms over 2 year period
abundTestCombo %>%
  ggplot(mapping = aes(x = collection_date, fill = data)) +
  geom_histogram(binwidth = 16,
                 color = "white",
                 alpha = 0.5,
                 position = 'dodge'
  ) +
  labs(x = "Collection Date", y = "Number of Samples",
       title = "2019-2020 Abundance & Testing Data Temporal Distribution") +
  #scale_y_continuous(breaks = seq(0,2000, 200)) +
  scale_x_date(breaks = date_breaks("2 month"),
               labels = date_format("%m/%y") )


# Abundance vs Testing Distribution by Week of Year

abundTestCombo %>%
  ggplot(mapping = aes(x = woy, fill = data)) +
  geom_histogram(binwidth = 1,
                 color = "white",
                 alpha = 0.5,
                 position = 'dodge'
  ) +
  labs(x = "Collection Date (Week of Year)", y = "Number of Samples",
       title = "Abundance + Testing Data Week of Year Distribution") +
  #scale_y_continuous(breaks = seq(0,2000, 200)) +
  scale_x_continuous(breaks = seq(0, 53, 2),
                     limits = c(0, 53))



# Find Start and End Dates ------------------------------------------------

# Create table of cluster : number of observations
densityFunction <- function(tibble) {

  topClusters <- tibble %>%
    group_by(clust) %>%
    # Find number of observations per cluster
    summarize(n = n()) %>%
    # Find top 10% largest clusters
    slice_max(order_by = n, prop = 0.1) %>%
    pull(clust)


  tibble %>%
    filter(clust %in% topClusters) %>%
    group_by(clust) %>%
    summarize(
      weeks = n_distinct(collection_date),
      obs = n()
    )



}




# Create Table of All possible Week Ranges --------------------------------



# Using woy vectors as input get number of observations remaining  when cluster
# must have all but 2 weeks as output
# with inputs in a table and plot

dateFunction <- function() {


  # Get table of all possible week range combinations
  woyVec <- seq(1, 52, 1)

  table <- expand.grid(woyVec, woyVec)

  # Filter to date ranges that make sense
  table %>%
    filter(Var1 < Var2) %>%
    mutate(
      start = Var1,
      end = Var2,
      .keep = "unused"
    )
}

# Test of above function
rawDates <- dateFunction()




# Filter Clusters by Week Ranges ------------------------------------------

# Amount of weeks permissible to be missing within range
tolerance <- 2L


# Filter clusters to high sampled within range
weekFilter <- function(start, end, tibble) {

  # Number of weeks (inclusive) in range
  cutoff <- end - start + 1

  # Clean observations to week range
  clean <- tibble %>%
    filter(woy >= start & woy <= end) %>%
    group_by(clust) %>%
    # Get number of weeks data present within range
    mutate(
      numWeeks = n_distinct(woy)
    ) %>%
    arrange(clust, woy)

  sumClean <- clean %>%
    distinct(clust) %>%
    nrow()

  # Filter to clusters with consistent sampling (only missing 2 weeks)
  filtered <- clean %>%
    filter(numWeeks >= (cutoff - tolerance))

  sumFiltered <- filtered %>%
    distinct(clust) %>%
    nrow()

  # Summary table
  out <- tribble(
    ~totalObs, ~consistentObs, ~totalClust, ~consistentClust,
    nrow(clean), nrow(filtered), sumClean, sumFiltered
  )

  # Add proportion of obs / clusters remaining as well as range info
  out %>%
    mutate(
      propObs = consistentObs / totalObs,
      propClust = consistentClust / totalClust,
      wStart = start,
      wEnd = end
    )

}

# # Test
# weekFilter(20,30, dataColTiny[[2]]) %>% view()
#
# test2 <- pmap(rawDates, weekFilter, dataColTiny[[2]])
#
# test2 %>%
#   bind_rows() %>%
#   filter(TotalObs > 0)


# Mapping function over both datasets

mapWeekFilter <- function(table) {

  # Map week filtering function
  out <- pmap(rawDates, weekFilter, table)

  # Bind rows to form tibble, filter out week ranges with no consistent clusters
  out %>%
    bind_rows() %>%
    filter(consistentObs > 0) %>%
    mutate(weeks = wEnd - wStart + 1)
}


weekData <- map(dataColTiny, mapWeekFilter)

# # Test
# weekData[[1]] %>% view()
#
# map(weekData, summary)


# Merge into 1 dataset ----------------------------------------------------



# Merge into 1 dataset
mergeFunction <- function(list) {

  ab <- list[[1]] %>%
    mutate(
      data = "Abundance"
    )

  tes <- list[[2]] %>%
    mutate(
      data = "Testing"
    )

  both <- bind_rows(ab, tes)


  both
  #weeksFunc(both)

}


merge <- mergeFunction(weekData)



# Plot With week filtering ------------------------------------------------




# Plot number of weeks vs total obs
weeksFunc <- function(tibble) {


  tibble %>%
    ggplot(mapping = aes(x = weeks, y = consistentObs, color = consistentClust)) +
    geom_point() +
    scale_x_continuous("Number of Weeks Sampled", breaks = seq(0, 52, 2)) +
    scale_y_continuous("Number of Observations from Clusters Consistently Sampled") +
    scale_color_viridis_c(
      breaks = seq(0, 500, 50),
      #option = "H"
    ) +
    #labs(title = data)+
    facet_grid(cols  = vars(data), scales = "free")


}

# Plot total obs vs Number of Clusters
weeksFunc1 <- function(tibble) {

  #filterDat <- filterFunc(tibble)


  tibble %>%
    ggplot(mapping = aes(y = consistentObs, x = consistentClust, color = weeks)) +
    geom_point() +
    scale_y_continuous(
      "Number of Observations from Clusters Consistently Sampled",
      #breaks = seq(0, 52, 2)
    ) +
    scale_x_continuous("Number of Consistent Clusters") +
    scale_color_viridis_c(
      #breaks = seq(0, 500, 50),
      #option = "H"
    ) +
    #labs(title = data)+
    facet_grid(cols  = vars(data), scales = "free")


}



# Plot for both datasets
weeksFunc(merge)




# Arrange Data for Viewing
arrangeFunc <- function(data) {
  data %>%
    arrange(desc(consistentObs), desc(weeks), desc(consistentClust))
}


map(weekData, arrangeFunc)




# Filter data to week ranges with most observations -----------------------


# Range with most observations
range <- merge %>%
  group_by(data) %>%
  arrange(desc(consistentObs)) %>%
  slice(1) %>%
  dplyr::select(wStart, wEnd) %>%
  mutate(
    start = wStart,
    end = wEnd,
    .keep = "unused"
  ) %>%
  ungroup()


# Range with most clusters
rangeClust <- merge %>%
  group_by(data) %>%
  arrange(desc(consistentClust)) %>%
  slice(1) %>%
  dplyr::select(wStart, wEnd) %>%
  mutate(
    start = wStart,
    end = wEnd,
    .keep = "unused"
  ) %>%
  ungroup()

# Range with longest time period

rangeWeeks <- merge %>%
  group_by(data) %>%
  arrange(desc(weeks)) %>%
  slice(1) %>%
  dplyr::select(wStart, wEnd) %>%
  mutate(
    start = wStart,
    end = wEnd,
    .keep = "unused"
  ) %>%
  ungroup()


rangeList <- list(obs = range, clust = rangeClust, weeks = rangeWeeks)




weekTable <- function(start, end, tibble) {

  # Number of weeks (inclusive) in range
  cutoff <- end - start + 1

  # Clean observations to week range
  clean <- tibble %>%
    filter(woy >= start & woy <= end) %>%
    group_by(clust) %>%
    # Get number of weeks data present within range
    mutate(
      numWeeks = n_distinct(woy)
    ) %>%
    arrange(clust, woy)



  # Filter to clusters with consistent sampling (only missing 2 weeks)
  filtered <- clean %>%
    filter(numWeeks >= (cutoff - tolerance))
}


# Format tibble for mapping pmap function

nest <- tibble(range, dataYearsFilter) %>%
  mutate(
    tibble = dataYearsFilter,
    .keep = "unused"
  ) %>%
  dplyr::select(-data)

nestMap <- function(weekRange) {

  tibble(weekRange, dataYearsFilter) %>%
    mutate(
      tibble = dataYearsFilter,
      .keep = "unused"
    ) %>%
    dplyr::select(-data)


}



# Format all date ranges and reduce data

weekTableMap <- function(list) {

  rangeMap <- map(list, nestMap)

  map(rangeMap, pmap, weekTable)

}


# All date ranges with filtered data
totalWeek <- weekTableMap(rangeList)


# Range with most observations
dataWeek <- totalWeek[[1]]#pmap(nest, weekTable)


# Get initial relationship models -----------------------------------------


# Testing aggregation

testWeek <- dataWeek[[2]]

# Find infection rates for each cluster / date
infectionRates <- testWeek %>%
  arrange(clust, desc(collection_date)) %>%
  # For every date and location get number of mosquitoes p / n
  group_by(clust, year, woy, collection_date, test_status) %>%
  summarise(
    traps = n(),
    count = sum(num_count)
  ) %>%
  # Get ratio of positive to negative
  summarize(
    total = sum(count),
    pos = sum(count[test_status == "Confirmed"]),
    testNegative = sum(traps[test_status == "Negative"]),
    testPositive = sum(traps[test_status == "Confirmed"]),
    totalTests = sum(traps)
  ) %>%
  mutate(
    rate = pos / total,
    posMosPerTest = pos / totalTests,
    posFreq = testPositive / totalTests
  )



#Abundance aggregation
abundWeek <- dataWeek[[1]]

abundanceRates <- abundWeek %>%
  ungroup() %>%
  # Sum mosquito numbers across columns
  mutate(
    # Set 0 trap night value to 1
    trap_nights = if_else(trap_nights == 0, 1, trap_nights),
    mos = rowSums(dplyr::select(., starts_with("culex"))),
    # Adjust mosquito numbers for number of traps / trap nights
    traps = num_trap * trap_nights,
    mosAdjust = mos / traps
  ) %>%
  # Get number of trapping nights, mosquitos by cluster
  group_by(clust, year, woy, collection_date) %>%
  summarise(
    points = n(),
    trapNights = sum(traps),
    mosAdjust = sum(mosAdjust)
  ) %>%
  # Final adjustment
  mutate(
    mosquitos = (mosAdjust / trapNights) %>% round(., 2)
  )



write_csv(testWeek, "./Data/7_dataDistribution/testData_weekFilter_mostObs.csv")
write_csv(abundWeek, "./Data/7_dataDistribution/abundData_weekFilter_mostObs.csv")


# # Some very high numbers for mosquito abundance but checks out with raw data
# abundanceRates %>%
#   filter(collection_date == mdy("08-06-2019")) %>%
#   view()




# Data Plots --------------------------------------------------------------

# Positive to Negative Test Frequency by Collection Date
infectionRates %>%
  group_by(collection_date) %>%
  summarize(
    posFreq = sum(testPositive, na.rm = TRUE) / sum(totalTests, na.rm = TRUE)
  ) %>%
  #slice_max(order_by = posFreq, prop = 0.1) %>%
  ggplot(mapping = aes(x = collection_date, y = posFreq)) +
  geom_point() +
  #geom_smooth() +
  scale_x_date(
    name = "Collection Date",
    date_breaks = "2 month",
    date_labels = "%m/%y",
    date_minor_breaks = "month"
  ) +
  scale_y_continuous(name = "Frequency of Positive Tests", #limits = c(0, 0.1),
                     breaks = seq(-1, 1, 0.05)) +
  labs(title = "Infection Rates Positive to Negative Test Frequency")



# Positive to Negative Test Frequency by  Week of Year
infectionRates %>%
  group_by(woy) %>%
  summarize(
    posFreq = sum(testPositive, na.rm = TRUE) / sum(totalTests, na.rm = TRUE),
    totalTests = sum(totalTests, na.rm = TRUE)
  ) %>%
  ggplot(mapping = aes(x = woy, y = posFreq, size = totalTests)) +
  geom_point(
    #shape = 1,
    color = "skyblue"
  ) +
  scale_x_continuous(name = "Week of Year",
                     breaks = seq(1, 52, 2)) +
  scale_y_continuous(name = "Frequency of Positive Tests", #limits = c(0, 0.1),
                     breaks = seq(-1, 1, 0.02)) +
  labs(title = "Infection Rates P/N Tests Frequency")

#
#
# abundanceRates %>%
#   mutate(
#     date = floor_date(collection_date, "2 week"),
#     week = week(date)
#   ) %>%
#   group_by(woy, date) %>%
#   summarize(mos = sum(mosquitos)) %>%
#   # filter(mos > 1000) %>%
#   # filter(mos < 5000) %>%
#   ggplot(mapping = aes(x = week, y = mos)) +
#   geom_line()
#
#





# Appendix ----------------------------------------------------------------
#
#
# # Plot Distribution -------------------------------------------------------
#
#
#
#
# plotDist <- function(data) {
#
#   data %>%
#     ggplot(mapping = aes(x = wStart,
#                          y = wEnd,
#                          size = consistentObs,
#                          color = consistentClust
#     )) +
#     geom_point(
#       #color = "skyblue",
#       #alpha = 0.4
#     ) +
#     #scale_size_area(TotalObs) +
#     scale_x_continuous(name = "Week of Year Start", breaks = seq(0, 52, 2)) +
#     scale_y_continuous(name = "Week of Year End", breaks = seq(0, 52, 2))
#
#
# }
# #
# #
# # abundance %>% filterFunc() %>% plotDist()
# #
# # testing %>% filterFunc() %>% plotDist()
#

#
# test <- map(dataColTiny, densityFunction)
#
#
# test[[2]]
#
# test[[2]] %>%
#   ggplot(mapping = aes(x = weeks)) +
#   geom_histogram(binwidth = 1,
#                  color = "white",
#                  fill = "skyblue",
#                  alpha = 0.7,
#                  position = 'dodge'
#   ) +
#   labs(x = "Number of weeks Sampled", y = "Number of Clusters",
#        title = "Testing Data top 10% of Clusters") +
#   scale_y_continuous(breaks = seq(0,10, 2)) +
#   scale_x_continuous(breaks = seq(0, 125, 5))
#
#
# test[[2]] %>%
#   ggplot(mapping = aes(x = weeks, y = obs)) +
#   geom_jitter(#binwidth = 1,
#                  #color = "white",
#                  #fill = "red",
#                  #alpha = 0.7,
#                  #position = 'dodge'
#   ) +
#   labs(x = "Number of weeks Sampled", y = "Number of Observations",
#        title = "Testing Data top 10% of Clusters") +
#   #scale_y_continuous(breaks = seq(0,10, 2)) +
#   scale_x_continuous(breaks = seq(0, 125, 5))

#
# # id : Cluster key
# clustKey <- read_csv("idClusterKey.csv")
#
#
#
# # Lat / lon file
# geoID <- read_csv("idPoints.csv") %>%
#   st_as_sf(coords = c("lon", "lat"))
#
#
#
# # Check for NA's
#
#
# naColFunction <- function(tibble) {
#
#   names(which(colSums(is.na(tibble)) > 0))
#
#
#
# }
#
#
# naColFunction(abund)
#
# naColFunction(test)
#
#
# # Spatial Join ------------------------------------------------------------
#
#
#
# addGeo <- function(tibble) {
#
#   tibble %>%
#     left_join(clustKey, by = "clust") %>%
#     left_join(geoID, by = "id") %>%
#     st_as_sf()
#
#
#
# }
#
#
# map(dataList, addGeo)
#
#
#
# spatial <- abundSpatial %>%
#   bind_rows(testSpatial) %>%
#   arrange(clust) %>%
#   filter(duplicated(clust) == FALSE)
#
#
#
# # Plot Spatial Distribution of weekrange Filters --------------------------
#
#
#
# # Clean and Filter spatial Cluster Data
#
# cleanGeoFunc <- function(table) {
#
#   clustNum <- table %>%
#     pull(clust) %>%
#     unique()
#
#
#   spatial %>%
#     filter(clust %in% clustNum)
#
# }
#
#
#
# # Name and clean function
# nameFunc <- function(list) {
#
#   # Filter spatial data
#   cleanGeo <- map(list, cleanGeoFunc)
#
#
#   # Name lists
#   list1 <- cleanGeo[[1]] %>%
#     mutate(
#       data = "Abundance"
#     )
#
#   list2 <- cleanGeo[[2]] %>%
#     mutate(
#       data = "Testing"
#     )
#
#   # Merge
#   newList <- list1 %>%
#     bind_rows(list2)
#
#
#
#   newList
#
#
# }
#
#
#
#
# # Clean plots and combine and map
# namePlotFunc <- function(list) {
#
#
#   cleanList <- map(list, nameFunc)
#
#
#   list1 <- cleanList[[1]] %>%
#     mutate(
#       data2 = "Most observations"
#     )
#
#   list2 <- cleanList[[2]] %>%
#     mutate(
#       data2 = "Most clusters"
#     )
#
#   list3 <- cleanList[[3]] %>%
#     mutate(
#       data2 = "Most weeks"
#     )
#
#   newList <- list1 %>%
#     bind_rows(list2) %>%
#     bind_rows(list3)
#
#
#
#   newList %>%
#     tm_shape() +
#     tm_dots() +
#     tm_facets(by = c("data2", "data"), sync = TRUE) +
#     tm_basemap(server = "OpenStreetMap")
#
# }
#
#
# # Map of comparisons
# namePlotFunc(totalWeek)
#
#





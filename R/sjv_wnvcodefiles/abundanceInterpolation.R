
library(zoo)
library(sf)
library(tidyverse)
library(tmap)
library(lubridate)
library(scales)
library(splines2)


# Read in Data ------------------------------------------------------------




# id : Cluster key
clustKey <- read_csv("Data/KMVCD/wnv/3_clusterExport/idClusterKey.csv")



# Abundance by species and lifecycle
abund <- read_csv("Data/KMVCD/wnv/5_summarizeClusters/abundSpeciesSummary.csv") %>%
  dplyr::select(-wy)

abundAll <- abund %>%
  group_by(clust, woy, year, collection_date) %>%
  summarize(
    points = sum(points),
    trapNights = sum(trapNights),
    totalMosquitoes = sum(totalMosquitoes),
    mosPerTrapNight = sum(mosPerTrapNight)
  )

abundPipiensAll <- abund %>%
  filter(species == "pipiens")

abundPipiensF <- abundPipiensAll %>%
  filter(lifeCycle == "female")


abundTarsalisAll <- abund %>%
  filter(species == "tarsalis")


abundTarsalisF <-  abundTarsalisAll %>%
  filter(lifeCycle == "female")


# WNV test by species and lifecycle
test <- read_csv("Data/KMVCD/wnv/5_summarizeClusters/testPIRMIRClusterTEST.csv")

# Split by year -----------------------------------------------------------




abundDataSets <- list(abundAll, #abundPipiensAll, abundPipiensF,
                      abundTarsalisAll, abundTarsalisF)






# Extrapolate Starting and Ending Values ----------------------------------



# Find start and end week range for positive results for each year
# Add these weeks as observations (where necessary) to extrapolate
# trailing / leading values before substituting 0's


# For abundance data
setSampleRangeAbundance <- function(tibble) {

  # Subset table to first weeks with non 0 results
  nonZero <- tibble %>%
    filter(totalMosquitoes > 0)

  # Get first week of year
  start <- min(nonZero$woy)
  # Last week of year
  end <- max(nonZero$woy)

  # Get table of min / max weeks for each cluster
  # If outside the range of positive results, remove
  # This yields table with
  toBeTailed <- tibble %>%
    group_by(clust) %>%
    # Get first and last week cluster sampled
    summarize(
      min = min(woy),
      max = max(woy)
    ) %>%
    # If outside pos range set to -99
    mutate(
      min = if_else(min <= start, -99, min),
      max = if_else(max >= end, -99, max)
    ) %>%
    # Pivot and filter to get table of minimum and maximum
    # that needs to be imputed per each cluster
    pivot_longer(
      cols = min:max,
      names_to = "tail",
      values_to = "woy"
    ) %>%
    filter(woy > 0) %>%
    dplyr::select(-woy)

  # Convert "min" and "max" to range woy
  tail <- tribble(
    ~tail, ~woy,
    "min", start,
    "max", end
  )

  # Join to main table with 0's for mosquito values at
  # boundaries
  joinTable <- toBeTailed %>%
    left_join(tail, by = "tail") %>%
    dplyr::select(-tail) %>%
    mutate(
      totalMosquitoes = 0,
      mosPerTrapNight = 0
    )


  tibble %>%
    bind_rows(joinTable) %>%
    fill(year, .direction = "updown") %>%
    arrange(clust, year, woy)



}



abundExtrapolate <- map(abundDataSets, setSampleRange) %>%
  bind_rows() %>%
  group_split(clust, year)

# For testing data

setSampleRangeTest <- function(tibble) {

  # Subset table to first weeks with non 0 results
  nonZero <- tibble %>%
    filter(Pos > 0)

  # Get first week of year
  start <- min(nonZero$woy)
  # Last week of year
  end <- max(nonZero$woy)

  # Get table of min / max weeks for each cluster
  # If outside the range of positive results, remove
  # This yields table with
  toBeTailed <- tibble %>%
    group_by(clust) %>%
    # Get first and last week cluster sampled
    summarize(
      min = min(woy),
      max = max(woy)
    ) %>%
    # If outside pos range set to -99
    mutate(
      min = if_else(min <= start, -99, min),
      max = if_else(max >= end, -99, max)
    ) %>%
    # Pivot and filter to get table of minimum and maximum
    # that needs to be imputed per each cluster
    pivot_longer(
      cols = min:max,
      names_to = "tail",
      values_to = "woy"
    ) %>%
    filter(woy > 0) %>%
    dplyr::select(-woy)

  # Convert "min" and "max" to range woy
  tail <- tribble(
    ~tail, ~woy,
    "min", start,
    "max", end
  )

  # Join to main table with 0's for mosquito values at
  # boundaries
  joinTable <- toBeTailed %>%
    left_join(tail, by = "tail") %>%
    dplyr::select(-tail) %>%
    mutate(
      totalMosquitoes = 0,
      mosPerTrapNight = 0
    )


  tibble %>%
    bind_rows(joinTable) %>%
    fill(year, .direction = "updown") %>%
    arrange(clust, year, woy)



}

testExtrapolate <- setSampleRangeTest(test) %>%
  bind_rows() %>%
  group_split(clust, year)


# Split into lists --------------------------------------------------------


# Split abundance by cluster and year

abundSplit <- abund %>%
  group_split(clust, year)


# Split test by cluster and year
testSplit <- test %>%
  group_split(clust, year)

# Interpolate -------------------------------------------------------------


# Create vector of every woy between start and end sample dates
minMaxFunc <- function(table) {

  min <- min(table$woy)

  max <- max(table$woy)

  seq(min, max, 1)
}


# Abundance Interpolations ------------------------------------------------


abundFunc <- function(table) {

  # table_distinct <- table %>%
  #   dplyr::dplyr::select(woy, totalMosquitoes) %>%
  #   distinct()

  # If else to account for clusters with only 1 obs in a year
  if(nrow(table) > 1) {

    # Create Spline functions

    # If there is an insufficient number of slopes m_i, then copy data over to interpolated columns

    out <- tryCatch(
      {
        # Get weeks of year in range for
        weeks <- table %>%
          minMaxFunc() %>%
          enframe(name = NULL, value = "woy")



        # Weeks with Data
        tableWeeks <- table %>%
          distinct(woy) %>%
          full_join(weeks, by = "woy")



        # Missing weeks of data
        missingWeeks <- weeks %>%
          setdiff(tableWeeks)


        # Join
        joinTable <- table %>%
          full_join(weeks, by = "woy") %>%
          arrange(woy) %>%
          fill(clust) %>%
          fill(year)

        # Linear Interpolations of Missing Data
        linTable <- joinTable %>%
          mutate(
            totalLinear = na.approx(totalMosquitoes),
            mosPerLinear = na.approx(mosPerTrapNight),
            #collection_date = na.approx(collection_date) %>% as.Date()
          )

        totalSplineFun <- splinefun(table$woy,
                                    table$totalMosquitoes,
                                    method = "monoH.FC")

        densitySplineFun <- splinefun(table$woy,
                                      table$mosPerTrapNight,
                                      method = "monoH.FC")

        # Spline interpolations of Missing Data
        totalFill <- linTable %>%
          mutate(
            totalSpline = totalSplineFun(woy) %>% round(digits = 1),
            mosPerSpline = densitySplineFun(woy) %>% round(digits = 2)
          )

        return(totalFill)
      },
      error = function(e) {

        clean <- table %>%
          mutate(
            totalLinear = totalMosquitoes,
            mosPerLinear = mosPerTrapNight,
            totalSpline = totalMosquitoes,
            mosPerSpline = mosPerTrapNight
          )
        return(clean)
      }
    )

    out

  } else {
    print(unique(table$clust))


    # If only one year copy data over to interpolated columns
    clean <- table %>%
      mutate(
        totalLinear = totalMosquitoes,
        mosPerLinear = mosPerTrapNight,
        totalSpline = totalMosquitoes,
        mosPerSpline = mosPerTrapNight
      )



    return(clean)

  }


}

# Interpolate for non-added tailing values


# Map over all clusters / years
abundInterpolate <- map(abundSplit, abundFunc)
# Bind to one dataset
abundIntFinal <- abundInterpolate %>%
  bind_rows() %>%
  unite(col = "wy", c(year, woy), remove = FALSE)



# Interpolate for added tailing values


# Map over all clusters / years
abundExt <- map(abundExtrapolate, abundFunc)
# Bind to one dataset
abundExtFinal <- abundExt %>%
  bind_rows()



# Testing Interpolation ---------------------------------------------------



testFunc <- function(table) {

  # If else to account for clusters with only 1 obs in a year
  if(nrow(table) > 1) {
    print(table$clust %>% unique())

    out <- tryCatch({
      # Get weeks of year in range for
      weeks <- table %>%
        minMaxFunc() %>%
        enframe(name = NULL, value = "woy")



      # Weeks with Data
      tableWeeks <- table %>%
        distinct(woy) %>%
        full_join(weeks, by = "woy")



      # Missing weeks of data
      missingWeeks <- weeks %>%
        setdiff(tableWeeks)


      # Join
      joinTable <- table %>%
        full_join(weeks, by = "woy") %>%
        arrange(woy) %>%
        fill(clust) %>%
        fill(year)

      # Linear Interpolations of Missing Data
      linTable <- joinTable %>%
        mutate(
          rateLinear = na.approx(MIR),
          mosPerTestLinear = na.approx(Pos),
          testFreqLinear = na.approx(posTestFreq),
          #collection_date = na.approx(collection_date) %>% as.Date()
        )

      # Create Spline functions
      rateSplineFun <- splinefun(table$woy,
                                 table$MIR,
                                 method = "monoH.FC")

      mosPerTestSplineFun <- splinefun(table$woy,
                                       table$Pos,
                                       method = "monoH.FC")

      testFreqSplineFun <- splinefun(table$woy,
                                     table$posTestFreq,
                                     method = "monoH.FC")

      # Spline interpolations of Missing Data
      totalFill <- linTable %>%
        mutate(
          rateSpline = rateSplineFun(woy), # %>% round(digits = 1),
          mosPerTestSpline = mosPerTestSplineFun(woy), # %>% round(digits = 2),
          testFreqSpline = testFreqSplineFun(woy), # %>% round(digits = 1),
        )

      return(totalFill)
    }, error = function(e) {
      clean <- table %>%
        mutate(
          rateLinear = MIR,
          mosPerTestLinear = Pos,
          testFreqLinear = posTestFreq,
          rateSpline = MIR,
          mosPerTestSpline = Pos,
          testFreqSpline = posTestFreq,
        )

      return(clean)
      }
    )

    out

  } else {


    # If only one year copy data over to interpolated columns
    clean <- table %>%
      mutate(
        rateLinear = MIR,
        mosPerTestLinear = Pos,
        testFreqLinear = posTestFreq,
        rateSpline = MIR,
        mosPerTestSpline = Pos,
        testFreqSpline = posTestFreq,
      )



    return(clean)

  }

}



# Interpolate for non-added tailing values



# Map over all clusters / years
testInterpolate <- map(testSplit, testFunc)
# Bind into one dataset
testIntFinal <- testInterpolate %>%
  bind_rows() %>%
  unite(col = "wy", c(year, woy), remove = FALSE)


# Interpolate for added tailing values


# Map over all clusters / years
testExt <- map(testExtrapolate, testFunc)
# Bind to one dataset
testExtFinal <- testExt %>%
  bind_rows()



# Infer 0's for rest of year ----------------------------------------------


# Make grid of all possible dates as well as clusters to fill 0's in
# extrapolated data

# Vector of all clusters
clustsAbund <- abund %>%
  #bind_rows(test) %>%
  pull(clust) %>%
  unique()

clustsTest <- test %>%
  pull(clust) %>%
  unique()

# Vector of all years
yrs <- seq(2010, 2023, 1)

# Vector of all weeks for year
wks <- seq(1, 52, 1)



# Create all possible combinations of year + week
allWeeks <- expand_grid(yrs, wks) %>%
  mutate(
    yrs = as.character(yrs),
    wks = as.character(wks)
  )



abundWClust <- allWeeks %>%
  unite(col = "wy", c(yrs, wks)) %>%
  expand_grid(clustsAbund)

testWClust <- allWeeks %>%
  unite(col = "wy", c(yrs, wks)) %>%
  expand_grid(clustsTest)






abundWeeks <- abundExtFinal %>%
  unite(col = "wy", c(year, woy)) %>%
  full_join(abundWClust, by = c("wy", "clust" = "clustsAbund")) %>%
  separate(wy, c("year", "woy"), remove = FALSE) %>%
  mutate(
    year = as.numeric(year),
    woy = as.numeric(woy)
  ) %>%
  dplyr::select(clust, year, woy, everything()) %>%
  mutate(
    across(totalLinear:mosPerSpline, ~replace_na(., 0))
  )


testWeeks <- testExtFinal %>%
  unite(col = "wy", c(year, woy)) %>%
  full_join(testWClust, by = c("wy", "clust" = "clustsTest")) %>%
  separate(wy, c("year", "woy"), remove = FALSE) %>%
  mutate(
    year = as.numeric(year),
    woy = as.numeric(woy)
  ) %>%
  dplyr::select(clust, year, woy, everything()) %>%
  mutate(
    across(rateLinear:testFreqSpline, ~replace_na(., 0))
  )


# Export ------------------------------------------------------------------


write_csv(abundIntFinal, "Data/KMVCD/wnv/6_interpolateData/abundClusterInterpolation.csv")
write_csv(testIntFinal, "Data/KMVCD/wnv/6_interpolateData/testClusterInterpolation.csv")


write_csv(abundWeeks, "Data/KMVCD/wnv/6_interpolateData/abundClusterInterpolationInferNA.csv")
write_csv(testWeeks, "Data/KMVCD/wnv/6_interpolateData/testClusterInterpolationInferNA.csv")



# Appendix ----------------------------------------------------------------

# # Convert year and week of year to date function
# convertDate <- function(yrs, wks) {
#
#   yr = str_c(yrs, "0101") %>% ymd()
#
#   week(yr) <- wks
#
#   yr
#
#
# }
# %>%
#   # Change to collection_date
#   mutate(
#     collection_date = convertDate(yrs, wks)
#   )

#
# # Why are 10 rows being added when expanding for all weeks? in abundance data
#
# we1 <- abundWClust %>%
#   distinct(wy, clustsAbund) %>%
#   mutate(clust = clustsAbund, .keep = "unused")
#
# glimpse(abundWClust)
#
#
# we2 <- abundWeeks %>%
#   distinct(wy, clust)
#
# # There are 53 weeks in 2017 for 10 clusters...why
# weirdClust <- setdiff(we2, we1) %>%
#   pull(clust)
#
#
#
# abundWeeks %>%
#   filter(clust %in% weirdClust) %>%
#   filter(woy == 53) %>%
#   view()
#
#
# abundWeeks %>%
#   filter(is.na(totalSpline)) %>%
#   arrange(clust, year, woy) %>% view()
#
# # Looked up online and some years can have 53 weeks, 2012 and 2017 within
# # our range.



# Compare Interpolations --------------------------------------------------
#
#
#
#
#
# viewData <- abundWeeks %>%
#   filter(clust == 230) %>%
#   arrange(year, woy)
#
#
#
#
# plotInterpolation <- function(cluster, yr) {
#
#   section <- abundWeeks %>%
#     filter(clust == cluster) %>%
#     filter(year == yr) %>%
#     pivot_longer(
#       cols = totalMosquitoes:mosPerSpline,
#       names_to = "Data",
#       values_to = "Values"
#     ) %>%
#     mutate(
#       method  = case_when(str_detect(Data, "Linear$") ~ "linear",
#                           str_detect(Data, "Spline$") ~ "spline",
#                           .default = "raw") %>%
#         as_factor(),
#       dat = case_when(str_detect(Data, "^total") ~ "totalMos",
#                       str_detect(Data, "^mos") ~ "mosPerTrapNight",
#                       .default = "ERROR") %>%
#         as_factor()
#     )
#
#
#   test <- section %>%
#     filter(dat == "totalMos")
#
#
#   # Create Spline functions
#   totalSplineFun <- splinefun(abundWeeks$woy,
#                               abundWeeks$totalMosquitoes,
#                               method = "monoH.FC")
#
#   densitySplineFun <- splinefun(abundWeeks$woy,
#                                 abundWeeks$mosPerTrapNight,
#                                 method = "monoH.FC")
#
#
#
#   test %>%
#     ggplot(mapping = aes(x = woy, y = Values)) +
#     #geom_point() #+
#     geom_point(data = filter(test, method == "raw")) +
#     geom_line(data = filter(test, method == "linear"), color = "red") +
#     #stat_function(fun =totalSplineFun, args = list(deriv = 0))
#     geom_line(data = filter(test, method == "spline"))
#
#
#   # section %>%
#   #   filter(dat == "totalMos") %>%
#   #   filter(method == "linear") %>%
#   #   arrange(woy)
#
# }
#
#
#
# plotInterpolation(230, 2011) #%>% view()
#
#
#
#
# plot2 <- function(table) {
#
#
#   # Get weeks of year in range for
#   weeks <- table %>%
#     minMaxFunc() %>%
#     enframe(name = NULL, value = "woy")
#
#
#
#   # Weeks with Data
#   tableWeeks <- table %>%
#     distinct(woy) %>%
#     full_join(weeks, by = "woy")
#
#
#
#   # Missing weeks of data
#   missingWeeks <- weeks %>%
#     setdiff(tableWeeks)
#
#
#   # Join
#   joinTable <- table %>%
#     full_join(weeks, by = "woy") %>%
#     arrange(woy) %>%
#     fill(clust) %>%
#     fill(year)
#
#   # Linear Interpolations of Missing Data
#   linTable <- joinTable %>%
#     mutate(
#       totalLinear = na.approx(totalMosquitoes),
#       mosPerLinear = na.approx(mosPerTrapNight),
#       collection_date = na.approx(collection_date) %>% as.Date()
#     )
#
#   # Create Spline functions
#   totalSplineFun <- splinefun(table$woy,
#                               table$totalMosquitoes,
#                               method = "monoH.FC")
#
#   densitySplineFun <- splinefun(table$woy,
#                                 table$mosPerTrapNight,
#                                 method = "monoH.FC")
#
#   # Spline interpolations of Missing Data
#   totalFill <- linTable %>%
#     mutate(
#       totalSpline = totalSplineFun(woy) %>% round(digits = 1),
#       mosPerSpline = densitySplineFun(woy) %>% round(digits = 2)
#     )
#
#
#
#   ggplot(data.frame(x = table$woy, y = table$totalMosquitoes), aes(x,y)) +
#     stat_function(fun = totalSplineFun, color = "orange") +
#     geom_point() +
#     labs(x = "week of year", y = "total Mosquitoes") +
#     scale_y_continuous(limits = c(0, 100))
#
#
# }
#
#
#
# abundSplit[[2168]] %>% plot2
#



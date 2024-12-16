
library(zoo)
library(sf)
library(tidyverse)
library(tmap)
library(lubridate)
library(scales)

library(splines2)

dir.create("Data/SJV/wnv/6_interpolateData/")


# Read in Data ------------------------------------------------------------


# Abundance by species and lifecycle
abund <- read_csv("Data/SJV/wnv/5_summarizeClusters/abundSpeciesSummary.csv") %>%
  dplyr::select(-wy)


# Sum Function
summarizeFunction <- function(tibble) {

  tibble %>%
    group_by(clust, woy, year, collection_date) %>%
    summarize(
      points = sum(points),
      trapNights = sum(trapNights),
      totalMosquitoes = sum(totalMosquitoes),
      mosPerTrapNight = totalMosquitoes/trapNights
    )

}

abundAll <- abund %>%
  summarizeFunction()

abund$species %>% unique()

abundPipiensAll <- abund %>%
  filter(species == "pipiensMix") %>%
  summarizeFunction()

abundPipiensF <- abund %>%
  filter(species == "pipiensMix") %>%
  filter(lifeCycle == "female") %>%
  summarizeFunction()


abundTarsalisAll <- abund %>%
  filter(species == "tarsalis") %>%
  summarizeFunction()


abundTarsalisF <-  abund %>%
  filter(species == "tarsalis") %>%
  filter(lifeCycle == "female") %>%
  summarizeFunction()

abundQuinquefasciatusAll <- abund %>%
  filter(species == "quinquefasciatus") %>%
  summarizeFunction()

abundQuinquefasciatusF <- abund %>%
  filter(species == "quinquefasciatus") %>%
  filter(lifeCycle == "female") %>%
  summarizeFunction()

abundAegyptiAll <- abund %>%
  filter(species == "aegypti") %>%
  summarizeFunction()

abundAegyptiF <- abund %>%
  filter(species == "aegypti") %>%
  filter(lifeCycle == "female") %>%
  summarizeFunction()




# Split by year -----------------------------------------------------------




abundDataSets <- list(abundAll, abundPipiensAll, abundPipiensF,
                      abundTarsalisAll, abundTarsalisF)






# Extrapolate Starting and Ending Values ----------------------------------



# Find start and end week range for positive results for each year
# Add these weeks as observations (where necessary) to extrapolate
# trailing / leading values before substituting 0's


# For abundance data
setSampleRange <- function(tibble) {

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




# Split into lists --------------------------------------------------------


# Split abundance by cluster and year

abundSplit <- abund %>%
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


  # If else to account for clusters with only 1 obs in a year
  if(nrow(table) > 1) {
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

    # Create Spline functions
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

  } else {


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




# Interpolate for added tailing values


# Map over all clusters / years
abundExt <- function(tibble) {

  out1 <- tibble %>%
    group_by(year) %>%
    group_split() %>%
    map(setSampleRange) %>%
    bind_rows() %>%
    group_split(clust, year) %>%
    map(abundFunc) %>%
    bind_rows()


  # Infer 0's for rest of year ----------------------------------------------


  # Make grid of all possible dates as well as clusters to fill 0's in
  # extrapolated data

  # Vector of all clusters
  clustsAbund <- tibble %>%
    #bind_rows(test) %>%
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



  abundWeeks <- out1 %>%
    unite(col = "wy", c(year, woy)) %>%
    full_join(abundWClust, by = c("wy", "clust" = "clustsAbund")) %>%
    separate(wy, c("year", "woy"), remove = FALSE) %>%
    mutate(
      year = as.numeric(year),
      woy = as.numeric(woy)
    ) %>%
    dplyr::select(clust, year, woy, everything())


  abundWeeks


}




allExt <- abundExt(abundAll)
pipAllExt <- abundExt(abundPipiensAll)
pipFExt <- abundExt(abundPipiensF)
tarAllExt <- abundExt(abundTarsalisAll)
tarFExt <- abundExt(abundTarsalisF)
# quinxAllExt <- abundExt(abundQuinquefasciatusAll)
# quinxFExt <- abundExt(abundQuinquefasciatusF)
# aegAllExt <- abundExt(abundAegyptiAll)
# aegFExt <- abundExt(abundAegyptiF)


#map(abundDataSets, abundInterpolate)


# Export ------------------------------------------------------------------



write_csv(allExt, "Data/SJV/wnv/6_interpolateData/allExt_NA.csv")
write_csv(pipAllExt, "Data/SJV/wnv/6_interpolateData/pipiensAllExt_NA.csv")
write_csv(pipFExt, "Data/SJV/wnv/6_interpolateData/pipiensFemaleExt_NA.csv")
write_csv(tarAllExt, "Data/SJV/wnv/6_interpolateData/tarsalisAllExt_NA.csv")
write_csv(tarFExt, "Data/SJV/wnv/6_interpolateData/tarsalisFemaleExt_NA.csv")
# write_csv(quinxAllExt, "Data/SJV/wnv/6_interpolateData/quinquefasciatusAllExt_NA.csv")
# write_csv(quinxFExt, "Data/SJV/wnv/6_interpolateData/quinquefasciatusFemaleExt_NA.csv")
# write_csv(aegAllExt, "Data/SJV/wnv/6_interpolateData/aegyptiAllExt_NA.csv")
# write_csv(aegFExt, "Data/SJV/wnv/6_interpolateData/aegyptiFemaleExt_NA.csv")

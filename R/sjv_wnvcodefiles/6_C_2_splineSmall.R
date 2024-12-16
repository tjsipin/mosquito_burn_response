
library(zoo)
library(sf)
library(tidyverse)
library(tmap)
library(lubridate)
library(scales)
library(hrbrthemes)
library(splines2)
library(sitar)


# Read in Data ------------------------------------------------------------


# Abundance by species and lifecycle
abund <- read_csv("Data/5_summarizeClusters/abundSpeciesSummary.csv") %>%
  dplyr::select(-wy)


# Sum Function
summarizeFunction <- function(tibble) {

  tibble %>%
    group_by(clust, woy, year) %>%
    summarize(
      points = sum(points),
      trapNights = sum(trapNights),
      totalMosquitoes = sum(totalMosquitoes),
      mosPerTrapNight = sum(mosPerTrapNight)
    )

}


abundAll <- abund %>%
  summarizeFunction()

# abundPipiensAll <- abund %>%
#   filter(species == "pipiens") %>%
#   summarizeFunction()

# abundPipiensF <- abund %>%
#   filter(species == "pipiens") %>%
#   filter(lifeCycle == "female") %>%
#   summarizeFunction()


abundTarsalisAll <- abund %>%
  filter(species == "tarsalis") %>%
  summarizeFunction()


abundTarsalisF <-  abund %>%
  filter(species == "tarsalis") %>%
  filter(lifeCycle == "female") %>%
  summarizeFunction()




# Split by year -----------------------------------------------------------




abundDataSets <- list(abundAll, #abundPipiensAll, abundPipiensF,
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

  table2 <- table %>%
    mutate(
      mosPerLog = log(1 + mosPerTrapNight)
    )


  # If else to account for clusters with only 1 obs in a year
  if(nrow(table2) > 1) {
    # Get weeks of year in range for
    weeks <- table %>%
      minMaxFunc() %>%
      enframe(name = NULL, value = "woy")



    # Weeks with Data
    tableWeeks <- table2 %>%
      distinct(woy) %>%
      full_join(weeks, by = "woy")



    # Missing weeks of data
    missingWeeks <- weeks %>%
      setdiff(tableWeeks)


    # Join
    joinTable <- table2 %>%
      full_join(weeks, by = "woy") %>%
      arrange(woy) %>%
      fill(clust) %>%
      fill(year)

    # Linear Interpolations of Missing Data
    linTable <- joinTable %>%
      mutate(
        #totalLinear = na.approx(totalMosquitoes),
        mosPerLinear = na.approx(mosPerTrapNight),

        #collection_date = na.approx(collection_date) %>% as.Date()
      )

    # Create Spline functions
    # totalSplineFun <- splinefun(table$woy,
    #                             table$totalMosquitoes,
    #                             method = "monoH.FC")

    densitySplineFun <- splinefun(table2$woy,
                                  table2$mosPerLog,
                                  method = "monoH.FC")


    #return(densitySplineFun)

    # Spline interpolations of Missing Data
    totalFill <- linTable %>%
      mutate(
        #totalSpline = totalSplineFun(woy) %>% round(digits = 1),
        mosPerSpline = densitySplineFun(woy) %>% round(digits = 2)
      )

    return(totalFill)

  } else {


    # If only one obs copy data over to interpolated columns
    clean <- table %>%
      mutate(
        totalLinear = totalMosquitoes,
        mosPerLinear = mosPerTrapNight,
        #totalSpline = totalMosquitoes,
        mosPerSpline = mosPerTrapNight
      )


    # Test of function
    return("poop")

    #return(clean)

  }


}





# Map over all clusters / years
abundExt <- function(tibble) {

  out1 <- tibble %>%
    group_by(year) %>%
    group_split() %>%
    map(setSampleRange) %>%
    bind_rows() %>%
    group_split(clust, year) %>%
    map(abundFunc) #%>%
  #bind_rows()

  return(out1)

  # # Infer 0's for rest of year ----------------------------------------------
  #
  #
  # # Make grid of all possible dates as well as clusters to fill 0's in
  # # extrapolated data
  #
  # # Vector of all clusters
  # clustsAbund <- tibble %>%
  #   #bind_rows(test) %>%
  #   pull(clust) %>%
  #   unique()
  #
  #
  # # Vector of all years
  # yrs <- seq(2010, 2020, 1)
  #
  # # Vector of all weeks for year
  # wks <- seq(1, 52, 1)
  #
  #
  #
  # # Create all possible combinations of year + week
  # allWeeks <- expand_grid(yrs, wks) %>%
  #   mutate(
  #     yrs = as.character(yrs),
  #     wks = as.character(wks)
  #   )
  #
  #
  #
  # abundWClust <- allWeeks %>%
  #   unite(col = "wy", c(yrs, wks)) %>%
  #   expand_grid(clustsAbund)
  #
  #
  #
  # abundWeeks <- out1 %>%
  #   unite(col = "wy", c(year, woy)) %>%
  #   full_join(abundWClust, by = c("wy", "clust" = "clustsAbund")) %>%
  #   separate(wy, c("year", "woy"), remove = FALSE) %>%
  #   mutate(
  #     year = as.numeric(year),
  #     woy = as.numeric(woy)
  #   ) %>%
  #   dplyr::select(clust, year, woy, everything()) %>%
  #   mutate(
  #     across(totalLinear:mosPerSpline, ~replace_na(., 0))
  #   )
  #
  #
  # abundWeeks
  #

}


# Check Interpolations ----------------------------------------------------

testData <- abundTarsalisF %>%
  filter(clust == 15, year == "2020")

testData %>%
  view()


df <- sitar::dfset(testData$woy,
                   testData$mosPerTrapNight,
                   table2,
                   FUN = BIC)


nsFit <- naturalSpline()



x <- seq.int(0, 1, 0.01)
knots <- c(0.3, 0.5, 0.6)

## natural spline basis
nsMat0 <- naturalSpline(x, knots = knots, intercept = TRUE)
## integrals
nsMat1 <- naturalSpline(x, knots = knots, intercept = TRUE, integral = TRUE)
## first derivatives
nsMat2 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 1)
## second derivatives
nsMat3 <- naturalSpline(x, knots = knots, intercept = TRUE, derivs = 2)

op <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
matplot(x, nsMat0, type = "l", ylab = "basis")
matplot(x, nsMat1, type = "l", ylab = "integral")
matplot(x, nsMat2, type = "l", ylab = "1st derivative")
matplot(x, nsMat3, type = "l", ylab = "2nd derivative")
par(op) # reset to previous plotting settings

## use the deriv method
all.equal(nsMat0, deriv(nsMat1), check.attributes = FALSE)
all.equal(nsMat2, deriv(nsMat0))
all.equal(nsMat3, deriv(nsMat2))
all.equal(nsMat3, deriv(nsMat0, 2))



pipFExt <- abundExt(testData)

asTibble <- pipFExt %>%
  bind_rows()


asTibble %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot() +
  geom_line(mapping = aes(x = woy, y = mosPerSpline, col = 'orange')) +
  geom_point(mapping = aes(x = woy, y = mosPerTrapNight, col = 'skyblue')) +
  #facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  xlab("Week of Year") +
  scale_x_continuous(breaks = seq(0, 54, 2)) +
  ylab("Mosquitoes Trap Night") +
  labs(title = "2020 Total Mosquitoes Cluster 15")






# Packages ----------------------------------------------------------------

library(zoo)
library(sf)
library(tidyverse)
library(tmap)
library(lubridate)
library(scales)
library(splines2)

tmap_mode("view")

dir.create("Data/SJV/wnv/6_interpolateData/")

# Read in Data ------------------------------------------------------------


# WNV testing data
test <- read_csv("Data/SJV/wnv/5_summarizeClusters/testPIRMIRClusterTEST.csv")

#test <- read_csv("./Data/5_summarizeClusters/testPIRMIRCluster.csv")

# Subset by species -------------------------------------------------------


glimpse(test)

test %>%
  distinct(species)

tar <- test %>%
  filter(species == "Culex tarsalis")

pip <- test %>%
  filter(species == "Culex pipiens")


all <- test %>%
  group_by(clust, woy, year) %>%
  summarize(
    PoolSize = sum(PoolSize),
    Pos = sum(Pos),
    NumPools = sum(NumPools),
    P = mean(P)
  ) %>%
  mutate(
    MIR = (Pos / PoolSize) * 1000,
    posTestFreq = Pos / NumPools
  ) %>%
  ungroup()

# Split by year -----------------------------------------------------------

allYear <- all %>%
  group_split(year)


tarYear <- tar %>%
  group_split(year)

pipYear <- pip %>%
  group_split(year)




# Extrapolate Starting and Ending Values ----------------------------------



# Find start and end week range for positive results for each year
# Add these weeks as observations (where necessary) to extrapolate
# trailing / leading values before substituting 0's



# For test data
setSampleRange2 <- function(tibble) {

  # Subset table to first weeks with non 0 results
  nonZero <- tibble %>%
    filter(MIR > 0)

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
      MIR = 0,
      P = 0,
      posTestFreq = 0
    )


  tibble %>%
    bind_rows(joinTable) %>%
    fill(year, .direction = "updown") %>%
    arrange(clust, year, woy)




}





# Map over testing data
allExt <- map(allYear, setSampleRange2) %>%
  bind_rows() %>%
  group_split(clust, year)

pipExt <- map(pipYear, setSampleRange2) %>%
  bind_rows() %>%
  group_split(clust, year)


tarExt <- map(tarYear, setSampleRange2) %>%
  bind_rows() %>%
  group_split(clust, year)


# Split into lists --------------------------------------------------------




glimpse(test)



# Split testing by cluster and year
allSplit <- all %>%
  group_split(clust, year)

pipSplit <- pip %>%
  group_split(clust, year)

tarSplit <- tar %>%
  group_split(clust, year)


# Interpolate -------------------------------------------------------------


# Create vector of every woy between start and end sample dates
minMaxFunc <- function(table) {

  min <- min(table$woy)

  max <- max(table$woy)

  seq(min, max, 1)
}



# Testing Interpolation ---------------------------------------------------



testFunc <- function(table) {

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
    joinTableRaw <- table %>%
      full_join(weeks, by = "woy") %>%
      arrange(woy) %>%
      fill(clust) %>%
      fill(year) %>%
      mutate(
        rowNum = row_number()
      )


    joinTable <- joinTableRaw %>%
      mutate(

        P = case_when(
          (rowNum == 1 & is.na(P)) ~ 0,
          (rowNum == max(joinTableRaw$rowNum) & is.na(P)) ~ 0,
          T ~ P
        )

      ) %>%
      dplyr::select(-rowNum)



    #return(joinTable)


    # Linear Interpolations of Missing Data
    linTable <- joinTable %>%
      mutate(
        MIRLinear = na.approx(MIR),
        PIRLinear = na.approx(P),
        posTestFreqLinear = na.approx(posTestFreq)
      )


    #return(linTable)


    # Create Spline functions
    MIRSplineFun <- splinefun(table$woy,
                              table$MIR,
                              method = "monoH.FC")



    ## I had to redefine this function to deal with errors
    ## Issue with allExtF in index 35
    PIRSplineFun <- tryCatch(
      expr = {
        splinefun(table$woy,
                  table$P,
                  method = "monoH.FC")
      },
      error = function(e) {
        return(function(woy){
          NA
        })
      },
      finally = {

      }
    )


    posSplineFun <- splinefun(table$woy,
                              table$posTestFreq,
                              method = "monoH.FC")

    # Spline interpolations of Missing Data
    totalFill <- linTable %>%
      mutate(
        MIRSpline = MIRSplineFun(woy), # %>% round(digits = 1),
        PIRSpline = PIRSplineFun(woy),
        posTestFreqSpline = posSplineFun(woy)

      )

    return(totalFill)

  } else {


    # If only one year copy data over to interpolated columns
    clean <- table %>%
      mutate(
        MIRLinear = MIR,
        MIRSpline = MIR,
        PIRLinear = P,
        PIRSpline = P,
        posTestFreqLinear = posTestFreq,
        posTestFreqSpline = posTestFreq
      )



    return(clean)

  }

}


allExtF <- map(allExt, testFunc) %>%
  bind_rows() %>%
  unite(col = "wy", c(year, woy), remove = FALSE) %>%
  rename_with(~ str_c(.x, "All"), starts_with(c("MIR", "PIR", "posTest")))



# pipExtF <- map(pipExt, testFunc) %>%
#   bind_rows()  %>%
#   unite(col = "wy", c(year, woy), remove = FALSE) %>%
#   rename_with(~ str_c(.x, "Pip"), starts_with(c("MIR", "PIR", "posTest")))

pipExtF <- map(pipExt, testFunc) %>%
  bind_rows() %>%
  unite(col = "wy", c(year, woy), remove = FALSE) %>%
  rename_with(~ str_c(.x, "pip"), starts_with(c("MIR", "PIR", "posTest")))

tarExtF <- map(tarExt, testFunc) %>%
  bind_rows() %>%
  unite(col = "wy", c(year, woy), remove = FALSE) %>%
  rename_with(~ str_c(.x, "Tar"), starts_with(c("MIR", "PIR", "posTest")))

# aegExtF <- map(aegExt, testFunc) %>%
#   bind_rows() %>%
#   unite(col = "wy", c(year, woy), remove = FALSE) %>%
#   rename_with(~ str_c(.x, "aeg"), starts_with(c("MIR", "PIR", "posTest")))



extPre <- allExtF %>%
  # left_join(pipExtF, by = c("clust", "wy", "woy", "year"), suffix = c("", ".pip")) %>%
  left_join(tarExtF, by = c("clust", "wy", "woy", 'year'), suffix = c("", ".tar")) %>%
  left_join(pipExtF, by = c("clust", "wy", "woy", 'year'), suffix = c("", ".pip")) %>%
  # left_join(aegExtF, by = c("clust", "wy", "woy", "year"), suffix = c("", ".aeg")) %>%
  dplyr::select(-ends_with(".x")) %>%
  dplyr::select(-ends_with(".y")) #%>%
#dplyr::select(-PoolSize:-NumPools)






glimpse(extPre)
# Infer 0's for rest of year ----------------------------------------------


# Make grid of all possible dates as well as clusters to fill 0's in
# extrapolated data

# Vector of all clusters

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


testWClust <- allWeeks %>%
  unite(col = "wy", c(yrs, wks)) %>%
  expand_grid(clustsTest)


testWeeks <- extPre %>%
  unite(col = "wy", c(year, woy)) %>%
  full_join(testWClust, by = c("wy", "clust" = "clustsTest")) %>%
  separate(wy, c("year", "woy"), remove = FALSE) %>%
  mutate(
    year = as.numeric(year),
    woy = as.numeric(woy)
  ) %>%
  dplyr::select(-ends_with(".x")) %>%
  dplyr::select(-ends_with(".y")) %>%
  dplyr::select(clust, year, woy, everything()) %>%
  mutate(
    across(starts_with(c("MIR", "PIR", "posTest")), ~replace_na(., 0))
  )


# Export ------------------------------------------------------------------


#write_csv(intFinal, "./Data/6_interpolateData/wnvMIRPIRInt.csv")

write_csv(testWeeks, "Data/SJV/wnv/6_interpolateData/wnvMIRPIRExt.csv")



# Appendix ----------------------------------------------------------------


# # Check posTestFreq against old data --------------------------------------
#
# oldData <- read_csv("./Data/old/testClusterInterpolationInferNA.csv")
#
#
#
# aV <- allExtF %>%
#   pull(clust) %>%
#   unique()
#
#
#
# oV <- oldData %>%
#   pull(clust) %>%
#   unique()
#
#
# setdiff(oV, aV)
#
#
# oldClean <- oldData %>%
#   dplyr::select(clust, year, woy, posTestFreq)
#
# allClean <- allExtF %>%
#   dplyr::select(clust, year, woy, posTestFreq)
#
#
# joined <- oldClean %>%
#   left_join(allClean, by = c("clust", "year", "woy"), suffix = c(".old", ".all")) %>%
#   mutate(
#     diff = posTestFreq.old - posTestFreq.all
#   )
#
#
# trouble <- joined %>%
#   filter(!is.na(diff)) %>%
#   filter(abs(diff) > 0.0001) %>%
#   distinct(clust, year, woy)
#
#
# write_csv(trouble, "./Data/positiveWNVIssues.csv")
#   view()
#   summary()
#
#
#
#




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



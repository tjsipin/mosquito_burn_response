library(tidyverse)
library(lubridate)

#options(scipen = 999)

dir.create("Data/SJV/wnv/9_joinGEE")

# Read in Data ------------------------------------------------------------


interpolatedPaths <- list.files("Data/SJV/wnv/6_interpolateData",
                                full.names = T) %>%
  grep(pattern = "_NA.csv", invert = T, value = T) %>%
  enframe(name = NULL, value = "data")


geeClean <- list.files("Data/SJV/extractions", pattern = "Clean", full.names = T) %>%
  enframe(name = NULL, value = "extract")


allCombo <- interpolatedPaths %>%
  expand_grid(geeClean)


lagWeek  <-  "0"


# Join --------------------------------------------------------------------


# Function to join GEE extractions to calSurv data
lagGEEJoin <- function(data, extract) {

  # Simple subtraction of date for GEE

  lagNum <- as.numeric(lagWeek)
  lagWeeks <- weeks(lagNum)


  datTable <- read_csv(data)


  fileNameTemplate <- data %>%
    str_replace(replacement = "9_joinGEE", pattern =  "6_interpolateData")

  #return(fileNameTemplate)

  geeLag <- read_csv(extract) %>%
    mutate(
      Date = Date - lagWeeks,
      woy = week(Date),
      year = year(Date)
    ) %>%
    unite("wy", c(year, woy), remove = FALSE)


  joined <- datTable %>%
    left_join(geeLag, by = c("clust", "wy", "woy", "year"))



  # bufferDist <- extract %>%
  #   str_extract(pattern = "[:digit:]{4}")
  bufferDist <- "1500" # 1500 for now, may need to try 2000 and 3000


  writeName2 <- fileNameTemplate %>%
    str_replace("Ext", str_c(bufferDist, "LagWeeks", lagWeek))


  #return(writeName2)

  write_csv(joined, writeName2)



}



pmap(allCombo, lagGEEJoin)


list.files("Data/SJV/wnv/9_joinGEE/")


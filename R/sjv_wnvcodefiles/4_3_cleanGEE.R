library(tidyverse)
library(lubridate)

buffers <- list("1500", "2000", "3000")

# no buffers for now: only 1500
buffers <- list("")

# Read in Data ------------------------------------------------------------


cleanFunction <- function(buffer) {

  chirpsReadName <- str_c("Data/SJV/extractions/cluster", buffer, "CHIRPS.csv")


  #Floor date down to julien weeks
  chirps <- read_csv(chirpsReadName) %>%
    mutate(
      Date = as_datetime(Date) %>% date() - ddays(3),
      chirpsMean = mean,
      .keep = "unused"
    )



  grid1ReadName <- str_c("Data/SJV/extractions/cluster", buffer, "MeanGrid.csv")


  gridmet1 <- read_csv(grid1ReadName) %>%
    mutate(
      Date = as_datetime(Date) %>% date - ddays(3),
      .keep = "unused"
    )


  grid2ReadName <- str_c("Data/SJV/extractions/cluster", buffer, "SumGrid.csv")


  gridmet2 <- read_csv(grid2ReadName) %>%
    mutate(
      Date = as_datetime(Date) %>% date - ddays(3),
      .keep = "unused"
    ) %>%
    dplyr::select(Date, clust, mean) %>%
    mutate(
      `Precip Mean mm / day` = mean,
      .keep = "unused"
    )


  gridmet <- gridmet1 %>%
    left_join(gridmet2, by = c("Date", "clust"))



  modisReadName <- str_c("Data/SJV/extractions/cluster", buffer, "MODIS.csv")


  #Convert seconds since epoch to date
  modis <- read_csv(modisReadName) %>%
    mutate(
      modDate = as_datetime(Date1/1000),
      .keep = "unused"
    )




  jrcReadName <- str_c("Data/SJV/extractions/cluster", buffer, "JRC.csv")


  #Format month/year and year as first day of time periods
  jrc <- read_csv(jrcReadName) %>%
    mutate(
      jrcDate = str_replace(Date, "_", "-") %>% str_c("-01") %>% as_datetime(),
      jrcStandingWater = sum,
      .keep = "unused"
    )


  irrReadName <- str_c("Data/SJV/extractions/cluster", buffer, "IRR.csv")


  # there is an issue with the irrMapper IC in GEE, this is a placeholder for now [08/16/2023]
  irr <- read_csv(irrReadName) %>%
    filter(Date %>% startsWith("CA")) %>%
    mutate(
      Year = str_replace(Date, '.+_(.+)', '\\1'),
      irrDate = str_c(Year, "-01-01") %>% as_date(),
      irrWater = sum,
      .keep = "unused"
    )




  droughtReadName <- str_c("Data/SJV/extractions/cluster", buffer, "Drought.csv")

  drought <- read_csv(droughtReadName) %>%
    mutate(
      droDate = as_date(Date),
      .keep = "unused"
    )

  #
  # # Ebird
  #
  # ebird <- read_csv("Data/SJV/extractions/clusterEbird.csv")

  # Bird competence
  competenceReadName <- str_c("Data/SJV/extractions/cluster", buffer, "BirdCompetence.csv")

  competence <- read_csv(competenceReadName) %>%
    mutate(
      compDate = ymd(paste(year, "01", "01", sep = "-")) + weeks(week - 1),
      .keep = "unused"
    ) %>%
    rename(
      competence = mean
    )


  # Pesticide
  ## Not cluster-specific; apply to each cluster-month-year
  # pesticide <- read.csv("Data/SJV/extractions/Kern_Mos_Public_Health_PUR_2010_2021.csv") %>%
  #   dplyr::select(POUNDS_PRODUCT_APPLIED, POUNDS_CHEMICAL_APPLIED, DATE_Month, Mos_Life_Stage_Target, n) %>%
  #   mutate(pestDate = as_date(DATE_Month)) %>%
  #   dplyr::select(-DATE_Month) %>%
  #   group_by(Mos_Life_Stage_Target, pestDate) %>%
  #   dplyr::summarise(
  #     POUNDS_PRODUCT_APPLIED = sum(POUNDS_PRODUCT_APPLIED, na.rm = T),
  #     # POUNDS_CHEMICAL_APPLIED = sum(POUNDS_CHEMICAL_APPLIED, na.rm = T),
  #     .groups = "keep"
  #   ) %>%
  #   ungroup() %>%
  #   pivot_wider(id_cols = pestDate, names_from = Mos_Life_Stage_Target, values_from = POUNDS_PRODUCT_APPLIED) %>%
  #   dplyr::rename(
  #     Adulticide.Product.Lb = Adulticide,
  #     Larvicide.Product.Lb = Larvicide
  #   )




  # Create Join Key Table of Dates ------------------------------------------

  #Join datasets with matching dates
  joinEasy <- chirps %>%
    left_join(gridmet, by = c("clust", "Date"))

  #Get column of unique dates
  dayDates <- joinEasy %>%
    dplyr::select(Date) %>%
    distinct()




  #Get dates from non matching datasets
  modDates <- modis %>%
    dplyr::select(modDate) %>%
    distinct()

  jrcDates <- jrc %>%
    dplyr::select(jrcDate) %>%
    distinct()

  irrDates <- irr %>%
    dplyr::select(irrDate) %>%
    distinct()

  droDates <- drought %>%
    dplyr::select(droDate) %>%
    distinct()

  compDates <- competence %>%
    dplyr::select(compDate) %>%
    distinct()

  # pestDates <- pesticide %>%
  #   dplyr::select(pestDate) %>%
  #   distinct()

  #Create vector of every day in time period, match Dates to
  #Every day and fill down through time period
  dateKey <- seq(as.Date("2009-12-03"), as.Date("2023-12-31"), by = "days") %>%
    as_tibble() %>%
    left_join(modDates, by = c("value" = "modDate"), keep = TRUE) %>%
    left_join(jrcDates, by = c("value" = "jrcDate"), keep = TRUE) %>%
    left_join(irrDates, by = c("value" = "irrDate"), keep = TRUE) %>%
    left_join(droDates, by = c("value" = "droDate"), keep = TRUE) %>%
    left_join(compDates, by = c("value" = "compDate"), keep = TRUE) %>%
    # left_join(pestDates, by = c("value" = "pestDate"), keep = TRUE) %>%
    mutate(year = year(value)) %>%
    group_by(year) %>%
    fill(modDate, jrcDate, irrDate, droDate, compDate) %>%
    ungroup() %>%
    dplyr::select(-year) %>%
    #Select only the first day of week
    semi_join(dayDates, by = c("value" = "Date"))




  # Join  Tables ------------------------------------------------------------


  totalJoin <- joinEasy %>%
    left_join(dateKey, by = c("Date" = "value")) %>%
    left_join(modis, by = c("clust", "modDate")) %>%
    left_join(jrc, by = c("clust", "jrcDate")) %>%
    left_join(irr, by = c("clust", "irrDate")) %>%
    left_join(drought, by = c('clust', 'droDate')) %>%
    left_join(competence, by = c('clust', 'compDate')) %>%
    # left_join(pesticide, by = c('pestDate')) %>%
    dplyr::select(-starts_with("system")) %>%
    dplyr::select(-starts_with(".geo")) %>%
    dplyr::select(-matches("...Date"))


  finalWriteName <- str_c("Data/SJV/extractions/extract", buffer, "Clean.csv")


  write_csv(totalJoin, finalWriteName)

}

map(buffers, cleanFunction)

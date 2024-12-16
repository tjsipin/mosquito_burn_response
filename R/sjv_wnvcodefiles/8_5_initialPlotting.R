library(tidyverse)
library(lubridate)

#options(scipen = 999)



dir.create("Data/8_initialPlotting")


# Pipiens -----------------------------------------------------------------

all <- read_csv("Data/6_interpolateData/allExt.csv") %>%
  mutate(
    data = "all"
  )

pipAll <- read_csv("Data/6_interpolateData/pipiensAllExt.csv")  %>%
  mutate(
    data = "pipA"
  )


pipF <- read_csv("Data/6_interpolateData/pipiensFemaleExt.csv")  %>%
  mutate(
    data = "pipF"
  )

tarAll <- read_csv("Data/6_interpolateData/tarsalisAllExt.csv")  %>%
  mutate(
    data = "tarA"
  )

tarF <- read_csv("Data/6_interpolateData/tarsalisFemaleExt.csv")  %>%
  mutate(
    data = "tarF"
  )


dataV <- c("all", "tarA", "tarF")


# Combine Datasets --------------------------------------------------------

total <- bind_rows(all, tarAll, tarF)


# 2018 Mosquitoes / Trap by Cluster
plot2018 <- function(dataID) {

  chartTitle <- str_c("2018 ", dataID, " mosquitoes by Cluster")


  total %>%
    filter(data == dataID) %>%
    filter(year == 2018) %>%
    mutate(
      clust = as.factor(clust)
    ) %>%
    ggplot(mapping = aes(x = woy, y = mosPerTrapNight, col = clust)) +
    geom_line() +
    #facet_grid(cols = vars(year)) +
    theme(legend.position = "none") +
    ylab("Mosquitoes / Trap Night") +
    xlab("Week of Year") +
    labs(title = chartTitle)

  exportTitle <- str_c("Data/8_initialPlotting/2018abund", dataID, ".jpeg")

  ggsave(exportTitle)
}


#plot2010(dataV[[1]])

map(dataV, plot2010)



# 2018, 2020, 2022 Mosquitoes / Trap by Cluster
plot3Year <- function(dataID) {


  chartTitle <- str_c("2018, 2020, 2022 ", dataID, " mosquitoes by Cluster")


  total %>%
    filter(data == dataID) %>%
    filter(year %in% c(2018, 2020, 2022)) %>%
    mutate(
      clust = as.factor(clust)
    ) %>%
    ggplot(mapping = aes(x = woy, y = mosPerTrapNight, col = clust)) +
    geom_line() +
    facet_grid(cols = vars(year)) +
    theme(legend.position = "none") +
    ylab("Mosquitoes / Trap Night") +
    xlab("Week of Year") +
    labs(title = chartTitle)



  exportTitle <- str_c("Data/8_initialPlotting/2018_2020_2022_abund", dataID, ".jpeg")

  ggsave(exportTitle)
}


map(dataV, plot3Year)



plotAllYears <- function(dataID) {


  chartTitle <- str_c("2018 - 2023 ", dataID, " mosquitoes by Cluster")


  # Mos / Trap by Cluster all years
  total %>%
    filter(data == dataID) %>%
    #filter(year == 2010) %>%
    mutate(
      clust = as.factor(clust),
      collection_date = ymd(str_c(year, "-01-01")) + weeks(woy - 1)
    ) %>%

    ggplot(mapping = aes(x = collection_date, y = mosPerTrapNight, col = clust)) +
    geom_line() +
    theme(legend.position = "none") +
    ylab("Mosquitoes / Trap Night") +
    xlab("Collection Date") +
    labs(title = chartTitle) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y")



  exportTitle <- str_c("Data/8_initialPlotting/2018_2023_abund", dataID, ".jpeg")

  ggsave(exportTitle)

}


map(dataV, plotAllYears)


# Total Mosquitoes --------------------------------------------------------





# 2018 Mosquitoes by Cluster
all %>%
  filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = totalMosquitoes, col = clust)) +
  geom_line() +
  #facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018 Total Mosquitoes by Cluster")


# 2018, 2020, 2022 Mosquitoes by Cluster
all %>%
  filter(year %in% c(2018, 2020, 2022)) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = totalMosquitoes, col = clust)) +
  geom_line() +
  facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018, 2020, 2023 Total Mosquitoes by Cluster")


# Total Mosquitoes by Cluster all years
all %>%
  #filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = year, y = totalMosquitoes, col = clust)) +
  geom_jitter() +
  theme(legend.position = "none") +
  labs(title = "2018 - 2020 Total Mosquitoes by Cluster") #+
  #scale_x_date(date_breaks = "1 year", date_labels = "%Y")






# Infection Rate ----------------------------------------------------------





# 2018 Infection Rate by Cluster
testExt %>%
  filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = MIRAll, col = clust)) +
  geom_jitter() +
  #facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018 WNV Infection Rate by Cluster")


# 2018, 2020, 2023 Infection Rate by Cluster
testExt %>%
  filter(year %in% c(2018, 2020, 2023)) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = MIRAll, col = clust)) +
  geom_jitter() +
  facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018, 2015, 2020 WNV Infection Rate by Cluster")


# Total Infection Rate by Cluster all years
testExt %>%
  #filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%

  ggplot(mapping = aes(x = year, y = MIRAll, col = clust)) +
  geom_jitter() +
  theme(legend.position = "none") +
  labs(title = "2018 - 2023 WNV Infection Rate by Cluster")






# Positive Mosquitoes per Test --------------------------------------------




# 2018 Positive Mosquitoes / Test by Cluster
testExt %>%
  filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = Pos, col = clust)) +
  geom_line() +
  #facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018 WNV Positive Mosquitoes / Test by Cluster")


# 2018, 2020, 2023 Positive Mosquitoes / Test by Cluster
testExt %>%
  filter(year %in% c(2018, 2020, 2023)) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = Pos, col = clust)) +
  geom_line() +
  facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018, 2020, 2023 WNV Positive Mosquitoes / Test by Cluster")


# Total Positive Mosquitoes / Test by Cluster all years
testExt %>%
  #filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%

  ggplot(mapping = aes(x = collection_date, y = Pos, col = clust)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title = "2018 - 2020 WNV Positive Mosquitoes / Test by Cluster") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")





# Positive Test Frequency per Test --------------------------------------------


# Raw

# 2018 Postive Test Frequency by Cluster
testExt %>%
  filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = posTestFreqAll, col = clust)) +
  geom_jitter() +
  #facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018 WNV Positive Test Frequency by Cluster")


# 2018, 2020, 2023 Positive Test Frequency by Cluster
testExt %>%
  filter(year %in% c(2018, 2020, 2023)) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = posTestFreqAll, col = clust)) +
  geom_jitter() +
  facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018, 2020, 2023 WNV Positive Test Frequency by Cluster")


# Total Positive Test Frequency by Cluster all years
testExt %>%
  #filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%

  ggplot(mapping = aes(x = collection_date, y = posTestFreqAll, col = clust)) +
  geom_jitter() +
  theme(legend.position = "none") +
  labs(title = "2018 - 2023 WNV Positive Test Frequency by Cluster") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")





# Spline Extrapolation

# 2018 Positive Test Frequency by Cluster
testExt %>%
  filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = posTestFreqSplineAll, col = clust)) +
  geom_line() +
  #facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018 WNV Positive Test Frequency by Cluster Spline Ext")


# 2018, 2020, 2023 Positive Test Frequency by Cluster
testExt %>%
  filter(year %in% c(2018, 2020, 2023)) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%
  ggplot(mapping = aes(x = woy, y = posTestFreqSplineAll, col = clust)) +
  geom_line() +
  facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018, 2020, 2023 WNV Positive Test Frequency by Cluster Spline Ext")


# Total Positive Test Frequency by Cluster all years
testExt %>%
  #filter(year == 2018) %>%
  mutate(
    clust = as.factor(clust)
  ) %>%

  ggplot(mapping = aes(x = collection_date, y = posTestFreqSplineAll, col = clust)) +
  geom_point() +
  theme(legend.position = "none") +
  labs(title = "2018 - 2023 WNV Positive Test Frequency by Cluster Spline Ext") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")





# Total Positive Tests ----------------------------------------------------


totalPos <- testExt %>%
  group_by(year, woy) %>%
  summarize(
    testPositive = sum(Pos, na.rm = TRUE),
    collection_date = first(collection_date)
  )


# 2018 Positive Test Frequency by Cluster
totalPos %>%
  filter(year == 2018) %>%
  # mutate(
  #   clust = as.factor(clust)
  # ) %>%
  ggplot(mapping = aes(x = woy,
                       y = testPositive,
                       #col = clust
                       )) +
  geom_line() +
  #facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018 Total Positive Tests")


# 2018, 2020, 2023 Positive Test Frequency by Cluster
totalPos %>%
  filter(year %in% c(2018, 2020, 2023)) %>%
  # mutate(
  #   clust = as.factor(clust)
  # ) %>%
  ggplot(mapping = aes(x = woy,
                       y = testPositive,
                       #col = clust
                       )) +
  geom_line() +
  facet_grid(cols = vars(year)) +
  theme(legend.position = "none") +
  labs(title = "2018, 2020, 2023 Total Positive Tests")


# Total Positive Test Frequency by Cluster all years
totalPos %>%
  #filter(year == 2018) %>%
  # mutate(
  #   clust = as.factor(clust)
  # ) %>%

  ggplot(mapping = aes(x = collection_date,
                       y = testPositive,
                       #col = clust
                       )) +
  geom_line() +
  theme(legend.position = "none") +
  labs(title = "2018 - 2023 Total Positive Tests") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")







# Join to GEE extractions -------------------------------------------------

lag <- weeks(4)


geeLag <- read_csv("Data/extractions/extractClean.csv") %>%
  mutate(
    Date = Date - lag,
    woy = week(Date),
    year = year(Date)
  ) %>%
  unite("wy", c(year, woy), remove = FALSE)


gee <- read_csv("Data/extractions/extractClean.csv") %>%
  mutate(
    woy = week(Date),
    year = year(Date)
  ) %>%
  unite("wy", c(year, woy), remove = FALSE)


glimpse(gee)



joinList <- map(dataList, left_join, gee, by = c("clust", "wy", "woy", "year"))

map(joinList, nrow)


joinListLag <- map(dataList, inner_join, geeLag, by = c("clust", "wy", "woy", "year"))

map(joinListLag, nrow)

# Reshape data ------------------------------------------------------------

joinList[[2]] %>%
  select(clust, wy, woy, year, everything()) %>%
  select(-Date) %>%
  select(-eddi14d:-z) %>%

  group_by(year, woy) %>%
  summarize(
    across(totalMosquitoes:mosPerSpline, ~ sum(., na.rm = TRUE)),
    across(chirpsMean:irrWater, ~ mean(., na.rm = TRUE))
  ) %>%
  pivot_longer(chirpsMean:irrWater, names_to = "x", values_to = "xVal") %>%
  pivot_longer(totalMosquitoes:mosPerSpline, names_to = "y", values_to = "yVal") %>%
  filter(!is.na(yVal) ) %>%
  #filter(year == 2018) %>%
  filter(y == "mosPerSpline") %>%



  ggplot(mapping = aes(x = xVal, y = yVal)) +
  geom_point() +
  facet_grid(rows = vars(year), cols = vars(x), scales = "free") +
  labs(y = "Mosquito per Trap Night")


# joinList[[1]] %>%
#   select(clust, wy, woy, year, everything()) %>%
#   select(-Date) %>%
#   select(-eddi14d:-z) %>%
#
#   group_by(year, woy) %>%
#   summarize(
#     across(totalMosquitoes:mosPerTrapNight, ~ sum(., na.rm = TRUE)),
#     across(chirpsMean:irrWater, ~ mean(., na.rm = TRUE))
#   ) %>%
#   pivot_longer(chirpsMean:irrWater, names_to = "x", values_to = "xVal") %>%
#   pivot_longer(totalMosquitoes:mosPerTrapNight, names_to = "y", values_to = "yVal") %>%
#   filter(!is.na(yVal) ) %>%
#   #filter(year == 2018) %>%
#   filter(y == "mosPerTrapNight") %>%
#
#
#
#   ggplot(mapping = aes(x = xVal, y = yVal)) +
#   geom_point() +
#   facet_grid(rows = vars(year), cols = vars(x), scales = "free") +
#   labs(y = "Mosquito per Trap Night")

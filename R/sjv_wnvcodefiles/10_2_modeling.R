library(tidyverse)
library(lubridate)
library(modelr)

#options(scipen = 999)




# Abundance Data ----------------------------------------------------------

abundRaw <- read_csv("Data/9_joinGEE/all1500LagWeeks0.csv")



abund <- abundRaw %>%
  mutate(
    tempMaxC = `Temp Max K` - 273.15,
    tempMinC = `Temp Min K` - 273.15,
    clust = as_factor(clust),
    .keep = "unused"
  )


glimpse(abund)




clusterWOY <- abund %>%
  group_by(woy, clust) %>%
  summarize(
    meanMos = mean(mosPerTrapNight, na.rm = TRUE),
    medianMos = median(mosPerTrapNight, na.rm = TRUE),

    meanMosLin = mean(mosPerLinear, na.rm = TRUE),
    medianMosLin = median(mosPerLinear, na.rm = TRUE),


    meanMaxT = mean(tempMaxC, na.rm = TRUE),
    #medianMaxT = median(tempMaxC, na.rm = TRUE)
  )



woyY <- abund %>%
  group_by(woy, year) %>%
  summarize(
    meanMos = mean(mosPerTrapNight, na.rm = TRUE),
    medianMos = median(mosPerTrapNight, na.rm = TRUE),

    meanMosLin = mean(mosPerLinear, na.rm = TRUE),
    medianMosLin = median(mosPerLinear, na.rm = TRUE),


    meanMaxT = mean(tempMaxC, na.rm = TRUE),
    #medianMaxT = median(tempMaxC, na.rm = TRUE)
  )


clusterWOY %>%
  # mutate(meanMosLin = log10(1 + meanMos)) %>%
  mutate(meanMosLin = meanMos) %>%
  ggplot(mapping = aes(meanMaxT, meanMosLin)) +
  geom_point(aes(color = clust)) +
  theme(legend.position = "none") +
  labs(title = "Max Temp vs Adjusted Abundance by Cluster, WOY") +
  xlab("Maximum Temp (C)") +
  ylab("Mosquitoes per Trap Night") +
  scale_x_continuous(breaks = seq(0, 50, 5))


clusterWOY %>%
  mutate(meanMosLin = log(1 + meanMos)) %>%
  #mutate(meanMosLin = meanMos) %>%
  ggplot(mapping = aes(meanMaxT, meanMosLin)) +
  geom_point(aes(color = clust)) +
  theme(legend.position = "none") +
  labs(title = "Max Temp vs Adjusted Abundance by Cluster, WOY") +
  xlab("Maximum Temp (C)") +
  ylab("Mosquitoes per Trap Night ln") +
  scale_x_continuous(breaks = seq(0, 50, 5))


abund %>%
  mutate(mosLin = log10(1 + mosPerTrapNight)) %>%
  ggplot(mapping = aes(tempMaxC, mosLin)) +
  geom_point(aes(color = clust)) +
  theme(legend.position = "none") +
  labs(title = "Max Temp vs Adjusted Abundance Log10") +
  xlab("Maximum Temp (C)") +
  ylab("Mosquitoes per Trap Night Log10") +
  scale_x_continuous(breaks = seq(0, 50, 5))



clusterWOY %>%
  ggplot(mapping = aes(meanMaxT, meanMos)) +
  geom_point(aes(color = clust)) +
  theme(legend.position = "none") +
  labs(title = "Max Temp vs Adjusted Abundance by Cluster, Year")

woyY %>%
  ggplot(mapping = aes(meanMaxT, meanMos)) +
  geom_point(aes(color = as.factor(year))) +
  theme(legend.position = "none") +
  #scale_y_continuous(limits = c(0, 40)) +
  labs(title = "Max Temp vs Adjusted Abundance by  WOY, Year")


newModel <- lm(mosPerTrapNight ~ tempMaxC,
               data = abund %>%
                 group_by(clust)
               )


newModel <- lm(mosPerTrapNight ~ tempMaxC + factor(clust),
               data = abund)



summary(newModel)

abundModel <- abund %>%
  add_predictions(newModel) %>%
  add_residuals(newModel) %>%
  select(-eddi14d:-z)


glimpse(abundModel)

abundModel %>%
  ggplot(aes(tempMaxC)) +
  geom_point(aes(y = mosPerTrapNight)) +
  geom_line(aes(x = tempMaxC, y = pred), data = abundModel, colour = "red")




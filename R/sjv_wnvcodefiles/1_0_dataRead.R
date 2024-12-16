library(tidyverse)
library(backports)


# Read in Abundance Data --------------------------------------------------



readAbund <- function(string){

  #Import data
  data <- suppressWarnings(read_csv(string))


  # #Save table of column reading problems
  # prob <- problems()
  #
  # #View by columns
  # prob %>%
  #   group_by(col, expected) %>%
  #   distinct(actual) %>%
  #   arrange(col) %>%
  #   view()
  #
  # #Get list of problem columns
  # probCol <- distinct(prob["col"]) %>%
  #   arrange(col) %>%
  #   pull(col)
  #
  #
  # data %>%
  #   select(probCol) %>%
  #   view()


  #Get data classes as tibble
  orig <- sapply(data, class) %>% tibble()

  #orig %>% view()

  #Build conversion matrix, for new import
  conversionMatrix <-
    tribble(
      ~old, ~new,
      "character", "c",
      "numeric", "n",
      "logical", "l",
      "Date", "D"
    )

  #Join conversion matrix to class matrix
  colVectRaw <- orig %>%
    left_join(conversionMatrix, by = c("." = "old")) %>%
    pull(new)




  #Create extra variable
  colVect <- colVectRaw
  #Modify problem columns
  colVect[c(8, 14)] <- "c"
  #colVect[23:24] <- "d"
  colVect[36:180] <- "n"
  #Collapse vector of character to single string
  stringCol <- paste(colVect, collapse =  "")
  #Read in with correct column types
  abund <- read_csv(string, col_types = stringCol)

  #rm(conversionMatrix, colVect, colVectRaw, data, orig, stringCol)

  #Return table and replace NA's with string equivalent
  abund %>%
    replace_na(list(name = "NA"))
}


#Test function
# abund <- readAbund("Data/KMVCD/Raw/Mosquito_Pools_2018_2023.csv")

# Read in Pools Data ------------------------------------------------------

readPools <- function(string){
  #Import data
  data <- read_csv(string)


  # #Save table of column reading problems
  # prob <- problems()
  #
  # #View by columns
  # prob %>%
  #   group_by(col, expected) %>%
  #   distinct(actual) %>%
  #   arrange(col) %>%
  #   view()
  #
  # #Get list of problem columns
  # probCol <- distinct(prob["col"]) %>%
  #   arrange(col) %>%
  #   pull(col)
  #
  #
  # data %>%
  #   select(probCol) %>%
  #   view()
  #
  #
  # str(data)

  #Get data classes as tibble
  orig <- sapply(data, class) %>% tibble()

  #orig %>% view()

  #Build conversion matrix, for new import
  conversionMatrix <-
    tribble(
      ~old, ~new,
      "character", "c",
      "numeric", "n",
      "logical", "l",
      "Date", "D"
    )

  #Join conversion matrix to class matrix
  colVectRaw <- orig %>%
    left_join(conversionMatrix, by = c("." = "old")) %>%
    pull(new)

  #Create extra variable
  colVect <- colVectRaw
  #Modify problem columns
  colVect[12] <- "c"
  colVect[18] <- "c"

  #Collapse vector of character to single string
  stringCol <- paste(colVect, collapse =  "")
  #Read in with correct column types
  pools <- read_csv(string , col_types = stringCol)

  #rm(conversionMatrix, colVect, colVectRaw, data, orig, stringCol)

  pools %>%
    replace_na(list(name = "NA"))
}



# Read in Testing Data ----------------------------------------------------


readTest <- function(string) {
  #Import data
  data <- read_csv(string)

  #
  # #Save table of column reading problems
  # prob <- problems()
  #
  # #View by columns
  # prob %>%
  #   group_by(col, expected) %>%
  #   distinct(actual) %>%
  #   arrange(col) %>%
  #   view()
  #
  # #Get list of problem columns
  # probCol <- distinct(prob["col"]) %>%
  #   arrange(col) %>%
  #   pull(col)
  #
  #
  # data %>%
  #   select(probCol) %>%
  #   view()
  #
  #
  # str(data)

  #Get data classes as tibble
  orig <- sapply(data, class) %>% tibble()

  #orig %>% view()

  #Build conversion matrix, for new import
  conversionMatrix <-
    tribble(
      ~old, ~new,
      "character", "c",
      "numeric", "n",
      "logical", "l",
      "Date", "D"
    )

  #Join conversion matrix to class matrix
  colVectRaw <- orig %>%
    left_join(conversionMatrix, by = c("." = "old")) %>%
    pull(new)

  #Create extra variable
  colVect <- colVectRaw
  #Modify problem columns
  colVect[c(9, 15, 17, 32)] <- "c"
  colVect[18] <- "n"

  #Collapse vector of character to single string
  stringCol <- paste(colVect, collapse =  "")
  #Read in with correct column types

  testResults <- read_csv(string, col_types = stringCol)


  #rm(conversionMatrix, colVect, colVectRaw, data, orig, stringCol)
  testResults %>%
    replace_na(list(name = "NA"))
}

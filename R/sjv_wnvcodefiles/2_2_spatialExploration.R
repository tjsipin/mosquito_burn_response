
# Packages ----------------------------------------------------------------

library(tidyverse)
library(sf)
library(rgdal)
library(geosphere)

# Notes -------------------------------------------------------------------
#This is the KMVCD version
#There is no test results data, so all test results code is commented out

# Run Import Script -------------------------------------------------------

#This will throw 3 parsing error messages
#But are corrected from initial export in the script
source("R/sjv_wnvcodefiles/1_0_dataRead.R") #this is just to import data with correct column data types; doesn't work so should modify



readFunc <- function(dataset) {
  dataset %>%
    mutate(
      agency_code = str_to_lower(str_squish(agency_code)),
      lName = str_to_lower(name),
      aName = str_c(agency_code, lName),
      code = as.character(code),
      aCode = str_c(agency_code, code)
    ) %>%
    filter(latitude != 0)
}



abund <- readAbund("Data/SJV/wnv/0_input/MosquitoAbundance2010.2020.csv")  %>%
  readFunc()


pools <- readPools("Data/SJV/wnv/0_input/MosquitoPools2010.2020.csv") %>%
  readFunc()


testResults <- readTest("Data/SJV/wnv/0_input/MosquitoPoolsTestResults2010.2020.csv") %>%
  readFunc()




# Create output directory
dir.create("Data/SJV/wnv/2_spatialExploration")




# Compare aCode with aName ------------------------------------------------

#Bind codes and names from the 3 datasets
codeNameRaw <- abund %>%
  # dplyr::select(aCode, aName, name, lName, code, latitude, longitude) %>%
  # bind_rows(pools %>%
  #             dplyr::select(aCode, aName, name, lName, code, latitude, longitude)
  # )
  # ) %>%
  bind_rows(testResults %>%
              dplyr::select(aCode, aName, name, lName, code, latitude, longitude)
  )

#Find distinct code / name combinations
codeName <- codeNameRaw %>%
  distinct(aCode, aName, latitude, longitude)

#Group by code, and find observastions with multiple names
aCodeMulti <- codeName %>%
  dplyr::select(-latitude, -longitude) %>%
  group_by(aCode) %>%
  distinct() %>%
  dplyr::count(aCode) %>%
  filter(n > 1) %>%
  pull(aCode)

aCodeRows <- codeName %>%
  filter(aCode %in% aCodeMulti)

#Group by name, and find observations with multiple codes
aNameMulti <- codeName %>%
  dplyr::select(-latitude, -longitude) %>%
  group_by(aName) %>%
  distinct() %>%
  count(aName) %>%
  filter(n > 1) %>%
  pull(aName)

aNameRows <- codeName %>%
  filter(aName %in% aNameMulti)

# #Rows with no lat / long
noGeo <- codeName %>%
  filter(latitude == 0 | longitude == 0)

codeNameRaw %>%
  semi_join(noGeo, by = c("aCode", "aName"))


mismatches <- aNameRows %>%
  full_join(aCodeRows, by = c("aCode", "aName"))

#Check how many observations in data are under mismatched codes / names
codeNameRaw %>%
  semi_join(mismatches, by = c("aCode", "aName"))


# Drop Mismatched names ---------------------------------------------------

dropMismatch <- codeNameRaw %>%
  anti_join(mismatches, by = c("aCode", "aName"))


mismatchNames <- mismatches %>% pull(aName) %>% unique()

abund %>%
  semi_join(mismatches, by = c("aCode", "aName"))

#Drop Mismatched names from calSurv Data

abundDrop <- abund %>%
  anti_join(mismatches, by = c("aCode", "aName"))

poolsDrop <- pools %>%
  anti_join(mismatches, by = c("aCode", "aName"))

testDrop <- testResults %>%
  anti_join(mismatches, by = c("aCode", "aName"))




# Rowbind Pools, Testing, Abundance Name + geo ----------------------------



createTable <- function() {



  #Names, lat long for all data
  id <- dropMismatch %>%
    dplyr::select(aName, latitude, longitude)

  #Check how many duplicate geometries there are
  distinctNames <- id %>%
    distinct(aName, latitude, longitude) %>%
    group_by(aName)


  # abund %>%
  #   distinct(!!col, latitude, longitude) %>%
  #   group_by(!!col)

  distinctNames

}



# Find Distance between mismatched points ---------------------------------

#Find distance between a centroid for grouped points with
#multiple lat/longs
dstFunction <- function(colName) {

  #Enquo column name to use as argument
  col <- enquo(colName)

  #Create table of lat longs for each column id
  colTable <- createTable()


  #Subset data with to colID with multiple lat/longs per id
  subsetID <- colTable %>% group_by(!!col) %>%
    summarize(
      n = n_distinct(longitude),
      nlat = n_distinct(latitude),
      count = n()
    ) %>%
    filter((n > 1) | (nlat > 1)) %>%
    pull(!!col)

  #Find distinct colID + geo for data
  geo <- colTable %>%
    filter(!!col %in% subsetID) %>%
    distinct(!!col, longitude, latitude) %>%
    mutate(
      lon = longitude,
      lat = latitude
    )

  #Import as sf object
  geoSF <- st_as_sf(geo, coords = c("longitude", "latitude"))
  st_crs(geoSF) <- 4269

  #Create centroids for groups
  dataCentroids <- geo %>%
    group_by(!!col) %>%
    summarize(
      x = mean(longitude),
      y = mean(latitude),
      lat.c = y,
      lon.c = x
    ) %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(geoSF))
  #Modify CRS
  st_crs(dataCentroids) <- 4269


  #Get distance from centroid for each point
  dataDist <- geoSF


  #Function to find distance between point and centroid
  newFunc <- function(x){

    #Get row of Ungrouped points by row number x
    xPoint <- dataDist %>%
      ungroup() %>%
      slice(x)

    #Extract id
    val <- xPoint %>%
      pull(!!col)

    #Find centroid for row number x
    xCentroid <- dataCentroids %>%
      filter(
        !!col == val
      )

    #Get distance between point and centroid
    st_distance(xPoint, xCentroid)
  }

  #Map distance function over every row
  dataDist$dst <- map_dbl(1:nrow(dataDist), newFunc)

  #Convert centroid sf to tibble
  centroidsGeo <- dataCentroids %>%
    st_drop_geometry() %>%
    dplyr::select(!!col, lat.c, lon.c)

  #Add centroid lat long to original tibble
  dataDist %>%
    left_join(centroidsGeo, by = quo_name(col))


  #dataDist
  #dataDist %>% slice(1)
  #newFunc(2)
}

#
# #Run names function on all data
aNameDist <- dstFunction(aName)



#Create function to display list of metrics for distance between grouped points
dstMetrics <- function(colName) {

  col <- enquo(colName)

  dat <- dstFunction(!!col)

  noGeo <- dat %>%
    st_drop_geometry()

  f <- c("mean", "median", "max", "min", "sd")

  params <- list(
    list(dat$dst)
  )

  firstMap <- invoke_map(f, params)

  setNames(firstMap, f)


}

#Run metrics function on column variables

dstMetrics(aName)




# Update Close Grouped Points -----------------------------------------------


#Check to see if there are groups with some spatial outliers
#Filter by distance to centroid
farPoints <- aNameDist %>%
  filter(dst > 1500) %>%
  pull(aName)

farPoints %>% unique()



#If in farPoints, change lat lon .c to old geometry
updatePoints <- aNameDist %>%
  mutate(
    lat.c = if_else(aName %in% farPoints, lat, lat.c),
    lon.c = if_else(aName %in% farPoints, lon, lon.c)
  ) %>%
  st_drop_geometry()




# Subset distinct names matching latlong ----------------------------------

multiName <- aNameDist %>%
  distinct(aName) %>%
  pull(aName)




uniqueGeoTable <- dropMismatch %>%
  distinct(aName, latitude, longitude) %>%
  group_by(aName) %>%
  filter(!(aName %in% multiName)) %>%
  mutate(
    lat = latitude,
    lon = longitude,
    .keep = "unused"
  ) %>%
  bind_rows(updatePoints) %>%
  dplyr::select(-dst) %>%
  mutate(
    lat.c = if_else(is.na(lat.c), lat, lat.c),
    lon.c = if_else(is.na(lon.c), lon, lon.c)
  ) %>%
  group_by(aName, lat.c, lon.c) %>%
  mutate(
    id = cur_group_id()
  )

#


# Export ID /Lat/Long CSV -------------------------------------------------



idExport <- uniqueGeoTable %>%
  ungroup %>%
  dplyr::select(id, lat.c, lon.c) %>%
  mutate(
    lat = lat.c,
    lon = lon.c,
    .keep = "unused"
  ) %>%
  distinct()

write_csv(idExport, "Data/SJV/1_DataProcessing/idExports/idPoints.csv")


idExport %>%
  filter(is.na(lat) | is.na(lon))


# Join ID to data ---------------------------------------------------------


#Create function to join id tibble to data set based on name, and original
#lat/long
joinId <- function(tibble) {
  tibble %>%
    #Inner join
    left_join(uniqueGeoTable, by = c("longitude" = "lon",
                                     "aName" = "aName",
                                     "latitude" = "lat")
    ) #%>%
  #select(aName, id, longitude, latitude, lon.c, lat.c, everything())
}



poolsID <- joinId(poolsDrop)
testID <- joinId(testDrop)
abundID <- joinId(abundDrop)


write_csv(poolsID, "Data/SJV/wnv/2_spatialExploration/poolsID.csv")
write_csv(testID, "Data/SJV/wnv/2_spatialExploration/testID.csv")
write_csv(abundID, "Data/SJV/wnv/2_spatialExploration/abundID.csv")




# Appendix ----------------------------------------------------------------
#
# # Number of points per ID -------------------------------------------------
#
#
#
# #Looking for unique identifier columns in the data
#
# #Think i decided not to use code because of the large number of points
# #There were for certain codes,
# abund %>% count(code) %>%
#   #filter(n > 10) %>%
#   {{mean(.$n)}}
#
# abund %>% count(collection_id)
#
# abund %>% count(name) %>%
#   #filter(n > 10) %>%
#   {{mean(.$n)}}
#
# abund %>% count(aName) %>%
#   {{mean(.$n)}}
#
# abund %>% count(lName) %>%
#   {{mean(.$n)}}
#
#
# #Find number of names with multiple lat longs
# abund %>% group_by(name) %>%
#   summarize(
#     n = n_distinct(longitude),
#     nlat = n_distinct(latitude),
#     count = n()
#   ) %>%
#   filter((n > 1) | (nlat > 1)) %>%
#   #filter(n != nlat)
#   {{sum(.$count)}}
#
# abund %>% group_by(aName) %>%
#   summarize(
#     n = n_distinct(longitude),
#     nlat = n_distinct(latitude),
#     count = n()
#   ) %>%
#   filter((n > 1) | (nlat > 1)) %>%
#   #filter(n != nlat)
#   {{sum(.$count)}}
#
#
# #Was thinking that code would uniquely identify a trap location
# #but address might be better
# #Name seems to be a location
#
#
# abund %>%
#   filter(is.na((latitude)))
#
#

#
#
# # Number of Unique ID's ---------------------------------------------------
#
#
#
#
# abund %>% distinct(coordinate_precision)
#
#
#
#
#
# spatialExplorationFunction <- function(colName) {
#
#   col <- enquo(colName)
#
#   #Try with names
#   abundCol <- abund %>% distinct(!!col) %>% pull()
#   poolsCol <- pools %>% distinct(!!col) %>% pull()
#   testCol <- testResults %>% distinct(!!col) %>% pull()
#
#
#   #Get number of IDs
#   abundLength <- abundCol %>% length()
#   poolsLength <- poolsCol %>% length()
#   testLength <- testCol %>% length()
#
#
#   # #Check NA's
#   # abundNA <- abund %>%
#   #   distinct(!!col) %>%
#   #   filter(is.na(!!col) | !!col == "NA") %>%
#   #   pull()
#   #
#   # poolsNA <- pools %>%
#   #   distinct(!!col) %>%
#   #   filter(is.na(!!col) | !!col == "NA") %>%
#   #   pull()
#   #
#   # testNA <- testResults %>%
#   #   distinct(!!col) %>%
#   #   filter(is.na(!!col) | !!col == "NA") %>%
#   #   pull()
#
#
#   #Number of differences between pool and test !!cols equals
#   #The difference between number of entries
#   #So !!col column matches up between these two
#   poolsTestDiff <- setdiff(poolsCol, testCol) %>% length() # 37 differences
#   testPoolsDiff <- setdiff(testCol, poolsCol) %>% length() # 0 differences
#
#
#   #Testing number of differences between abundance and pools
#   abundPoolsDiff <- setdiff(abundCol, poolsCol) %>% length() #958 differences
#   poolsAbundDiff <- setdiff(poolsCol, abundCol) %>% length() #839 differences
#
#   #Create List to output
#   display <- tibble(rlang::as_name(col),
#                     abundLength, poolsLength, testLength,
#                     #abundNA, poolsNA, testNA,
#                     poolsTestDiff, testPoolsDiff,
#                     abundPoolsDiff, poolsAbundDiff)
#
#   #Name list elements for readability
#   set_names(display,
#             c("colName",
#               "Abundance ID's", "Pool ID's", "Test ID's",
#               #"Abundance NA's", "Pool NA", "Test NA",
#               "Pools Test Difference", "Test Pools Difference",
#               "Abundance Pools Difference", "Pools Abundance Difference")
#   )
#
# }
#
#
# try1 <- spatialExplorationFunction(name)
# try2 <- spatialExplorationFunction(lName)
# try3 <- spatialExplorationFunction(aName)
# try4 <- spatialExplorationFunction(code)
# try5 <- spatialExplorationFunction(aCode)
#
# bind_rows(try1, try2, try3, try4, try5) %>%
#   view()
#
#
# abund %>%
#   filter(
#     coordinate_precision == "Exact" | coordinate_precision == "Approximate"
#   )
#

#
# uniqueGeo <- createTable()
#
# disGeoRaw <- uniqueGeo %>%
#   st_as_sf(coords = c("longitude", "latitude")) %>%
#   ungroup() %>%
#   distinct(geometry) %>%
#   mutate(
#     geoID = row_number()
#   )
#
#
# disGeoRawJoin <- uniqueGeo %>%
#   st_as_sf(coords = c("longitude", "latitude")) %>%
#   st_join(disGeoRaw)
#
#
# disGeoRawJoin %>%
#   group_by(geoID) %>%
#   mutate(
#     n = n_distinct(aName)
#   ) %>%
#   filter(n > 1) %>%
#   view()

#
# abundID %>%
#   filter(is.na(id)) %>%
#   pull(aName) %>%
#   setdiff(mismatches %>% pull(aName))
# #view()

#
#
# poolsID %>%
#   filter(is.na(lat.c) | is.na(lon.c))
#
# abundID %>%
#   filter(is.na(lat.c) | is.na(lon.c))

#
# disGeo <- uniqueGeoTable %>%
#   distinct(aName, lat.c, lon.c) %>%
#   #select(-aName) %>%
#   st_as_sf(coords = c("lon.c", "lat.c")) %>%
#   ungroup() %>%
#   distinct(geometry) %>%
#   mutate(geoID = row_number())
#
#
# disGeoJoin <- uniqueGeoTable %>%
#   st_as_sf(coords = c("lon.c", "lat.c")) %>%
#   st_join(disGeo)
#
#
# dupes <- disGeoJoin %>%
#   group_by(geoID) %>%
#   mutate(n = n_distinct(id)) %>%
#   filter(n > 1)
#
# #
# # dupes %>%
# #   inner_join(updatePoints, by = c("aName.x" = "aName")) %>% view()
# #
#
#
# #
# # # Check for NA's
# # uniqueGeoTable  %>%
# #   filter(is.na(lat.c) | is.na(lon.c))
# #
# #
# # # Find locations with same geometry
# # uniqueGeoTable %>%
# #   filter((id == 3891) | (id == 3954))
# #
# #
# # library(tmap)
# # tmap_mode("view")
# #
# # dropMismatch %>%
# #   filter((aName == "suya000108 plumas lake") | (aName == "suyaplumas lake")) %>%
# #   distinct(aName, latitude, longitude) %>%
# #   st_as_sf(coords = c("longitude", "latitude")) %>%
# #   tm_shape() +
# #   tm_dots()
#
#
#
#
#
# geoTableAName <- uniqueGeoTable %>%
#   distinct(aName) %>%
#   pull(aName)
#
# abundAName <- abund %>% distinct(aName) %>% pull(aName)
#
# poolsAName <- pools %>% distinct(aName) %>% pull(aName)
# testAName <-  testResults %>% distinct(aName) %>% pull(aName)
#
# #Function to test different tables names
# nameSetDiffFunc <- function(tibble) {
#   newName <- tibble %>% distinct(aName) %>% pull(aName)
#   setdiff(abundAName, newName)
# }
#
# nameSetDiffFunc(dropMismatch) %>%
#   setdiff(nameSetDiffFunc(aNameRows))
#
#
# nameSetDiffFunc(aNameRows) %>% setdiff(nameSetDiffFunc(dropMismatch))
#
#
# nameSetDiffFunc(uniqueGeoTable) %>%
#
#   setdiff(
#     mismatches %>% pull(aName) %>% unique()
#   )
#
#
# nameSetDiffFunc(updatePoints)
# nameSetDiffFunc(codeName)

SJV_abundance = read_csv("Data/SJV/wnv/0_input/MosquitoAbundance2010.2020.csv")
SJV_pools = read_csv("Data/SJV/wnv/0_input/MosquitoPools2010.2020.csv")
SJV_test = read_csv("Data/SJV/wnv/0_input/MosquitoPoolsTestResults2010.2020.csv")



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

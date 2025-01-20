#' Register-based CoMorbidities in a SWEdish setting (rcmswe)
#'
#' Extract Co-Morbidities from Swedish Population Health Registers
#' @param search_df A data-frame with two columns, column one should be study ID's and column two should be an index date.
#' @param sqlite_path A character vector specifying the path to an SQLite3 database containing registers.
#' @param sqlite_NPR_name A character vector specifying the name of the NPR inside the SQLite database (currently not functioning)
#' @param sqlite_LMED_name A character vector specifying the name of the Prescribed Drug Register inside the SQLite database (currently not functioning)
#' @param NPR A logical parameter dictating if the returning dataset should contain individual groups of co-morbidities (default = TRUE)
#' @param LMED A logical parameter dictating if the returning dataset should be expanded with data from LMED (default = FALSE, currently non-functioning)
#' @param CCI A logical parameter dictating if the returning dataset should contain columns for weighted and unweighted CCI (currently solely based on NPR data)
#' @import tidyverse
#' @import DBI
#' @import magrittr
#' @export
rcmswe <- function(search_df, sqlite_path, sqlite_NPR_name = "PAR", sqlite_LMED_name = "LMED", NPR = TRUE, LMED = TRUE, CCI = TRUE) {

  # Error checking input ----

  # Check that Patient ID is in a numeric format
  try(if(!is.numeric(search_df[[1]]))
    return("Error: Patient ID column is not numeric.")
  )

  # Check that the date can be parsed and parse it to ymd
  if(is.character(search_df[[2]])) {
    options(warn=2)
    search_df[[2]] <- lubridate::ymd(stringr::str_replace_all(search_df[[2]], "-", ""))}

  # Split the data frame into one single row per patient per index date ----

  # Only the first two columns are used and thus only these are retained

  # split_rcmswe() also converts dates into character strings for subsequent analysis

  search_df_split <- rcmswe:::split_rcmswe(search_df[1:2])

  # The function group_split() contained in rcmswe:::split_rcmswe() produces a
  # vctrs::list_of object which is difficult to use. Reformatting this into a
  # normal list is preferred. Following this, apply the extract_comorbs()-function
  # across the split list.

  rcmswe_result <- seq(1:length(search_df_split)) %>%
    map(~search_df_split[[.]]) %>%
    map(~rcmswe:::extract_comorbs(., sqlite_path, sqlite_NPR_name, sqlite_LMED_name, NPR, LMED, CCI))

  # Return the resulting list reduced into a single df ----

  return(seq(1:length(search_df_split)) %>%
           map(~left_join(rcmswe_result[[.]], search_df_split[[.]][,1:2])) %>%
           reduce(full_join))


}

#' Add an ICD-10-based comorbidity
#' @param database_extract A data-frame of the extract from PAR
#' @param resulting_df The current resulting df to append the comorbidity to
#' @param dx_vector A vector of diagnoses to collapse into a searchable reprex string
#' @param comorbidity A string to use as the name of the comorbidity
#' @noRd
add_comorb <- function(database_extract, resulting_df, dx_vector, comorbidity){

  # Sanitize comorbidity name: replace spaces with underscores
  safe_name <- stringr::str_replace_all(comorbidity, " ", "_")

  # Create a usable regex string for ICD codes
  reprex <- paste0("|\<", paste0(dx_vector, collapse = '|\<'))

  # Filter patients with matching diagnoses and earliest date
  ptnts <- database_extract[database_extract$datum >= 19970000,][grep(reprex, database_extract[database_extract$datum >= 19970000,]$diagnos),] %>%
    dplyr::group_by(group) %>%
    dplyr::filter(dplyr::row_number(datum) == 1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(!!paste0('date.', safe_name) := datum,
                  !!paste0('diagnos.', safe_name) := diagnos)

  # Join with resulting_df and create indicator column
  Matrix <- resulting_df %>%
    dplyr::left_join(ptnts, by = c('group' = 'group'), copy = TRUE) %>%
    dplyr::mutate(!!safe_name := dplyr::if_else(!is.na(get(paste0('date.', safe_name))), 1, 0, missing = 0))

  return(Matrix)
}


split_rcmswe <- function(df) {
  df %>%
    rename(LopNr = 1, index_date = 2) %>%
    mutate(index_date = stringr::str_replace_all(as.character(index_date), "-", "")) %>%
    select(LopNr, index_date) %>%
    group_by(LopNr) %>%
    arrange(index_date, .by_group = TRUE) %>%
    mutate(id = row_number()) %>%
    group_by(id) %>%
    group_split() %>%
    return()
}

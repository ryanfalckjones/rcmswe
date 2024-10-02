#' Add an ICD-10-based comorbidity
#' @param database_extract A data-frame of the extract from PAR
#' @param resulting_df The current resulting df to append the comorbidity to
#' @param dx_vector A vector of diagnoses to collapse into a searchable reprex string
#' @param comorbidity A string to use as the name of the comorbidity
#' @noRd
add_comorb <- function(database_extract, resulting_df, dx_vector, comorbidity){

  # Create a usable reprex string
  reprex <- paste0("\\<",
                   paste0(dx_vector, collapse = '|\\<')
  )

  ptnts <- database_extract[database_extract$datum >= 19970000,][grep(reprex,database_extract[database_extract$datum >= 19970000,]$diagnos),] %>%
    group_by(group) %>%
    filter(row_number(datum)==1) %>%
    ungroup() %>%
    rename(!!paste0('date.',comorbidity) := datum,
           !!paste0('diagnos.', comorbidity) := diagnos)

  Matrix <- resulting_df %>%
    left_join(ptnts, by = c('group' = 'group'), copy = T) %>%
    mutate(!!comorbidity := if_else(!is.na(get(paste0('date.', comorbidity))),1,0,missing=0))

  return(Matrix)

}

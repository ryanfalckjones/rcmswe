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

db_search <- function(search_tbl, sqlite_path){
  db <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path, extended_types = TRUE) # Connect to DB
  DBI::dbWriteTable(db, 'temptable', search_tbl, temporary = TRUE) # Write temporary table to SQLite DB
  sql_pt_df <- DBI::dbGetQuery(db,"SELECT P.LopNr, P.UTDATUMA, P.DIAGNOS, temptable.index_date
                                              FROM PAR P
                                              INNER JOIN temptable ON P.LopNr = temptable.LopNr
                                              WHERE P.UTDATUMA < temptable.index_date") %>% 
    as_tibble() %>%
    rename(group = LopNr, datum = UTDATUMA, diagnos = DIAGNOS) %>%
    select(-index_date)
  DBI::dbDisconnect(db)
  return(as_tibble(sql_pt_df))
}

df_search <- function(search_tbl, par_df){
  par_df %>%
    select(LopNr, UTDATUMA, DIAGNOS) %>%
    mutate(UTDATUMA = lubridate::ymd(UTDATUMA)) %>%
    semi_join(search_tbl[["LopNr"]]) %>%
    left_join(search_tbl) %>%
    filter(UTDATUMA < index.date) %>%
    select(-index.date) %>%
    rename(group = LopNr, datum = UTDATUMA, diagnos = DIAGNOS) %>%
    return()
}

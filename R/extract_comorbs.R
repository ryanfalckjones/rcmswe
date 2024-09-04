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
#' @export

# Script for querying the Swedish National Patient Register for pre-existing comorbidities ----

extract_comorbs <- function(search_df, sqlite_path, sqlite_NPR_name = "PAR", sqlite_LMED_name = "LMED", NPR = TRUE, LMED = FALSE, CCI = TRUE){

  # Rename columns for ID and date if misspelled and convert dates to character strings
  temptable <- rename(search_df, LopNr = 1, index_date = 2) %>%
    mutate(search_df, index_date = stringr::str_replace_all(as.character(index_date), "-", ""))

  db <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path, extended_types = TRUE) # Connect to DB

  DBI::dbWriteTable(db, 'temptable', temptable, temporary = TRUE) # Write temporary table to SQLite DB


  patients <- as_tibble(DBI::dbGetQuery(db,"SELECT P.LopNr, P.UTDATUMA, P.DIAGNOS, temptable.index_date
                                              FROM PAR P
                                              INNER JOIN temptable ON P.LopNr = temptable.LopNr
                                              WHERE P.UTDATUMA < temptable.index_date") %>%
                          dplyr::rename(group = LopNr, datum = UTDATUMA, diagnos = DIAGNOS) %>%
                          dplyr::select(group, datum, diagnos)
  )

  Matrix <- distinct(patients, group) # Create object to store the Charlson score with one line per patient/ID.

  #######################################################################################################
  ## CCI

  # Myocardial_infarction
  icd7  <- "\\<420,1"
  icd8  <- "\\<410|\\<411|\\<412,01|\\<412,91"
  icd9  <- "\\<410|\\<412"
  icd10 <- "\\<I21|\\<I22|\\<I252"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9  <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10 <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Myocardial_infarction=datum,diagnos.Myocardial_infarction=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Myocardial_infarction=if_else(!is.na(date.Myocardial_infarction),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Congestive_heart_failure
  icd7  <- "\\<422,21|\\<422,22|\\<434,1|\\<434,2"
  icd8  <- "\\<425,08|\\<425,09|\\<427,0|\\<427,1|\\<428"
  icd9  <- paste(c("\\<402A", "402B", "402X", "404A","404B","404X","425E","425F","425H","425W","425X","428"),collapse="|\\<")
  icd10 <- "\\<I110|\\<I130|\\<I132|\\<I255|\\<I420|\\<I426|\\<I427|\\<I428|\\<I429|\\<I43|\\<I50"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10 <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Congestive_heart_failure=datum,diagnos.Congestive_heart_failure=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Congestive_heart_failure=if_else(!is.na(date.Congestive_heart_failure),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Peripheral_vascular_disease
  icd7  <- "\\<450,1|\\<451|\\<453"
  icd8  <- "\\<440|\\<441|\\<443,1|\\<443,9"
  icd9  <- "\\<440|\\<441|\\<443B|\\<443X|\\<447B|\\<557"
  icd10 <- "\\<I70|\\<I71|\\<I731|\\<I738|\\<I739|\\<I771|\\<I790|\\<I792|\\<K55"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Peripheral_vascular_disease=datum,diagnos.Peripheral_vascular_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Peripheral_vascular_disease=if_else(!is.na(date.Peripheral_vascular_disease),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Cerebrovascular_disease
  icd7  <- paste(c("\\<330",331:334),collapse="|\\<")
  icd8  <- "\\<430|\\<431|\\<432|\\<433|\\<434|\\<435|\\<436|\\<437|\\<438"
  icd9  <- "\\<430|\\<431|\\<432|\\<433|\\<434|\\<435|\\<436|\\<437|\\<438"
  icd10 <- "\\<G45|\\<I60|\\<I61|\\<I62|\\<I63|\\<I64|\\<I67|\\<I69"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Cerebrovascular_disease=datum,diagnos.Cerebrovascular_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Cerebrovascular_disease=if_else(!is.na(date.Cerebrovascular_disease),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Chronic_obstructive_pulmonary_disease
  icd7  <- "\\<502|\\<527,1"
  icd8  <- "\\<491|\\<492"
  icd9  <- "\\<491|\\<492|\\<496"
  icd10 <- "\\<J43|\\<J44"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Chronic_obstructive_pulmonary_disease=datum,diagnos.Chronic_obstructive_pulmonary_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Chronic_obstructive_pulmonary_disease=if_else(!is.na(date.Chronic_obstructive_pulmonary_disease),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Chronic_other_pulmonary_disease
  icd7  <- paste(c("\\<241",501,523:526),collapse="|\\<")
  icd8  <- paste(c("\\<490",493,515:518),collapse="|\\<")
  icd9  <- paste(c("\\<490",493:495,500:508,516,517),collapse="|\\<")
  icd10 <- paste(c("\\<J45",41,42,46,47,60:70),collapse="|\\<J")

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Chronic_other_pulmonary_disease=datum,diagnos.Chronic_other_pulmonary_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Chronic_other_pulmonary_disease=if_else(!is.na(date.Chronic_other_pulmonary_disease),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Rheumatic_disease
  icd7  <- paste(c("\\<722,00","722,01","722,10","722,20","722,23","456,0","456,1","456,2","456,3"),collapse="|\\<")
  icd8  <- paste(c("\\<446",696,"712,0","712,1","712,2","712,3","712,5", 716, "734,0", "734,1", "734,9"),collapse="|\\<")
  icd9  <- paste(c("\\<446","696A","710A","710B","710C","710D","710E",714,"719D",720,725),collapse="|\\<")
  icd10 <- paste(c("\\<M05","06",123,"070","071","072","073","08",13,30,313:316,32:34,350:351,353,45:46),collapse="|\\<M")

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Rheumatic_disease=datum,diagnos.Rheumatic_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Rheumatic_disease=if_else(!is.na(date.Rheumatic_disease),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Dementia
  icd7  <- "\\<304|\\<305"
  icd8  <- "\\<290"
  icd9  <- "\\<290|\\<294B|\\<331A|\\<331B|\\<331C|\\<331X"
  icd10 <- "\\<F00|\\<F01|\\<F02|\\<F03|\\<F051|\\<G30|\\<G311|\\<G319"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Dementia=datum,diagnos.Dementia=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Dementia=if_else(!is.na(date.Dementia),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Hemiplegia
  icd7 <- "\\<351|\\<352|\\<357,00"
  icd8 <- "\\<343|\\<344"
  icd9 <- "\\<342|\\<343|\\<344A|\\<344B|\\<344C|\\<344D|\\<344E|\\<344F"
  icd10 <- "\\<G114|\\<G80|\\<G81|\\<G82|\\<G830|\\<G831|\\<G832|\\<G833|\\<G838"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Hemiplegia=datum,diagnos.Hemiplegia=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Hemiplegia=if_else(!is.na(date.Hemiplegia),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Diabetes_without_chronic_complication
  icd7 <- "\\<260,09"
  icd8 <-  "\\<250,00|\\<250,07|\\<250,08"
  icd9 <- "\\<250A|\\<250B|\\<250C"
  icd10 <- "\\<E100|\\<E101|\\<E110|\\<E111|\\<E120|\\<E121|\\<E130|\\<E131|\\<E140|\\<E141"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Diabetes_without_chronic_complication=datum,diagnos.Diabetes_without_chronic_complication=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Diabetes_without_chronic_complication=if_else(!is.na(date.Diabetes_without_chronic_complication),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Diabetes_with_chronic_complication
  icd7 <- "\\<260,2|\\<260,21|\\<260,29|\\<260,3|\\<260,4|\\<260,49|\\<260,99"
  icd8 <- "\\<250,01|\\<250,02|\\<250,03|\\<250,04|\\<250,05"
  icd9 <- "\\<250D|\\<250E|\\<250F|\\<250G"
  icd10 <- "\\<E102|\\<E103|\\<E104|\\<E105|\\<E107|\\<E112|\\<E113|\\<E114|\\<E115|\\<E116|\\<E117|\\<E122|\\<E123|\\<E124|\\<E125|\\<E126|\\<E127|\\<E132|\\<E133|\\<E134|\\<E135|\\<E136|\\<E137|\\<E142|\\<E143|\\<E144|\\<E145|\\<E146|\\<E147"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Diabetes_with_chronic_complication=datum,diagnos.Diabetes_with_chronic_complication=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Diabetes_with_chronic_complication=if_else(!is.na(date.Diabetes_with_chronic_complication),1,0,missing=0))
  Matrix <- Matrix %>% mutate(Diabetes_without_chronic_complication=if_else(Diabetes_with_chronic_complication==1,0,Diabetes_without_chronic_complication))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Renal_disease
  icd7 <- "\\<592|\\<593|\\<792"
  icd8 <- "\\<582|\\<583|\\<584|\\<792|\\<593|\\<403,99|\\<404,99|\\<792,99|\\<Y29,01"
  icd9 <- "\\<403A|\\<403B|\\<403X|\\<582|\\<583|\\<585|\\<586|\\<588A|\\<V42A|\\<V45B|\\<V56"
  icd10 <- "\\<I120|\\<I131|\\<N032|\\<N033|\\<N034|\\<N035|\\<N036|\\<N037|\\<N052|\\<N053|\\<N054|\\<N055|\\<N056|\\<N057|\\<N11|\\<N18|\\<N19|\\<N250|\\<Q611|\\<Q612|\\<Q613|\\<Q614|\\<Z49|\\<Z940|\\<Z992"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Renal_disease=datum,diagnos.Renal_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Renal_disease=if_else(!is.na(date.Renal_disease),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Mild_liver_disease
  icd7  <- "\\<581"
  icd8  <- "\\<070|\\<571|\\<573"
  icd9 <-  "\\<070|\\<571C|\\<571E|\\<571F|\\<573"
  icd10 <- "\\<B15|\\<B16|\\<B17|\\<B18|\\<B19|\\<K703|\\<K709|\\<K73|\\<K746|\\<K754"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Mild_liver_disease=datum,diagnos.Mild_liver_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Mild_liver_disease=if_else(!is.na(date.Mild_liver_disease),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # liver special
  icd8  <- "\\<785,3"
  icd9 <- "\\<789F"
  icd10 <- "\\<R18"

  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.liver_special=datum,diagnos.liver_special=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Liver_special=if_else(!is.na(date.liver_special),1,0,missing=0))
  rm(icd8, icd9, icd10, ICD8, ICD9, ICD10, ptnts)

  # moderate severe liver disease
  icd7 <- "\\<462,1"
  icd8 <- "\\<456,0|\\<571,9|\\<573,02"
  icd9 <- "\\<456A|\\<456B|\\<456C|\\<572C|\\<572D|\\<572E"
  icd10 <-  "\\<I850|\\<I859|\\<I982|\\<I983"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Severe_liver_disease=datum,diagnos.Severe_liver_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Severe_liver_disease=if_else(!is.na(date.Severe_liver_disease),1,0,missing=0))
  Matrix <- Matrix %>% mutate(Severe_liver_disease=if_else(Mild_liver_disease==1 & Liver_special==1,1,Severe_liver_disease))
  Matrix <- Matrix %>% mutate(Mild_liver_disease=if_else(Severe_liver_disease==1,0,Mild_liver_disease))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Peptic_ulcer_disease
  icd7  <- "\\<540|\\<541|\\<542"
  icd8  <- "\\<531|\\<532|\\<533|\\<534"
  icd9 <- "\\<531|\\<532|\\<533|\\<534"
  icd10 <-"\\<K25|\\<K26|\\<K27|\\<K28"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Peptic_ulcer_disease=datum,diagnos.Peptic_ulcer_disease=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Peptic_ulcer_disease=if_else(!is.na(date.Peptic_ulcer_disease),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Malignancy
  icd7   <- paste(paste("\\<",paste(140:190,collapse = "|\\<"),sep=""), paste("|\\<",paste(192:197,collapse = "|\\<"),sep=""), paste("|\\<",paste(200:204,collapse = "|\\<"),sep=""),sep="")
  icd8   <- paste(paste("\\<",paste(c(140:172,174),collapse = "|\\<"),sep=""), paste("|\\<",paste(c(180:207,209),collapse = "|\\<"),sep=""),sep="")
  icd9   <- paste(paste("\\<",paste(140:172,collapse = "|\\<"),sep=""), paste("|\\<",paste(174:208,collapse = "|\\<"),sep=""),sep="")
  icd10  <- paste("\\<C00|\\<C0",paste(1:9,collapse = "|\\<C0",sep=""),paste("|\\<C",paste(c(10:41,43,45:58,60:76,81:86,88:97),collapse = "|\\<C"),sep=""),sep="")

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.malignancy=datum,diagnos.malignancy=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Malignancy=if_else(!is.na(date.malignancy),1,0,missing=0))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Metastatic_cancer
  icd7 <- "\\<156,91|\\<198|\\<199"
  icd8 <- "\\<196|\\<197|\\<198|\\<199"
  icd9 <- "\\<196|\\<197|\\<198|\\<199A|\\<199B"
  icd10 <- "\\<C77|\\<C78|\\<C79|\\<C80"

  ICD7  <- patients[patients$datum<19690000,][grep(icd7,patients[patients$datum<19690000,]$diagnos),]
  ICD8  <- patients[patients$datum >= 19690000 & patients$datum < 19870000,][grep(icd8,patients[patients$datum >= 19690000 & patients$datum < 19870000,]$diagnos),]
  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD7,ICD8,ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Metastatic_solid_tumor=datum,diagnos.Metastatic_solid_tumor=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Metastatic_solid_tumor=if_else(!is.na(date.Metastatic_solid_tumor),1,0,missing=0))
  Matrix <- Matrix %>% mutate(Malignancy=if_else(Metastatic_solid_tumor==1,0,Malignancy))
  rm(icd7, icd8, icd9, icd10, ICD7, ICD8, ICD9, ICD10, ptnts)

  # Aids
  icd9  <- "\\<079J|\\<279K"
  icd10 <- "\\<B20|\\<B21|\\<B22|\\<B23|\\<B24|\\<F024|\\<O987|\\<R75|\\<Z219|\\<Z717"

  ICD9 <- patients[patients$datum >= 19870000 & patients$datum < 19980000,][grep(icd9,patients[patients$datum >= 19870000 & patients$datum < 19980000,]$diagnos),]
  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- bind_rows(ICD9,ICD10) %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Aids=datum,diagnos.Aids=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Aids=if_else(!is.na(date.Aids),1,0,missing=0))
  rm(icd9, icd10, ICD9, ICD10, ptnts)

  ###
  ### Non-CCI comorbidities
  ###

  # Hypertension
  icd10 <- "\\<I10|\\<I11|\\<I12|\\<I13|\\<I14|\\<I15"

  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- ICD10 %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Hypertension=datum,diagnos.Hypertension=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(Hypertension=if_else(!is.na(date.Hypertension),1,0,missing=0))
  rm(icd10, ICD10, ptnts)

  # Atrial Fibrillation
  icd10 <- "\\<I48"

  ICD10  <- patients[patients$datum >= 19970000,][grep(icd10,patients[patients$datum >= 19970000,]$diagnos),]
  ptnts <- ICD10 %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.AFib=datum,diagnos.AFib=diagnos)
  Matrix <- left_join(Matrix,ptnts,by=c("group"="group"),copy=T)
  Matrix <- Matrix %>% mutate(AFib=if_else(!is.na(date.AFib),1,0,missing=0))
  rm(icd10, ICD10, ptnts)

  # Calculate CCI if requested

  if(CCI){
    # Calculate the unweighted comorbidity index
    Matrix$CCIunw <- Matrix$Myocardial_infarction + Matrix$Congestive_heart_failure + Matrix$Peripheral_vascular_disease +
      Matrix$Cerebrovascular_disease + Matrix$Chronic_obstructive_pulmonary_disease + Matrix$Chronic_other_pulmonary_disease +
      Matrix$Rheumatic_disease + Matrix$Dementia + Matrix$Hemiplegia + Matrix$Diabetes_without_chronic_complication +
      Matrix$Diabetes_with_chronic_complication + Matrix$Renal_disease + Matrix$Mild_liver_disease + Matrix$Severe_liver_disease +
      Matrix$Peptic_ulcer_disease + Matrix$Malignancy + Matrix$Metastatic_solid_tumor + Matrix$Aids

    # Calculate the weighted comorbidity index
    Matrix$CCIw <- Matrix$Myocardial_infarction + Matrix$Congestive_heart_failure + Matrix$Peripheral_vascular_disease +
      Matrix$Cerebrovascular_disease + Matrix$Chronic_obstructive_pulmonary_disease + Matrix$Chronic_other_pulmonary_disease +
      Matrix$Rheumatic_disease + Matrix$Dementia + 2*Matrix$Hemiplegia + Matrix$Diabetes_without_chronic_complication +
      2*Matrix$Diabetes_with_chronic_complication + 2*Matrix$Renal_disease + Matrix$Mild_liver_disease + 3*Matrix$Severe_liver_disease +
      Matrix$Peptic_ulcer_disease + 2*Matrix$Malignancy + 6*Matrix$Metastatic_solid_tumor + 6*Matrix$Aids
  }

  # Delete date and diagnos information in case not needed
  Matrix <- select(Matrix, -contains(".")) %>%
    rename(LopNr = group) # Rename the ID column to LopNr


  # Extend with LMED if requested
  if(LMED){

    # LMED stores dates in a Date format rather than yyyymmdd
    temptable2 <- temptable %>%
      mutate(index_date = ymd(index_date))

    DBI::dbWriteTable(db, 'temptable2', temptable2, temporary = TRUE) # Write temporary table to SQLite DB

    drugs <- as_tibble(DBI::dbGetQuery(db,"SELECT L.LopNr, L.ATC, L.FDATUM, temptable2.index_date
                                           FROM LMED L
                                           INNER JOIN temptable2 ON L.LopNr = temptable2.LopNr
                                           WHERE L.FDATUM < temptable2.index_date") %>%
                         dplyr::mutate(group = LopNr, datum = stringr::str_replace_all(as.character(FDATUM), "-", ""), drug = ATC) %>%
                         dplyr::select(group, datum, drug)
    )

    Matrix_drugs <- distinct(drugs, group) # Create object to store the Charlson score with one line per patient/ID.

    # Diabetes
    atc <- "\\<A10"
    ATC  <- drugs[drugs$datum >= 19970000,][grep(atc,drugs[drugs$datum >= 19970000,]$drug),]
    drgs <- ATC %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Diabetes=datum,drug.Diabetes=drug)
    Matrix_drugs <- left_join(Matrix_drugs,drgs,by=c("group"="group"),copy=T)
    Matrix_drugs <- Matrix_drugs %>% mutate(Diabetes=if_else(!is.na(date.Diabetes),1,0,missing=0))
    rm(atc, ATC, drgs)

    # Hypertension
    atc <- "\\<C08C|\\<C03A"
    ATC  <- drugs[drugs$datum >= 19970000,][grep(atc,drugs[drugs$datum >= 19970000,]$drug),]
    drgs <- ATC %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Hypertension=datum,drug.Hypertension=drug)
    Matrix_drugs <- left_join(Matrix_drugs,drgs,by=c("group"="group"),copy=T)
    Matrix_drugs <- Matrix_drugs %>% mutate(Hypertension=if_else(!is.na(date.Hypertension),1,0,missing=0))
    rm(atc, ATC, drgs)

    # Hyperlipidaemia
    atc <- "\\<C10"
    ATC  <- drugs[drugs$datum >= 19970000,][grep(atc,drugs[drugs$datum >= 19970000,]$drug),]
    drgs <- ATC %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Hyperlipidaemia=datum,drug.Hyperlipidaemia=drug)
    Matrix_drugs <- left_join(Matrix_drugs,drgs,by=c("group"="group"),copy=T)
    Matrix_drugs <- Matrix_drugs %>% mutate(Hyperlipidaemia=if_else(!is.na(date.Hyperlipidaemia),1,0,missing=0))
    rm(atc, ATC, drgs)

    # CKD
    atc <- "\\<B05D|\\<B05Z"
    ATC  <- drugs[drugs$datum >= 19970000,][grep(atc,drugs[drugs$datum >= 19970000,]$drug),]
    drgs <- ATC %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.CKD=datum,drug.CKD=drug)
    Matrix_drugs <- left_join(Matrix_drugs,drgs,by=c("group"="group"),copy=T)
    Matrix_drugs <- Matrix_drugs %>% mutate(CKD=if_else(!is.na(date.CKD),1,0,missing=0))
    rm(atc, ATC, drgs)

    # Dementia
    atc <- "\\<N06D"
    ATC  <- drugs[drugs$datum >= 19970000,][grep(atc,drugs[drugs$datum >= 19970000,]$drug),]
    drgs <- ATC %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Dementia=datum,drug.Dementia=drug)
    Matrix_drugs <- left_join(Matrix_drugs,drgs,by=c("group"="group"),copy=T)
    Matrix_drugs <- Matrix_drugs %>% mutate(Dementia=if_else(!is.na(date.Dementia),1,0,missing=0))
    rm(atc, ATC, drgs)

    # Ischaemic Heart Disease
    atc <- "\\<C01D"
    ATC  <- drugs[drugs$datum >= 19970000,][grep(atc,drugs[drugs$datum >= 19970000,]$drug),]
    drgs <- ATC %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.IHD=datum,drug.IHD=drug)
    Matrix_drugs <- left_join(Matrix_drugs,drgs,by=c("group"="group"),copy=T)
    Matrix_drugs <- Matrix_drugs %>% mutate(IHD=if_else(!is.na(date.IHD),1,0,missing=0))
    rm(atc, ATC, drgs)

    # Congestive Heart Failure
    atc <- "\\<C01A"
    ATC  <- drugs[drugs$datum >= 19970000,][grep(atc,drugs[drugs$datum >= 19970000,]$drug),]
    drgs <- ATC %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.Congestive_heart_failure=datum,drug.Congestive_heart_failure=drug)
    Matrix_drugs <- left_join(Matrix_drugs,drgs,by=c("group"="group"),copy=T)
    Matrix_drugs <- Matrix_drugs %>% mutate(Congestive_heart_failure=if_else(!is.na(date.Congestive_heart_failure),1,0,missing=0))
    rm(atc, ATC, drgs)

    # Treated ADHD
    atc <- "\\<N06BA"
    ATC  <- drugs[drugs$datum >= 19970000,][grep(atc,drugs[drugs$datum >= 19970000,]$drug),]
    drgs <- ATC %>% group_by(group) %>% filter(row_number(datum)==1) %>% ungroup %>% rename(date.ADHD=datum,drug.ADHD=drug)
    Matrix_drugs <- left_join(Matrix_drugs,drgs,by=c("group"="group"),copy=T)
    Matrix_drugs <- Matrix_drugs %>% mutate(ADHD=if_else(!is.na(date.ADHD),1,0,missing=0))
    rm(atc, ATC, drgs)

    # Delete date and other information in case not needed
    Matrix_drugs <- select(Matrix_drugs, -contains(".")) %>%
      rename(LopNr = group)

    # Expand the NPR-based comorbidities with LMED-data
    LMED_Matrix <- Matrix %>%
      full_join(Matrix_drugs %>% select(LopNr, Hypertension, Congestive_heart_failure)) %>%
      group_by(LopNr) %>%
      summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))

    # Retun the data set (temporary while testing)
    return(LMED_Matrix)
  }

  # Return the requested dataset
  if(NPR & !LMED){return(Matrix)} # NPR only (with/without CCI depending on if it is requested)
  if(!NPR & CCI){return(Matrix %>% select(LopNr, CCIunw, CCIw))} # Only CCI (based on NPR only)
}

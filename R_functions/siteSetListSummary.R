library(RSQLite)
library(TFBSTools)

siteSetListSummary = function(query, siteSetList, updateProgress) {
  x = siteSetList
  # x = scannedRange
  
  if (class(query) == 'GRanges') {
    print('GR') 
  } else {
    stop('`query` must be a GenomicRanges object')
  }

  # queryLength = length(query)
  ########### STEP 1 ###############
  if (is.function(updateProgress)) {
    text <- ('gathering results...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print('gathering results...')
  }
  
  motif_ID = unlist(lapply(x, function(z) {ID(z@pattern)}))
  motif_lengths = unlist(lapply(x, function(z) {length(z@views)}))
  motifSummary = as.data.frame(t(rbind(motif_ID, motif_lengths)), row.names = F) %>%
    dplyr::mutate(motif_ID.version = motif_ID) %>%
    dplyr::mutate(motif_version = stringr::str_sub(motif_ID, start = 8L, end = 9L)) %>%
    dplyr::mutate(motif_ID = stringr::str_sub(motif_ID, start = 1L, end = 6L)) %>%
    dplyr::mutate(motif_lengths = as.numeric(motif_lengths))
  
  
  # get results straight from JASPAR SQLite db
  drv = RSQLite::SQLite()
  con = dbConnect(drv, dbname = JASPAR2016::JASPAR2016@db)
  res = gsub(capture.output(cat(motifSummary$motif_ID, sep = ',')), pattern = ',', replacement = '","')
  resu = dbGetQuery(con, paste0('SELECT * FROM MATRIX WHERE BASE_ID IN ("', res, '")')) %>%
    dplyr::select(-ID, -COLLECTION) %>%
    tidyr::unite(col = MOTIF_ID, BASE_ID, VERSION, sep = '.') %>%
    # dplyr::rename(MOTIF_ID = BASE_ID) %>%
    dplyr::filter(MOTIF_ID %in% motifSummary$motif_ID.version)

  # queryLength = length(query)
  ########### STEP 1 ###############
  if (is.function(updateProgress)) {
    text <- ('summarizing...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print('summarizing...')
  }
  
  motifSummary = motifSummary %>%
    group_by(motif_ID.version) %>%
    summarise(n = sum(motif_lengths))
  
  result = resu %>%
    inner_join(., motifSummary, by = c('MOTIF_ID' = 'motif_ID.version')) %>%
    arrange(-n)
  result
}

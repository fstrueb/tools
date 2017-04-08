library(RSQLite)

siteSetListSummary = function(query, siteSetList, ) {
  x = siteSetList
  # x = scannedRange
  
  if (class(query) == 'GRanges') {
    print('GR') 
  } else {
    stop('`query` must be a GenomicRanges object')
  }
  
  # first set the GRanges query frame
  # queryLength = length(query)
  
  motif_ID = unlist(lapply(x, function(z) {ID(z@pattern)}))
  motif_lengths = unlist(lapply(x, function(z) {length(z@views)}))
  motifSummary = as.data.frame(t(rbind(motif_ID, motif_lengths)), row.names = F) %>%
    # dplyr::mutate(motif_ID = stringr::str_sub(motif_IDs, start = 1L, end = 6L)) %>%
    dplyr::mutate(motif_lengths = as.numeric(motif_lengths)) %>%
    group_by(motif_ID) %>%
    summarise(n = sum(motif_lengths))
  
  # get results straight from JASPAR SQLite db
  drv = RSQLite::SQLite()
  con = dbConnect(drv, dbname = JASPAR2016::JASPAR2016@db)
  res = gsub(capture.output(cat(motifSummary$motif_ID, sep = ',')), pattern = ',', replacement = '","')
  resu = dbGetQuery(con, paste0('SELECT * FROM MATRIX WHERE BASE_ID IN ("', res, '")')) %>%
    dplyr::select(-ID, -COLLECTION) %>%
    tidyr::unite(col = MOTIF_ID, BASE_ID, VERSION, sep = '.')
  
  result = resu %>%
    inner_join(., motifSummary, by = c('MOTIF_ID' = 'motif_ID'))
  result
}

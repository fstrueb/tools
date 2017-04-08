# function to list all JASPAR 2016 TFBS matrices
library(RSQLite)

unlistJASPAR = function(species, collection) {
  drv = RSQLite::SQLite()
  con = dbConnect(drv, dbname = JASPAR2016::JASPAR2016@db)
  
  # dbListTables(con)
  # dbListFields(con, "MATRIX")
  # dbGetQuery(con, "SELECT * FROM MATRIX")
  
  ######## info if needed
  # collections = dbGetQuery(con, "SELECT COLLECTION FROM MATRIX") %>% unique()
  # names = dbGetQuery(con, "SELECT NAME FROM MATRIX") %>% unique()
  # matrixIDs = dbGetQuery(con, "SELECT BASE_ID FROM MATRIX")
  # matrixVersion = dbGetQuery(con, "SELECT VERSION FROM MATRIX")
  # matrices = cbind(matrixIDs, matrixVersion) %>% tidyr::unite(motif_ID, BASE_ID, VERSION, sep = '.')
  # matrices = matrices$motif_ID
  # species = dbGetQuery(con, "SELECT * FROM MATRIX_SPECIES") %>% group_by(TAX_ID) %>% tally() %>% arrange(-n)
  
  ######## create result objects
  # if (species == 'mouse') {
  #   x = 10090
  # } else if (species == 'human') {
  #   x = 9606
  # } else if (is.numeric(species)) {
  #   x = species
  # } else {
  #   stop('x is neither mouse, human, or a species ID')
  # }
  
  # if (is.numeric(species)) {
  #   x = species
  # } else {
  #   warning('`species` must be a numeric vector with taxon IDs')
  # }
  
  ####### collection should be a character vector
  if (is.character(collection)) {
    y = collection
  } else {
    warning('`collection` should be a character vector')
  }
  x = as.numeric(unlist(species))
  cat(x)
  # motifs = lapply(y, function(z) {
  #   dbGetQuery(con, paste('SELECT MATRIX.ID, MATRIX.BASE_ID, MATRIX.VERSION, MATRIX.NAME, MATRIX_SPECIES.TAX_ID FROM MATRIX, MATRIX_SPECIES WHERE MATRIX.ID = MATRIX_SPECIES.ID AND MATRIX_SPECIES.TAX_ID = :x AND MATRIX.COLLECTION = "', z, '"', sep=""), params = list(x = as.character(x)))
  # })
  # result = do.call(rbind, motifs)
  
  # must transform input to strange string objects for use in SQLite query
  # y = c('CORE', 'PBM')
  # x = list('9606', '10090')
  xx = gsub(capture.output(cat(x, sep = ',')), pattern = ',', replacement = '","')
  yy = gsub(capture.output(cat(y, sep = ',')), pattern = ',', replacement = '","')
  result = dbGetQuery(con, paste0('SELECT MATRIX.ID, MATRIX.BASE_ID, MATRIX.VERSION, MATRIX.NAME, MATRIX_SPECIES.TAX_ID FROM MATRIX, MATRIX_SPECIES WHERE MATRIX.ID = MATRIX_SPECIES.ID AND MATRIX_SPECIES.TAX_ID IN ("', xx, '") AND MATRIX.COLLECTION IN ("', yy, '")'))
  
  result %>% tidyr::unite(motif_ID, BASE_ID, VERSION, sep = '.', remove = F)
}

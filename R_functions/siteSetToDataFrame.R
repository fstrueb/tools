siteSetToDataFrame = function(input, query, return.sequence = F, return.p.val = F, updateProgress, ...) {
scannedRange = input
############ STEP 1 ##############
# insert status message how many tfbs were found
if (is.function(updateProgress)) {
  text <- paste('Gathering results for', length(scannedRange), 'TFBS...')
  updateProgress(detail = text)
} else {
  # extract DNA seq from promoters, return a DNAStringSet object
  print(paste('Gathering results for', length(scannedRange), 'TFBS...'))
}

####### try different types of tranlating searchSeq output to human-readable format
if (return.sequence) {
  result = TFBSTools::writeGFF3(scannedRange, scoreType = 'absolute')
} else {
  result = list()
  for (x in 1:length(scannedRange)) {
    xx = scannedRange[[x]]
    # extract elements from SiteSetList
    result[[x]] = list(
      seqname = as.character(xx@seqname),
      start = as.numeric(xx@views@ranges@start),
      end = as.numeric(xx@views@ranges@start + xx@views@ranges@width),
      score = as.numeric(xx@score),
      strand = as.character(xx@strand),
      motif_name = as.character(TFBSTools::tags(xx@pattern)$symbol),
      motif_ID = as.character(TFBSTools::ID(xx@pattern)),
      motif_species = as.character(TFBSTools::tags(xx@pattern)$species),
      motif_evidence = as.character(TFBSTools::tags(xx@pattern)$type)
    )
  }
  
  ############ STEP 2 ##############
  # insert status message how many tfbs were found
  if (is.function(updateProgress)) {
    text <- paste('TFBS extraction completed, converting objects...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print(paste('TFBS extraction completed, converting objects...'))
  }
  
  # filter out NULL results
  result2 = result[!sapply(result, is.null)] 
  result3 = as.data.frame(do.call(rbind, result2)) %>%
    dplyr::filter(start != 'numeric(0)')
  # MAKE A CHARACTER DATA FRAME 
  y = as.data.frame(lapply(result3, as.character), stringsAsFactors = F) %>%
    dplyr::mutate(mergeID = row_number())
  fun = function(z) {
    maxCols = max(sapply(strsplit(as.character(z$start),','), length))
    unnestedCols = tidyr::separate(z, start, paste0('start_', 1:maxCols), sep = ',') %>%
      tidyr::separate(end, paste0('end_', 1:maxCols), sep = ',') %>%
      tidyr::separate(score, paste0('score_', 1:maxCols), sep = ',') %>%
      tidyr::separate(strand, paste0('strand_', 1:maxCols), sep = ',') %>%
      tidyr::gather(key, value, dplyr::starts_with('start'), dplyr::starts_with('end'), dplyr::starts_with('score'), dplyr::starts_with('strand')) %>%
      dplyr::mutate(value = gsub(x = value, pattern = 'c\\(', replacement = '')) %>%
      dplyr::mutate(value = gsub(x = value, pattern = '\\)', replacement = '')) %>%
      dplyr::mutate(value = gsub(x = value, pattern = '\"', replacement = '')) %>%
      tidyr::separate(key, c('key', 'motif_hit'), sep = '_') %>%
      tidyr::spread(key, value) %>%
      dplyr::select(motif_ID, motif_hit, seqname, start, end, strand, dplyr::everything())
    return(unnestedCols)
  }
  ############ STEP 3 ##############
  # insert status message how many tfbs were found
  if (is.function(updateProgress)) {
    text <- paste('Character conversion completed, ordering rows...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print(paste('Character conversion completed, ordering rows...'))
  }
  
  result = purrrlyr::by_row(y, fun, .collate = 'rows')[12:22]
  ############ STEP 4 ##############
  # insert status message how many tfbs were found
  if (is.function(updateProgress)) {
    text <- paste('Concatenation successful, please wait...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print(paste('Concatenation successful, please wait...'))
  }
}

## implement later start ----
## scan using pairwise alignments
# chain = import.chain('mm9.hg19.all.chain')
# scannedRangePhylo = searchPairBSgenome(motif, anno, BSgenome.Hsapiens.NCBI.GRCh38, chr1="chr1", chr2="1", min.score="90%", strand="*", chain=chain)
# test = writeGFF3(scannedRangePhylo)
## implement later end ----

# # return the combined dataframe
# stopifnot(length(dna) == length(query))

if (return.sequence) {
  result = with(result, cbind(
    result,
    reshape::colsplit(
      attributes,
      split = '\\;',
      names = c('TF', 'TF_class', 'TF_sequence')
    )
  )) %>%
    dplyr::mutate(seqname = as.numeric(seqname)) %>%
    dplyr::select(-source, -feature, -frame, -attributes) %>%
    dplyr::rename(
      motif_start = start,
      motif_end = end,
      motif_strand = strand,
      motif_score = score
    )
} else {
  result = result %>%
    dplyr::rename(
      motif_start = start,
      motif_end = end,
      motif_strand = strand,
      motif_score = score
    ) %>%
    dplyr::mutate(seqname = as.numeric(seqname))
  
}

# caveat --> naming convention of column 'geneSymbol'
helper = query[result$seqname]@elementMetadata$geneSymbol
seqnames = query[result$seqname]@seqnames
seqranges = as.data.frame(query[result$seqname]@ranges)
names(seqranges) = c('query_start', 'query_end', 'query_width')
strand = query[result$seqname]@strand
if (is.null(helper)) {
  helper = rep('undefined', times = length(strand))
}
result2 = cbind(helper, seqnames, seqranges, strand, result)

if (return.p.val) {
  ############ STEP 5, optional ##############
  # insert status message how many tfbs were found
  if (is.function(updateProgress)) {
    text <- 'calculating p-values...'
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print('calculating p-values...')
  }
  scannedRange.pval = TFBSTools::pvalues(scannedRange)
  pval = unlist(scannedRange.pval)
  result2 = cbind(result2, pval)
}

# cosmetic changes
result2 = result2 %>%
  dplyr::select(-seqname) %>%
  dplyr::rename(
    promoter = helper,
    query_strand = strand,
    query_chromosome = seqnames
  ) %>%
  dplyr::mutate(motif_start = as.numeric(motif_start), 
    motif_end = as.numeric(motif_end)) %>%
  dplyr::mutate(motif_strand = stringr::str_trim(motif_strand, side = 'both')) %>%
  dplyr::mutate(motif_hit = as.numeric(motif_hit))

return(result2)
}
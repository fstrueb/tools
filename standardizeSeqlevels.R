# standardizeSeqlevels

standardizeSeqlevels = function(range,
                                seqlevels = 'mouse',
                                type = 'UCSC') {
  # check input
  if (class(range) != 'GRanges') {
    stop('You must supply a GenomicRanges object for range')
  }
  
  # define standard seqlevels
  if (seqlevels == 'mouse') {
    seqlevelsNew = c(seq(from = 1, to = 19), 'X', 'Y', 'M')
  } else if (seqlevels == 'human') {
    seqlevelsNew = c(seq(from = 1, to = 21), 'X', 'Y', 'M')
  }
  
  # check which seqlevels were in the input
  if (type == 'UCSC') {
    if (any(grepl(x = as.character(range@seqnames), pattern = 'chr'))) {
      seqlevelsNew = paste0('chr', seqlevelsNew)
      range = GenomeInfoDb::keepSeqlevels(x = range, value = seqlevelsNew)
    }
    else {
      range = GenomeInfoDb::keepSeqlevels(x = range, value = seqlevelsNew)
      rangeDF = as.data.frame(range) %>%
        dplyr::mutate(seqnames = paste0('chr', seqnames))
      range = GenomicRanges::makeGRangesFromDataFrame(rangeDF, keep.extra.columns = T)
    }
  }
  
  if (type == 'ensembl') {
    if (any(grepl(x = as.character(range@seqnames), pattern = 'chr'))) {
      rangeDF = as.data.frame(range) %>%
        dplyr::mutate(seqnames = gsub(seqnames, pattern = 'chr', replacement = ''))
      range = GenomicRanges::makeGRangesFromDataFrame(rangeDF, keep.extra.columns = T)
      range = GenomeInfoDb::keepSeqlevels(x = range, value = seqlevelsNew)
    }
    
    else {
      range = GenomeInfoDb::keepSeqlevels(x = range, value = seqlevelsNew)
    }
  
  }
  return(range)
}
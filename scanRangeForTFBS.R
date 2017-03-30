# function to scan promoters defined by CAGE expression of dataset of interest for TFBS
# the regulatory elements to scan for TFBS should come in a GRanges format, where every row/region of interest should correspond to exactly one regulatory element (promoter/enhancer/genomic anchor)

scanRangeForTFBS = function(query,
                                  motif.list,
                                  input.assembly,
                                  return.p.val = FALSE,
                                  return.sequence = FALSE,
                                  ...) {
  options(stringsAsFactors = FALSE)
  library(magrittr)
  
  anno = vector()
  if (input.assembly == 'mm10') {
    anno = BSgenome.Mmusculus.UCSC.mm10::Mmusculus
  }
  if (input.assembly == 'mm9') {
    anno = BSgenome.Mmusculus.UCSC.mm9::Mmusculus
  }
  
  if (class(query) != 'GRanges') {
    stop('Please provide a valid GenomicRanges object for query')
  }
  
  ########### STEP 1 ###############
  if (is.function(updateProgress)) {
    text <- ('fetching DNA...')
    updateProgress(detail = text)
  } else {
  # extract DNA seq from promoters, return a DNAStringSet object
  print('fetching DNA...')
  }
  dna = Biostrings::getSeq(anno, query)
  
  # motif.list can either be a vector with JASPAR IDs or character string 'all'
  opts = list()
  if (motif.list == 'all') {
    opts['all'] = T
  } else {
    opts[['ID']] = motif.list
  }
  ########### STEP 2 ###############
  if (is.function(updateProgress)) {
    text <- ('fetching TFBS matrices...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print('fetching TFBS matrices...')
  }
  motif = TFBSTools::toPWM(TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, opts))
  
  ########### STEP 3 ###############
  if (is.function(updateProgress)) {
    text <- ('scanning the sequence...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print('scanning the sequence...')
  }
  # scan the sequences!
  scannedPromoters = TFBSTools::searchSeq(
    x = motif,
    subject = dna,
    min.score = '90%',
    mc.cores = 8L
  )
  
  ############ STEP 4 ##############
  # insert status message how many tfbs were found
  if (is.function(updateProgress)) {
    text <- paste('found a total of', length(scannedPromoters), 'TFBS, gathering results...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print(paste('found a total of', length(scannedPromoters), 'TFBS, gathering results...'))
  }

  
  if (return.sequence) {
    result = TFBSTools::writeGFF3(scannedPromoters, scoreType = 'absolute')
  } else {
    result = list()
    for (x in 1:length(scannedPromoters)) {
      xx = scannedPromoters[[x]]
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
    
    
    result = purrr::by_row(y, fun, .collate = 'rows')[12:22]
  }

  ## implement later start ----
  ## scan using pairwise alignments
  # chain = import.chain('mm9.hg19.all.chain')
  # scannedPromotersPhylo = searchPairBSgenome(motif, anno, BSgenome.Hsapiens.NCBI.GRCh38, chr1="chr1", chr2="1", min.score="90%", strand="*", chain=chain)
  # test = writeGFF3(scannedPromotersPhylo)
  ## implement later end ----
  
  # return the combined dataframe
  stopifnot(length(dna) == length(query))
  
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
  if(is.null(helper)) {
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
    scannedPromoters.pval = TFBSTools::pvalues(scannedPromoters)
    pval = unlist(scannedPromoters.pval)
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

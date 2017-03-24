# function to scan promoters defined by CAGE expression of dataset of interest for TFBS

predictPromoterBinding = function(query,
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
  
  # extract DNA seq from promoters, return a DNAStringSet object
  print('fetching DNA...')
  dna = Biostrings::getSeq(anno, query)
  
  # motif.list can either be a vector with JASPAR IDs or character string 'all'
  opts = list()
  if (motif.list == 'all') {
    opts['all'] = T
  } else {
    opts[['ID']] = motif.list
  }
  print('fetching TFBS matrices...')
  motif = TFBSTools::toPWM(TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, opts))
  
  # scan the sequences!
  print('scanning the sequence...')
  scannedPromoters = TFBSTools::searchSeq(
    x = motif,
    subject = dna,
    min.score = '90%',
    mc.cores = 8L
  )
  # insert status message how many tfbs were found
  print(paste('found a total of', length(scannedPromoters), 'TFBS...'))
  
  print('gathering results...')
  
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
      unnestedCols = separate(z, start, paste0('start_', 1:maxCols), sep = ',') %>%
        separate(end, paste0('end_', 1:maxCols), sep = ',') %>%
        separate(score, paste0('score_', 1:maxCols), sep = ',') %>%
        separate(strand, paste0('strand_', 1:maxCols), sep = ',') %>%
        gather(key, value, starts_with('start'), starts_with('end'), starts_with('score'), starts_with('strand')) %>%
        mutate(value = gsub(x = value, pattern = 'c\\(', replacement = '')) %>%
        mutate(value = gsub(x = value, pattern = '\\)', replacement = '')) %>%
        mutate(value = gsub(x = value, pattern = '\"', replacement = '')) %>%
        separate(key, c('key', 'motif_hit'), sep = '_') %>%
        spread(key, value) %>%
        select(motif_ID, motif_hit, seqname, start, end, strand, everything())
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
  result2 = cbind(helper, seqnames, seqranges, strand, result)
  
  if (return.p.val) {
    print('calculating p-values...')
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
    dplyr::mutate(motif_strand = stringr::str_trim(motif_strand, side = 'both'))
  
  return(result2)
}

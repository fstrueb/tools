# function to scan promoters defined by CAGE expression of dataset of interest for TFBS

predictPromoterBinding = function(query,
                                  motif.list,
                                  input.assembly,
                                  p.val,
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
  motif = TFBSTools::toPWM(TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, opts))
  
  # scan the sequences!
  print('scanning the sequence...')
  scannedPromoters = TFBSTools::searchSeq(
    x = motif,
    subject = dna,
    min.score = '80%',
    mc.cores = 8L
  )
  print('gathering results...')
  result = TFBSTools::writeGFF3(scannedPromoters, scoreType = 'absolute')
  
  ## implement later start ----
  ## scan using pairwise alignments
  # chain = import.chain('mm9.hg19.all.chain')
  # scannedPromotersPhylo = searchPairBSgenome(motif, anno, BSgenome.Hsapiens.NCBI.GRCh38, chr1="chr1", chr2="1", min.score="90%", strand="*", chain=chain)
  # test = writeGFF3(scannedPromotersPhylo)
  ## implement later end ----
  
  # return the combined dataframe
  stopifnot(length(dna) == length(query))
  
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
  # caveat --> naming convention of column 'geneSymbol'
  helper = query[result$seqname]@elementMetadata$geneSymbol
  seqnames = query[result$seqname]@seqnames
  seqranges = as.data.frame(query[result$seqname]@ranges)
  names(seqranges) = c('query_start', 'query_end', 'query_width')
  strand = query[result$seqname]@strand
  result2 = cbind(helper, seqnames, seqranges, strand, result)
  
  if (p.val == T) {
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
    )
  
  return(result2)
}

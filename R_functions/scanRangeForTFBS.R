# function to scan promoters defined by CAGE expression of dataset of interest for TFBS
# the regulatory elements to scan for TFBS should come in a GRanges format, where every row/region of interest should correspond to exactly one regulatory element (promoter/enhancer/genomic anchor)
# source('standardizeSeqlevels.R')

scanRangeForTFBS = function(
  query,
  motif.list,
  input.assembly,
  #return.p.val = FALSE,
  #return.sequence = FALSE,
  updateProgress = NULL,
  ...) {
  options(stringsAsFactors = FALSE)
  library(magrittr)
  
  anno = vector()
  if (input.assembly == 'mm10') {
    anno = BSgenome.Mmusculus.UCSC.mm10::Mmusculus
  } else if (input.assembly == 'mm9') {
    anno = BSgenome.Mmusculus.UCSC.mm9::Mmusculus
  } else if (input.assembly == 'hg38') {
    anno = BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens
  } else {
    warning(paste('BSGenome object', input.assembly, 'is not available'))
    # anno = 'mm10'
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
  
  # need to standardize seqlevels according to anno seqlevels
  anno.seqlevels = grepl(seqlevels(anno)[1], pattern = 'chr')
  
  if (grepl(organism(anno), pattern = 'Homo')) {
    anno.species = 'human' 
  } else if (grepl(organism(anno), pattern = 'Mus')) {
    anno.species = 'mouse'
  } else {
    anno.species = 'mouse'
    warning('anno.species: nothing defined, returning mouse')
  }
  
  if (anno.seqlevels) {
    dna = Biostrings::getSeq(anno, query)
  } else {
    query = standardizeSeqlevels(query, seqlevels = anno.species, type = 'ensembl')
    dna = Biostrings::getSeq(anno, query)
  }
  
  # motif.list can either be a vector with JASPAR IDs or character string 'all'
  opts = list()
  # if (motif.list == 'all') {
  #   opts['all'] = T
  # } else {
  opts[['ID']] = motif.list
  # }
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
    text <- paste('found a total of', length(motif), 'motifs, scanning the sequence...')
    updateProgress(detail = text)
  } else {
    # extract DNA seq from promoters, return a DNAStringSet object
    print('scanning the sequence...')
  }
  # scan the sequences!
  scannedRange = TFBSTools::searchSeq(
    x = motif,
    subject = dna,
    min.score = '90%',
    mc.cores = 8L
  )
  
  scannedRange
}

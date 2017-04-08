# define range of interest

makeRange = function(chromosome, from, to) {
  result = GenomicRanges::GRanges(seqnames = chromosome, ranges = IRanges::IRanges(start = from, end = to), strand = '*')
  cat('worked')
  result
}

# makeRangeFromCsv = function(csvpath) {
#   csvfile = readr::read_csv(csvpath)
#   apply(csvfile, 1, function(x,y,z) makeRange(csvfile$Chr, csvfile$From, csvfile$To))
# }
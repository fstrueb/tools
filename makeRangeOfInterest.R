# define range of interest

makeRangeOfInterest = function(chromosome, from, to) {
  result = GenomicRanges::GRanges(seqnames = chromosome, ranges = IRanges::IRanges(start = from, end = to), strand = '*')
  result
}
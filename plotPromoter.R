# function to plot and visualize a promoter track

## xample
# test.gr = sox11InBdnf %>%
#   mutate(start = query_start+motif_start, end = query_start+motif_end, strand = motif_strand) %>%
#   makeGRangesFromDataFrame(., keep.extra.columns = T, seqnames.field = 'query_chromosome')
# chr = 'chr2'

# range = bdnf %>%
#   mutate(start = query_start+motif_start, end = query_start+motif_end, strand = motif_strand)


plotPromoter = function(range, chr, from = NULL, to = NULL) {
  ## check input
  # if (class(range) != 'GRanges') {
  #   stop('You must supply a GenomicRanges object for range')
  # }
  
  ## first remove unwanted scaffolds and convert seqnames to UCSC style
  #range = standardizeSeqlevels(range, seqlevels = 'mouse', type = 'UCSC')
  
  obj = as.data.frame(range) %>%
    dplyr::mutate(midpoint = (end - start) / 2)
  
  if (is.null(chr)) {
    stop('You must supply a chromosome')
  }
  
  if(!(grepl('chr', chr))) {
    chr = paste0('chr', chr)
  }
  
  if (is.null(from) | is.null(to)) {
    start = min(obj$start)
    end = max(obj$end)
  }
  
  # note: single hits for motif_ID are not supported
  
  bw = 'SJ' 
  n = 512
  adjust = 0.1
  
  m = c('MA0596.1', 'MA0075.1')
  
  range = range %>%
    dplyr::filter(motif_ID %in% m)
  
  nonSingles = range %>%
    dplyr::group_by(motif_ID) %>%
    dplyr::tally() %>%
    dplyr::filter(n > 1)
  
  ## constructed with smoothed linear heatmap
  # tfbsCoverage = range %>% 
  #   dplyr::filter(motif_ID %in% nonSingles$motif_ID) %>%
  #   ggplot(aes(x = start, y = factor(motif_ID)))+
  #   stat_density(aes(fill = ..density..), geom = "raster", position = "identity") +
  #   facet_wrap(~motif_strand) +
  #   scale_fill_continuous(guide = F) +
  #   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  #   xlim(start, end) +
  #   labs(y = '')
  
  ## constructed with ggbio
  ## not pretty
  # range %>%
  #   ggplot(aes(xmin = start, xmax = end+50, ymin = motif_hit, ymax = motif_hit + 4, fill = motif_name)) +
  #   geom_rect(stat = 'identity') +
  #   facet_wrap(~motif_strand)
  
  
  geneRegion = ggbio::autoplot(TxDb.Mmusculus.UCSC.mm10.knownGene, which = range(x = GRanges(seqnames = chr, IRanges(start = start, end = end)))) +
    xlim(start, end)
  geneRegion = ggplotGrob(geneRegion)
  
  #gridExtra::grid.arrange(tfbsCoverage, geneRegion, nrow = 2)
  
  plot(geneRegion)

  }
# 
# test = plotPromoter(range, 'chr12')

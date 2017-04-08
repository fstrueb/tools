# function to plot and visualize a promoter track

plotPromoter = function(range, chr, motif.selection, from = NULL, to = NULL) {
  ## check input
  # if (class(range) != 'GRanges') {
  #   stop('You must supply a GenomicRanges object for range')
  # }
  
  ## first remove unwanted scaffolds and convert seqnames to UCSC style
  #range = standardizeSeqlevels(range, seqlevels = 'mouse', type = 'UCSC')
  
  
  
  if (is.null(chr)) {
    stop('You must supply a chromosome')
  }
  
  if(!(grepl('chr', chr))) {
    chr = paste0('chr', chr)
  }
  
  
  
  # note: single hits for motif_ID are not supported
  
  # bw = 'SJ' 
  # n = 512
  # adjust = 0.1
  m = motif.selection
  
  obj = as.data.frame(range) %>%
    dplyr::mutate(plot_start = query_start+motif_start) %>%
    dplyr::mutate(plot_end = query_start+motif_end) %>%
    dplyr::filter(motif_ID %in% m)
  
  nonSingles = obj %>%
    dplyr::group_by(motif_ID) %>%
    dplyr::tally() %>%
    dplyr::filter(n > 1)
  
  
  if (is.null(from) | is.null(to)) {
    start = min(obj$query_start)
    end = max(obj$query_end)
  }
  
  # constructed with ggplot smoothed linear heatmap
tfbsCoverage = obj %>%
  # currently does not support single motif_ID results (must be more than 2 of the same motif)
    dplyr::filter(motif_ID %in% nonSingles$motif_ID) %>%
    # dplyr::filter(motif_ID %in% m) %>%
    ggplot2::ggplot(ggplot2::aes(x = plot_start, y = as.factor(motif_ID))) +
    ggplot2::stat_density(ggplot2::aes(fill = ..density..), geom = "raster", position = "identity", adjust = 0.05) +
    #stat_bin2d(aes(fill = ..count..), binwidth = c(3,1)) +
    #geom_point(stat = 'identity', size = 10) +
    ggplot2::geom_rug(sides = 't') +
    #facet_wrap(~motif_strand, nrow = 2) +
    ggplot2::scale_fill_continuous(guide = F) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::xlim(start, end) +
    ggplot2::labs(y = '', x = '')
  
  ## constructed with ggbio, not pretty
  # range %>%
  #   ggplot(aes(xmin = start, xmax = end+50, ymin = motif_hit, ymax = motif_hit + 4, fill = motif_name)) +
  #   geom_rect(stat = 'identity') +
  #   facet_wrap(~motif_strand)
  
  geneRegion = ggplot2::ggplotGrob(ggbio::autoplot(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, which = range(x = GenomicRanges::GRanges(seqnames = chr, IRanges::IRanges(start = start, end = end)))) +
    ggplot2::xlim(start, end))
  
  # import retina dnase bigwig
  dnase = rtracklayer::import.bw(con = '../ext/DNaseSeqRetina_SRX191033.bw', selection = rtracklayer::BigWigSelection(ranges = GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = start, end = end))))
  
  dnasePlot = ggplot2::ggplotGrob(ggplot2::ggplot(dnase, ggplot2::aes(x = start, y = score)) +
    ggplot2::geom_area() +
    ggplot2::labs(y = '', x = '') +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()))
  
  plot = gridExtra::grid.arrange(tfbsCoverage, geneRegion, dnasePlot, nrow = 3)
  plot

  }


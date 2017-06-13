# script to predict regulatory regions from CAGE data

library(CAGEr)
library(dplyr)
library(ggplot2)
library(stringr)
library(rtracklayer)
library(JASPAR2016)
library(TFBSTools)
library(GenomicFeatures)

getCage = function(samples, species, assembly, updateProgress, nrCores = 4) {
  
  if (species == 'mouse') {
    data("FANTOM5mouseSamples")
  } else if (species == 'human') {
    data('FANTOM5humanSamples')
  } else {
    stop('`Species` must be one of `mouse`, `human`')
  }
  
  ########### STEP 1 ###############
  if (is.function(updateProgress)) {
    text <- ('Downloading CAGE sets from FANTOM5...')
    updateProgress(detail = text)
  } else {
    print('Downloading CAGE sets from FANTOM5...')
  }
  
  # sample should return a character vector identical to the description field in the FANTOM5 data table
  cageData = importPublicData(source = 'FANTOM5', dataset = species, sample = samples) 
  
  ########### STEP 2 ###############
  if (is.function(updateProgress)) {
    text <- ('Clustering consensus transcription start sites, this can take some time...')
    updateProgress(detail = text)
  } else {
    print('Clustering consensus transcription start sites, this can take some time...')
  }
  
  if (length(samples) > 1) {
    normalizeTagCount(cageData, method = 'simpleTpm')
    clusterCTSS(cageData, method = 'distclu', maxDist = 10, useMulticore = T, nrCores = nrCores)
    aggregateTagClusters(cageData, tpmThreshold = 0.001, excludeSignalBelowThreshold = F, maxDist = 50)
    clusters = cageData@consensusClusters
  } else { 
    normalizeTagCount(cageData, method = 'simpleTpm')
    clusterCTSS(cageData, method = 'distclu', maxDist = 10, useMulticore = T, nrCores = nrCores)
    clusters = tagClusters(cageData, samples)
  }
  
  # rename cluster column to consensus cluster if input sample length was 1
  if (length(samples) == 1) {
      names(clusters)[1] = 'consensus.cluster'
    }
  
  
  ########### STEP 3 ###############
  if (is.function(updateProgress)) {
    text <- ('Clustering finished, predicting enhancers...')
    updateProgress(detail = text)
  } else {
    print('Clustering finished, predicting enhancers...')
  }
  
  clusters.gr = makeGRangesFromDataFrame(clusters, keep.extra.columns = T)
  cat('clusters present, length', dim(clusters)[1], '\n')
  
  # filter putative enhancers with bidirectional transcription in 400 bp window (see FANTOM5 CAGE paper, Nature 2012)
  # first remove overlapping antisense tss
  antisenseTss = findOverlaps(clusters.gr, drop.self = T, ignore.strand = T)
  clusters.gr = clusters.gr[-(queryHits(antisenseTss))]
  clusters.gr = clusters.gr[-(subjectHits(antisenseTss))]
  
  # separate strands, then find nearest cluster
  clusters.gr.minus = clusters.gr[strand(clusters.gr) == '-']
  clusters.gr.plus = clusters.gr[strand(clusters.gr) == '+']
  nearestCluster = distanceToNearest(clusters.gr.minus, clusters.gr.plus, select = 'arbitrary', ignore.strand = T)
  # filter for transcribed pairs separated by < 400 bp
  filteredClusters = nearestCluster[nearestCluster@elementMetadata$distance <= 400]
  # merge bidirectional clusters
  bidirectClusters = GenomicRangesList()
  for (i in 1:length(filteredClusters)) {
    bidirectClusters[[i]] = c(clusters.gr.minus[queryHits(filteredClusters[i])], clusters.gr.plus[subjectHits(filteredClusters[i])])
  }
  
  
  # calculate directionality score, return GRangesList object
  # D=(F-R)/ (F+R)
  biClu = GenomicRangesList()
  a = 1
  for (i in 1:length(bidirectClusters)) {
    
    if (start(bidirectClusters[[i]][1]) < start(bidirectClusters[[i]][2])) {
      biClu[[a]] = GRanges(
        seqnames = seqnames(bidirectClusters[[i]][1]),
        ranges = IRanges(
          start = min(start(bidirectClusters[[i]]@ranges)),
          end = max(end(bidirectClusters[[i]]@ranges))
        ),
        strand = '*',
        DataFrame(directionalityScore = (bidirectClusters[[i]]@elementMetadata$tpm[2] - bidirectClusters[[i]]@elementMetadata$tpm[1]) / (bidirectClusters[[i]]@elementMetadata$tpm[2] + bidirectClusters[[i]]@elementMetadata$tpm[1]), 
          origClusterMinus = bidirectClusters[[i]]@elementMetadata$consensus.cluster[1],
          origClusterPlus = bidirectClusters[[i]]@elementMetadata$consensus.cluster[2],
          origExprMinus = bidirectClusters[[i]]@elementMetadata$tpm[1],
          origExprPlus = bidirectClusters[[i]]@elementMetadata$tpm[2]
        )
      )
      a = a + 1
    } 
  }
  
  
  enhancers = unlist(biClu)
  enhancers = enhancers[abs(enhancers@elementMetadata$directionalityScore) <= 0.8]
  
  ########### STEP 4 ###############
  if (is.function(updateProgress)) {
    text <- ('Annotating clusters...')
    updateProgress(detail = text)
  } else {
    print('Annotating clusters...')
  }
  
  # find nearest genes for annotation
  if (species == 'mouse') {
    transcripts = exons(TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene)
  } else if (species == 'human') {
    # !!!!CHECK ANNOTATION FOR CORRECT CAGE VERSION!!!! (in human)
    # transcripts = transcripts(TxTxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
  } else {
    warning('you messed up')
  }
  
  enhDist = distanceToNearest(enhancers, transcripts, ignore.strand = T)
  enhFil = enhDist[enhDist@elementMetadata$distance > 1]
  enh = enhancers[queryHits(enhFil)]
  enh@elementMetadata$distToNearest = enhDist[queryHits(enhFil)]@elementMetadata$distance
  enhaMinus = data.frame(enh@elementMetadata$origClusterMinus, enh@elementMetadata$directionalityScore, enh@elementMetadata$distToNearest, enh@elementMetadata$origClusterPlus)
  enhaPlus = data.frame(enh@elementMetadata$origClusterPlus, enh@elementMetadata$origClusterMinus, enh@elementMetadata$directionalityScore, enh@elementMetadata$distToNearest)
  
  clustersResult = left_join(clusters, enhaMinus, by = c('consensus.cluster' = 'enh.elementMetadata.origClusterMinus')) %>%
    left_join(., enhaPlus, by = c('consensus.cluster' = 'enh.elementMetadata.origClusterPlus')) %>%
    dplyr::mutate(bidirectional = ifelse(is.na(enh.elementMetadata.directionalityScore.x) & is.na(enh.elementMetadata.directionalityScore.y), F, T)) %>%
    tidyr::unite(oppositeCluster, enh.elementMetadata.origClusterPlus, enh.elementMetadata.origClusterMinus, sep = '') %>%
    dplyr::mutate(oppositeCluster = as.numeric(gsub(oppositeCluster, pattern = 'NA', replacement = ''))) %>%
    tidyr::unite(distanceToNearestExon, enh.elementMetadata.distToNearest.x, enh.elementMetadata.distToNearest.y, sep = '') %>%
    dplyr::mutate(distanceToNearestExon = as.numeric(gsub(distanceToNearestExon, pattern = 'NA', replacement = ''))) %>%
    tidyr::unite(directionalityScore, enh.elementMetadata.directionalityScore.x, enh.elementMetadata.directionalityScore.y, sep = '') %>%
    dplyr::mutate(directionalityScore = as.numeric(gsub(directionalityScore, pattern = 'NA', replacement = '')))
  
  clustersResult.gr = makeGRangesFromDataFrame(clustersResult, keep.extra.columns = T)
  
  # find closest gene symbol
  # extract all transcripts from TxDb
  transcriptsMm9 = transcripts(TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene)
  transcriptsMm9Df = as_data_frame(transcriptsMm9) %>%
    dplyr::mutate(rowNum = row_number()) %>%
    dplyr::select(tx_name, rowNum)
  # # find the nearest transcript to the CTSSs
  # allCountsFilteredGr = makeGRangesFromDataFrame(allCountsFilteredDf, keep.extra.columns = T)
  nearestGene = nearest(clustersResult.gr, transcriptsMm9)
  # import mapping table for UCSC gene identifiers to standard symbols
  kgxref = read.csv('../ext/kgxref_mm9.csv', header = T, stringsAsFactors = F, sep = ';')
  kgxrefLight = dplyr::select(kgxref, kgID, geneSymbol)
  # sanity check
  length(nearestGene) == dim(clustersResult)[1]
  # annotate CTSSs with closest gene, respects strand
  clustersResult$closestGene = nearestGene
  helper1 = inner_join(clustersResult, transcriptsMm9Df, by = c('closestGene' = 'rowNum'))
  helper2 = left_join(helper1, kgxrefLight, by = c('tx_name' = 'kgID'), copy = T)
    # dplyr::select(1:10, 13)
  result = helper2
  result = result %>%
    dplyr::rename(closestGeneSymbol = geneSymbol, consensusClusterID = consensus.cluster) %>%
    dplyr::select(chr, start, end, strand, tpm, closestGeneSymbol, bidirectional, directionalityScore, distanceToNearestExon, oppositeCluster, consensusClusterID)
  # remove junk variables
  rm(helper1, helper2)
  
  result
}



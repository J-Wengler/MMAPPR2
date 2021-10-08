#' Generate potential causative mutations and consequences in peak regions
#'
#' Follows the \code{\link{peakRefinement}} step and produces a
#' \code{\linkS4class{MmapprData}} object ready for
#' \code{\link{outputMmapprData}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{candidates}
#'   slot filled with a \code{\link[GenomicRanges]{GRanges}} object for each
#'   peak chromosome containing variants and predicted consequences from
#'   Ensembl's Variant Effect Predictor.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)
#'         & all(Sys.which(c("samtools", "vep")) != "")) {
#'     mmappr_param <- MmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
#' \dontrun{
#' md <- new('MmapprData', param = mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#'
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' }

generateCandidates <- function(mmapprData) {
  
  saveRDS(mmapprData, file="mmappr_data_1.RDS") #testing
  #get GRanges representation of peak
  peakGRanges <- lapply(mmapprData@peaks, .getPeakRange)
  
  saveRDS(peakGRanges, file="peakGranges.RDS") #testing
  
  #call variants in peak
  mmapprData@candidates$snps <- lapply(peakGRanges,
                                       FUN=.getVariantsForRange,
                                       param=mmapprData@param)
  saveRDS(mmapprData, file="mmappr_data_2.RDS") #testing
  
  #remove NULL values from snps
  allPeaks <- length(mmapprData@candidates$snps)
  mmapprData@candidates$snps <- mmapprData@candidates$snps[!sapply(mmapprData@candidates$snps,is.null)]
  noNullPeaks <- length(mmapprData@candidates$snps)
  
  if(allPeaks != noNullPeaks) {
    peaksRemoved <- allPeaks - noNullPeaks
    .messageAndLog("Warning: peaks without valid snps were called but then removed", oF)
    .messageAndLog(paste("Number of peaks removed:" , peaksRemoved), oF)
  }
  
  
  
  #run VEP
  mmapprData@candidates$effects <- lapply(mmapprData@candidates$snps,
                                          FUN=.runVEPForVariants,
                                          param=mmapprData@param)
  saveRDS(mmapprData, file="mmappr_data_3.RDS") #testing
  
  #filter out low impact effects
  mmapprData@candidates$effects <- lapply(mmapprData@candidates$effects, 
                                          FUN=.filterVariants, 
                                          impact = mmapprData@param@vepImpact)
  
  #add diff expressed genes
  mmapprData@candidates$diff <- lapply(peakGRanges,
                                       FUN=.addDiff,
                                       param=mmapprData@param)
  
  #density score and order variants
  mmapprData@candidates <- lapply(mmapprData@candidates,
                                  FUN=.scoreVariants,
                                  mmapprData@peaks)
  
  print("Success: generateCandidates") #debug
  
  return(mmapprData)
}


.getPeakRange <- function(peakList) {
  ir <- IRanges::IRanges(start=as.numeric(peakList$start),
                         end=as.numeric(peakList$end),
                         names=peakList$seqname)
  
  gr <- GRanges(seqnames=names(ir),
                ranges=ir)
  print("Success: .getPeakRange") #debug
  return(gr)
}


.getVariantsForRange <- function(inputRange, param) {
  # merge files in desired region if there are multiple
  print("Start getVariantsForRange")
  mergedBam <- file.path(outputFolder(param), 'merged.tmp.bam')
  print("mergedBam")
  if (length(param@mutFiles) < 2) {
    print("length less than 2")
    mutBam <- param@mutFiles[[1]]
    print("got mutBam")
  }else{
    mutBam <- mergeBam(param@mutFiles,
                       destination=mergedBam,
                       region=inputRange)
    print("mergedBam")
  }
  
  # create param for variant calling
  print("Creating tallyParam")
  print(inputRange)
  tallyParam <- TallyVariantsParam(genome=param@refGenome,
                                   which=inputRange,
                                   indels=TRUE)
  print("tallyParam created")
  
  resultVr <- callVariants(mutBam, tally.param=tallyParam)
  print("callVariants complete")
  resultVr <- resultVr[altDepth(resultVr)/totalDepth(resultVr) > 0.8]
  print("filtered variants >.8")
  if (file.exists(mergedBam)) file.remove(mergedBam)
  
  print("removed Bam")
  if (length(resultVr) > 0) {
    # need sampleNames to convert to VCF; using mutant file names
    print("length resultVr > 0")
    Biobase::sampleNames(resultVr) <-
      paste0(names(param@mutFiles),
             collapse = " -- ")
    
    print("samplenames done")
    S4Vectors::mcols(resultVr) <- NULL
    print("Success: .getVariantsForRange") #debug
    return(resultVr)
  }
  else {    #originally return(null) on same line
    print(".getVariantsForRange returned NULL") #debug
    #return(NULL)
    return()
  }
}


.runVEPForVariants <- function(inputVariants, param){
  vepFlags <- vepFlags(param)
  stopifnot(is(vepFlags, "VEPFlags"))
  #stopifnot(is(inputVariants, 'VRanges'))
  #check that there are variants
  if(!length(inputVariants) == 0) {
  
  vcf <- file.path(outputFolder(param), paste0(seqnames(inputVariants)[1], 'peak.vcf'))

  tryCatch({
    writeVcf(inputVariants, vcf)
    resultGRanges <- ensemblVEP(vcf, vepFlags)
  }, error=function(e) {
    stop(e)
  },finally={
    # if (file.exists(vcf)) file.remove(vcf)
  })
  
  print("Success: .runVEPForVariants") #debug
  return(resultGRanges)
  }
  print("No variants")
}


.filterVariants <- function(candidateGRanges, impact) {
  impactLevels <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  filter <- mcols(candidateGRanges)$IMPACT %in% impactLevels[1:impact]
  filter[is.na(filter)] <- TRUE
  print("Success: .filterVariants") #debug
  return(candidateGRanges[filter])
}


.addDiff <- function(peakGRange, param) {
  #prep data
  suppressMessages(genes <- fread(cmd=paste("zcat", param@gtf), 
                                  showProgress = FALSE))
  genes <- genes[V3 == "gene"
               ][, gene_id := gsub(".*gene_id \"(.*?)\";.*", "\\1", V9)
               ][, gene_name := gsub(".*gene_name \"(.*?)\";.*", "\\1", V9)
               ][, .(seqnames = V1, start = V4, end = V5, strand = V7, 
                     gene_id, gene_name)]
  genes <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
  genes <- subsetByOverlaps(x=genes, ranges = peakGRange)
  
  #get and process counts
  readfiles <- c(param@wtFiles, param@mutFiles)
  counts <- GenomicAlignments::summarizeOverlaps(features=genes,
                                                 reads=readfiles,
                                                 param=ScanBamParam(which = peakGRange))
  num_wt <- length(param@wtFiles)
  countDF <- as.data.table(assays(counts)$counts)
  countDF[, ave_wt := rowMeans(countDF[, 1:(num_wt + 1)])    
        ][, ave_mt := rowMeans(countDF[, (num_wt + 1):length(countDF)])
        ][, log2FC := round(log2(ave_mt/ave_wt), 3)]
  
  #merge counts with genes
  mcols(genes) <- cbind(mcols(genes), countDF)
  
  #filter for differentially expressed genes
  print("Success: .addDiff") #debug
  return(genes[(abs(genes$log2FC) > 1 | is.na(genes$log2FC)) 
               & (genes$ave_wt > 10 | genes$ave_mt > 10)])
}


.scoreVariants <- function(candList, peaks) {
  returnData = list()
  for (seqname in names(candList)) {
    #density calculation
    densityFunction <- peaks[[seqname]]$densityFunction
    stopifnot(!is.null(densityFunction))
    positions <- BiocGenerics::start(candList[[seqname]]) +
      ((BiocGenerics::width(candList[[seqname]]) - 1) / 2)
    densityCol <- vapply(positions, densityFunction, FUN.VALUE=numeric(1))
    mcols(candList[[seqname]])$peakDensity <- densityCol
    
    #re-order
    returnData[[seqname]] <-
      .orderVariants(candList[[seqname]], densityCol)
  }
  print("Success: .scoreVariants") #debug
  return(returnData)
}


.orderVariants <- function(candidateGRanges, densityCol) {
  if(!(class(candidateGRanges) %in% c("VRanges", "GRanges"))) { 
    print("Success: .orderVariants") #debug
    return(candidateGRanges)
  }
  if(!is.null(candidateGRanges$IMPACT)) {
    impactLevels <- c("MODIFIER", "LOW", "MODERATE", "HIGH")
    orderVec <- order(match(candidateGRanges$IMPACT, impactLevels), 
                      densityCol, decreasing=TRUE)
  } else {
    orderVec <- order(densityCol, decreasing=TRUE)
  }
  candidateGRanges <- candidateGRanges[orderVec]
  print("Success: .orderVariants") #debug
  return(candidateGRanges)
}


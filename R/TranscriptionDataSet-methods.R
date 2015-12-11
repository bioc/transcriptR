####------------------- constructTDS -------------------####
#' constructTDS
#'
#' The function constructs an object of class \code{\link{TranscriptionDataSet}},
#' which is a container for holding processed sequencing data and the results of
#' all downstream analyses. All the slots of the created object are filled
#' during the workflow by applying specific functions to the object directly.
#'
#' @name constructTDS
#'
#' @param file A path to a BAM file with sequencing reads.
#' @param region \code{\link[GenomicRanges]{GRanges}}. Genomic region(s) to
#'     extract reads from. If not supplied, all the reads from a BAM file are
#'     extracted.
#' @param fragment.size \code{Numeric}. Extend read length to the fragment size.
#'     Default: 250.
#' @param unique \code{Logical}. Whether to remove duplicated reads (based on the
#'     genomic coordinates). Default: FALSE.
#' @param paired.end \code{Logical}. Whether to treat a BAM file as paired-end.
#'     Default: FALSE.
#' @param swap.strand \code{Logical}. Whether to reverse the strand of the read.
#'     Default: FALSE.
#' @param param \code{\link[Rsamtools]{ScanBamParam}} object influencing what
#'     fields and which records (reads) are imported from a BAM file. Default:
#'     NULL.
#'
#' @details
#'     The slots \code{fragments}, \code{fragmentSize}, \code{region},
#'     \code{coveragePlus}, \code{coverageMinus} are filled during the object
#'     construction. The \code{fragments} holds information about genomic
#'     coordinates of the sequenced fragments (reads extended to the fragmento
#'     size). \code{coveragePlus} and \code{coverageMinus} for each position in
#'     the genome counts the number of fragments that cover it (for the details,
#'     see \code{\link[IRanges]{coverage}}). \code{region} holds information about
#'     the region used for fragments extraction.
#'
#' @return An object of class \code{\link{TranscriptionDataSet}}.
#'
#' @author Armen R. Karapetyan
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### View a short summary of the object
#' tds
#'
#' @export
constructTDS <- function(file, region, fragment.size = 250, unique = FALSE,
                         paired.end = FALSE, swap.strand = FALSE,
                         param = NULL){

        if (fragment.size > 500) {
                warning("'fragment.size' seems to be rather large. Please, verify that the 'fragment.size' is set correctly.")
        }

        message("[INFO] Extracting fragments...", appendLF = FALSE)
        suppressWarnings(fragments <- .readBam(file, region, fragment.size, unique,
                                               paired.end, swap.strand, param))
        # trim out-of-bound ranges
        fragments <- trim(fragments)

        message("Done!")

        # Create coverage profile for each strand
        message("[INFO] Creating coverage profile...", appendLF = FALSE)

        cov.plus <- coverage(fragments[strand(fragments) == "+"])
        cov.minus <- coverage(fragments[strand(fragments) == "-"])

        # Discard empy chromosomes
        cov.plus <- cov.plus[which(unlist(lapply(cov.plus, function(x) length(S4Vectors::runValue(x)))) > 1)]
        cov.minus <- cov.minus[which(unlist(lapply(cov.minus, function(x) length(S4Vectors::runValue(x)))) > 1)]

        message("Done!")
        if (missing(region)){
                region <- GRanges()
        }

        return(new("TranscriptionDataSet", bamFile = file, fragments = fragments,
                   fragmentSize = fragment.size, region = region,
                   coveragePlus = cov.plus, coverageMinus = cov.minus))
}


####------------------- estimateBackground -------------------####
#' @rdname estimateBackground-methods
setMethod("estimateBackground",
          signature = signature("TranscriptionDataSet"),
          definition = function(object, fdr.cutoff = 0.05){

                  if (length(object@coveragePlus) == 0 | length(object@coverageMinus) == 0){
                          stop("[ERROR] 'object' doesn't contain coverage information.
                               Please, use constructTDS function to initialize an object of class 'TranscriptionDataSet'.")
                  }

                  if (is.numeric(fdr.cutoff)){
                          if (fdr.cutoff <= 0 | fdr.cutoff > 1){
                                  stop("[ERROR] 'fdr.cutoff' needs to be in a range (0, 1].")
                          }
                  } else {
                          stop("[ERROR] 'fdr.cutoff' needs to be numeric.")
                  }

                  objName <- deparse(substitute(object))
                  cutoff <- mean(unlist(lapply(list(object@coveragePlus, object@coverageMinus),
                                               function(x){
                                                       chipseq::peakCutoff(cov = x, fdr.cutoff = fdr.cutoff)
                                               })))
                  object@coverageCutoff <- round(cutoff, 3)
                  object@coverageCutoffFdr <- fdr.cutoff
                  assign(objName, object, envir = parent.frame())
          }
)


####------------------- estimateGapDistance -------------------####
.dissectTranscribedRegions <- function(cov, coverage.cutoff, strand){
        irl <- as(slice(cov, lower = coverage.cutoff), "IRangesList")
        trr <- as(irl, "GRanges")
        strand(trr) <- strand
        return(trr)
}

.rpkm <- function(C, S, N){
        (10^9 * C) / (as.numeric(N) * S)
}

.estimateFpkm <- function(object, reads, total.reads){
        # count number of reads per transcripts
        counts <- countOverlaps(query = object,
                                subject = reads,
                                ignore.strand = FALSE)
        object$fragments <- counts
        # estiamte rpkm/fpkm
        object$fpkm <- round(.rpkm(C = counts, S = width(object), N = total.reads), 3)
        return(object)
}

.estimateErrorRate <- function(trx, annot){
        trx.per.annot <- countOverlaps(annot, trx, ignore.strand = FALSE)
        annot.per.trx <- countOverlaps(trx, annot, ignore.strand = FALSE)
        error.dissected <- sum(trx.per.annot >= 2)/sum(trx.per.annot >= 1)
        error.merged <- sum(annot.per.trx >= 2)/sum(annot.per.trx >= 1)
        sum.two.errors <- error.dissected + error.merged
        error.df <- data.frame(error.dissected, error.merged, sum.two.errors)
}

#' @rdname estimateGapDistance-methods
setMethod("estimateGapDistance",
          signature = signature("TranscriptionDataSet", "GRanges"),
          definition = function(object,
                                annot,
                                coverage.cutoff,
                                filter.annot = TRUE,
                                fpkm.quantile = 0.25,
                                gap.dist.range = seq(from = 0, to = 10000, by = 100)){

                  # Tests
                  if (length(object@coveragePlus) == 0 | length(object@coverageMinus) == 0){
                          stop("[ERROR] 'object' doesn't contain coverage information.
                               Please, use constructTDS function to initialize an object of class 'TranscriptionDataSet'.")
                  }

                  if (missing(coverage.cutoff)){
                          if (length(object@coverageCutoff) == 0){
                                  stop("[ERROR] 'coverageCutoff' slot is empty.
                                       Please, specify 'coverage.cutoff' manually or call estimateBackground function
                                       on the 'object'.")
                          } else {
                                  coverage.cutoff <- object@coverageCutoff
                          }
                  } else {
                          if (is.numeric(coverage.cutoff)){
                                  if (coverage.cutoff <= 0) {
                                          stop("[ERROR] 'coverage.cutoff' needs to be above 0.")
                                  }
                          } else {
                                  stop("[ERROR] 'coverage.cutoff' needs to be numeric.")
                          }
                  }

                  if (!is.logical(filter.annot)){
                          stop("[ERROR] 'filter.annot' needs to be a boolean.")
                  }

                  if (is.numeric(fpkm.quantile)){
                          if (!(fpkm.quantile > 0 & fpkm.quantile < 1)){
                                  stop("[ERROR] 'fpkm.quantile' needs to be in a range (0, 1).")
                          }
                  } else {
                          stop("[ERROR] 'fpkm.quantile' needs to be numeric.")
                  }

                  if (is.numeric(gap.dist.range)){
                          if (any(gap.dist.range < 0)){
                                  stop("[ERROR] 'gap.dist.range' needs to be a vector of positive numbers.")
                          }
                  } else{
                          stop("[ERROR] 'gap.dist.range' needs to be a vector of positive numbers.")
                  }


                  objName <- deparse(substitute(object))

                  # Slice regions exceeding cutoff
                  message("[INFO] Dissecting transcribed regions...", appendLF = FALSE)

                  trr.plus <- .dissectTranscribedRegions(object@coveragePlus, coverage.cutoff = coverage.cutoff, strand = "+")
                  trr.minus <- .dissectTranscribedRegions(object@coverageMinus, coverage.cutoff = coverage.cutoff, strand = "-")
                  trr <- c(trr.plus, trr.minus)

                  message("Done!")

                  # Estimate gap distance minimazing the sum of two errors
                  message("[INFO] Estimating gap distance minimizing sum of two errors...", appendLF = FALSE)

                  # Keep annot. overlapping region
                  if (length(object@region) > 0){
                          annot <- subsetByOverlaps(query = annot,
                                                    subject = object@region,
                                                    ignore.strand = TRUE)
                  }

                  # Merge overlaping annot on the same strand
                  annot <- reduce(annot, ignore.strand = FALSE)

                  if (filter.annot){
                          # Estimate expression level of annotations
                          # Discard lowly expressed annotations
                          annot.stats <- .estimateFpkm(annot, object@fragments, total.reads = length(object@fragments))
                          fpkm.quant <- stats::quantile(x = annot.stats$fpkm, probs = fpkm.quantile)
                          annot <- annot[which(annot.stats$fpkm  > fpkm.quant)]
                  }

                  # data.frame to store output values
                  error.df <- data.frame()

                  for (i in seq_along(gap.dist.range)){
                          trx <- GenomicRanges::reduce(trr, min.gapwidth = gap.dist.range[i], ignore.strand = FALSE)
                          error <- .estimateErrorRate(trx, annot)
                          error$gap.distance <- gap.dist.range[i]
                          error.df <- rbind(error.df, error)
                  }
                  message("Done!")

                  object@gapDistanceTest <- error.df
                  object@gapDistanceTestCovCutoff <- coverage.cutoff
                  assign(objName, object, envir = parent.frame())
          }
)


####------------------- detectTranscripts -------------------####
.estimateDensity <- function(object, cov.plus, cov.minus, coverage.cutoff){
        chrs <- levels(droplevels(seqnames(object)))
        object.per.chr <- lapply(chrs,
                                 function(chr){
                                         # split by chromosome and strand
                                         object.plus <- object[seqnames(object) == chr & strand(object) == "+"]
                                         object.minus <- object[seqnames(object) == chr & strand(object) == "-"]

                                         # covearge profile per transcript
                                         views.plus <- Views(subject = cov.plus[[chr]],
                                                             start = start(object.plus),
                                                             end = end(object.plus))

                                         views.minus <- Views(subject = cov.minus[[chr]],
                                                              start = start(object.minus),
                                                              end = end(object.minus))

                                         # number of bases covered by fragments
                                         object.plus$bases.covered <- viewApply(views.plus,
                                                                                function(x){
                                                                                        sum(x > coverage.cutoff)
                                                                                }
                                         )

                                         object.minus$bases.covered <- viewApply(views.minus,
                                                                                 function(x){
                                                                                         sum(x > coverage.cutoff)
                                                                                 }
                                         )

                                         object.chr <- c(object.plus, object.minus)
                                         object.chr$coverage <- round(object.chr$bases.covered / width(object.chr), 3)
                                         object.chr
                                 }
        )

        object.comb <- do.call(c, object.per.chr)
        return(object.comb)
}

#' @rdname detectTranscripts-methods
setMethod("detectTranscripts",
          signature = signature("TranscriptionDataSet"),
          definition = function(object,
                                coverage.cutoff,
                                gap.dist,
                                estimate.params = TRUE,
                                total.reads,
                                combine.by.annot = FALSE,
                                annot){

                  if (length(object@coveragePlus) == 0 | length(object@coverageMinus) == 0){
                          stop("[ERROR] 'object' doesn't contain coverage information.
                               Please, use constructTDS function to initialize an object of class 'TranscriptionDataSet'.")
                  }

                  if (missing(coverage.cutoff)){
                          if (length(object@coverageCutoff) == 0){
                                  stop("[ERROR] 'coverageCutoff' slot is empty.
                                       Please, specify 'coverage.cutoff' manually or call estimateBackground function
                                       on the 'object'.")
                          } else {
                                  coverage.cutoff <- object@coverageCutoff
                          }
                  } else {
                          if (is.numeric(coverage.cutoff)){
                                  if (coverage.cutoff <= 0) {
                                          stop("[ERROR] 'coverage.cutoff' needs to be above 0.")
                                  }
                          } else {
                                  stop("[ERROR] 'coverage.cutoff' needs to be numeric.")
                          }
                  }

                  if (missing(gap.dist)){
                          if (length(object@gapDistanceTest) == 0){
                                  stop("[ERROR] 'gapDistanceTest' slot is empty.
                                       Please, specify 'gap.dist' manually or call estimateGapDistance function
                                       on the 'object'.")
                          } else {
                                  idx <- which.min(object@gapDistanceTest$sum.two.errors)
                                  gap.dist <- object@gapDistanceTest$gap.distance[idx]
                                  sum.two.errors <- object@gapDistanceTest$sum.two.errors[idx]
                          }
                  } else {
                          if (is.numeric(gap.dist)){
                                  if (gap.dist < 0) {
                                          stop("[ERROR] 'gap.dist' needs to be a positive number.")
                                  }
                          } else {
                                  stop("[ERROR] 'gap.dist' needs to be numeric.")
                          }
                  }

                  if (missing(total.reads)){
                          total.reads <- length(object@fragments)
                  } else {
                          if (is.numeric(total.reads)){
                                  if (total.reads < 0){
                                          stop("[ERROR] 'total.reads' needs to be a positive number.")
                                  }
                          } else {
                                  stop("[ERROR] 'total.reads' needs to be numeric.")
                          }
                  }

                  if (combine.by.annot){
                          if (missing(annot)){
                                  stop("[ERROR] Please, provide annotations.")
                          } else {
                                  if(!is(annot, "GRanges")){
                                          stop("[ERROR] 'annot' needs to be a GRanges object.")
                                  }
                          }
                  }

                  objName <- deparse(substitute(object))

                  # Slice regions exceeding coverage cutoff
                  message("[INFO] Dissecting transcripts...", appendLF = FALSE)

                  trr.plus <- .dissectTranscribedRegions(object@coveragePlus,
                                                         coverage.cutoff = coverage.cutoff,
                                                         strand = "+")
                  trr.minus <- .dissectTranscribedRegions(object@coverageMinus,
                                                          coverage.cutoff = coverage.cutoff,
                                                          strand = "-")
                  trx <- GenomicRanges::reduce(c(trr.plus, trr.minus),
                                               min.gapwidth = gap.dist,
                                               ignore.strand = FALSE)

                  message("Done!")

                  if (combine.by.annot){
                          message("[INFO] Merging transcripts by the annotation overlap...", appendLF = FALSE)

                          # overlap betweem annotations and transcripts
                          hits <- findOverlaps(query = reduce(annot, ignore.strand = FALSE),
                                               subject = trx, ignore.strand = FALSE)

                          # split transcripts by annotation overlap
                          trx.by.annot <- split(subjectHits(hits), queryHits(hits))
                          trx.per.annot <- lapply(trx.by.annot, length)

                          # select annotations overlaping multiple transcripts
                          annot.mult.trx <- trx.by.annot[trx.per.annot > 1]

                          # Merge transcripts overlaping the same annotation
                          trx.comb <- lapply(annot.mult.trx,
                                             function(x){
                                                     trx.per.annot <- trx[x]
                                                     s <- min(start(trx.per.annot))
                                                     e <- max(end(trx.per.annot))
                                                     GRanges(seqnames = as.character(S4Vectors::runValue(seqnames(trx.per.annot))),
                                                             ranges = IRanges(start = s, end = e),
                                                             strand = as.character(S4Vectors::runValue(strand(trx.per.annot))))
                                             }
                          )
                          trx.comb <- suppressWarnings(BiocGenerics::Reduce(c, trx.comb))

                          # Discard duplicated transcripts
                          trx <- trx[-unlist(annot.mult.trx)]

                          # Combine merged transcritps with the rest of the transcripts
                          trx <- c(trx.comb, trx)
                          trx <- sort.GenomicRanges(trx, decreasing = FALSE, ignore.strand = TRUE)

                          message("Done!")
                  }

                  trx <- GenomicRanges::sort(x = trx, decreasing = FALSE, ignore.strand = TRUE)
                  trx$id <- paste0("trx_", seq_along(trx))
                  trx$length <- width(trx)

                  if (estimate.params){
                          # Calculate transcripts coverage (proportion of bases covered by fragments/reads) and fpkm
                          message("[INFO] Calculating fpkm and coverage...", appendLF = FALSE)
                          trx <- .estimateDensity(trx, object@coveragePlus, object@coverageMinus, coverage.cutoff)
                          trx <- .estimateFpkm(trx, object@fragments, total.reads)
                          trx <- GenomicRanges::sort(x = trx, decreasing = FALSE, ignore.strand = TRUE)
                          message("Done!")
                  }

                  object@transcripts <- trx
                  object@transcriptsCovCutoff <- coverage.cutoff
                  object@transcriptsGapDist <- gap.dist
                  object@transcriptsNormalization <- total.reads
                  assign(objName, object, envir = parent.frame())
          }
)


####------------------- breakTranscriptsByPeaks -------------------####
#' @rdname breakTranscriptsByPeaks-methods
setMethod("breakTranscriptsByPeaks",
          signature = signature("TranscriptionDataSet", "ChipDataSet"),
          definition = function(tdsObj, cdsObj, estimate.params = TRUE){

                  objName <- deparse(substitute(tdsObj))
                  msg.alarm.plus <- FALSE
                  msg.alarm.minus <- FALSE

                  # Prepare data
                  trx <- tdsObj@transcripts
                  peaks <- cdsObj@peaks
                  predicted.tssOverlap <- cdsObj@tssOverlapPrediction$predicted.tssOverlap$predicted.tssOverlap
                  predicted.strand <- cdsObj@strandPrediction$predicted.strand
                  trans.start.pos.plus <- cdsObj@strandPrediction$results.plus$q1q2.sepline.coord
                  trans.start.pos.minus <- cdsObj@strandPrediction$results.minus$q1q2.sepline.coord
                  values(peaks) <- cbind(values(peaks), S4Vectors::DataFrame(predicted.tssOverlap,
                                                                             trans.start.pos.plus,
                                                                             trans.start.pos.minus,
                                                                             predicted.strand))

                  # Select peaks predicted to be gene associated
                  peaks <- peaks[ peaks$predicted.tssOverlap == "yes"]

                  # Split transcripts by strand
                  trx.plus <- trx[strand(trx) == "+"]
                  trx.minus <- trx[strand(trx) == "-"]

                  # Split peaks by the predicted strand
                  peaks.plus <- peaks[ peaks$predicted.strand == "+" | peaks$predicted.strand == "bi"]
                  peaks.minus <- peaks[ peaks$predicted.strand == "-" | peaks$predicted.strand == "bi"]

                  message("[IFNO] Breaking transcripts by peaks...", appendLF = FALSE)

                  # Process plus strand
                  # Select peaks overlaping with the transcripts starts
                  hits.plus <- findOverlaps(promoters(trx.plus, upstream = 1, downstream = 1), peaks.plus, ignore.strand = TRUE)
                  # Select peaks not ovelaping with the transcripts starts (located in the transcript body)
                  idx.plus <- setdiff(seq_along(peaks.plus), unique(subjectHits(hits.plus)))
                  if (length(idx.plus) == 0){
                          msg.alarm.plus <- TRUE
                  } else {
                          # Prepare break region for plus strand
                          s <- peaks.plus[idx.plus]$trans.start.pos.plus - 1
                          e <- peaks.plus[idx.plus]$trans.start.pos.plus - 1
                          chr <- seqnames(peaks.plus[idx.plus])
                          breaks.plus <- GRanges(seqnames = chr, ranges = IRanges(start = s, end = e), strand = "+")
                          # Truncate transcritps by transcription start position inside the peak (estimate by predictStrand function)
                          trx.plus <- setdiff(trx.plus, breaks.plus)
                  }

                  # Process minus strand
                  # Select peaks overlaping with the transcripts starts
                  hits.minus <- findOverlaps(promoters(trx.minus, upstream = 1, downstream = 1), peaks.minus, ignore.strand = TRUE)
                  # Select peaks not ovelaping with the transcripts starts (located in the transcript body)
                  idx.minus <- setdiff(seq_along(peaks.minus), unique(subjectHits(hits.minus)))
                  if (length(idx.minus) == 0){
                          msg.alarm.minus <- TRUE
                  } else {
                          # Prepare break region for plus strand
                          s <- peaks.minus[idx.minus]$trans.start.pos.minus + 1
                          e <- peaks.minus[idx.minus]$trans.start.pos.minus + 1
                          chr <- seqnames(peaks.minus[idx.minus])
                          breaks.minus <- GRanges(seqnames = chr, ranges = IRanges(start = s, end = e), strand = "-")
                          # Truncate transcritps by transcription start position inside the peak (estimate by predictStrand function)
                          trx.minus <- setdiff(trx.minus, breaks.minus)
                  }

                  # Combine transcripts
                  trx <- GenomicRanges::sort.GenomicRanges(c(trx.plus, trx.minus), decreasing = FALSE, ignore.strand = TRUE)
                  trx$id <- paste0("trx_", seq_along(trx))
                  trx$length <- width(trx)

                  message("Done!")

                  if (msg.alarm.plus){ message("[INFO] No peaks to break transcripts on the plus strand.") }
                  if (msg.alarm.minus){ message("[INFO] No peaks to break transcripts on the minus strand.") }

                  if (estimate.params){
                          # Calculate transcripts coverage (proportion of bases covered by fragments/reads) and fpkm
                          message("[INFO] Calculating fpkm and coverage ...", appendLF = FALSE)

                          chrs <- levels(droplevels(seqnames(trx)))
                          trx <- do.call(c, lapply(chrs,
                                                   function(chr){
                                                           trx.plus <- trx[seqnames(trx) == chr & strand(trx) == "+"]
                                                           trx.minus <- trx[seqnames(trx) == chr & strand(trx) == "-"]

                                                           # covearge profile per transcript
                                                           views.plus <- Views(subject = tdsObj@coveragePlus[[chr]],
                                                                               start = start(trx.plus),
                                                                               end = end(trx.plus))

                                                           views.minus <- Views(subject = tdsObj@coverageMinus[[chr]],
                                                                                start = start(trx.minus),
                                                                                end = end(trx.minus))

                                                           # amount of bases covered by fragments
                                                           trx.plus$bases.covered <- viewApply(views.plus,
                                                                                               function(x){
                                                                                                       sum(x > tdsObj@transcriptsCovCutoff)
                                                                                               }
                                                           )
                                                           trx.minus$bases.covered <- viewApply(views.minus,
                                                                                                function(x){
                                                                                                        sum(x > tdsObj@transcriptsCovCutoff)
                                                                                                }
                                                           )
                                                           trx.chr <- c(trx.plus, trx.minus)
                                                           trx.chr$coverage <- round(trx.chr$bases.covered / width(trx.chr), 3)
                                                           trx.chr
                                                   }
                          ))

                          # Calculate fpkm
                          counts <- countOverlaps(query = trx,
                                                  subject = tdsObj@fragments,
                                                  ignore.strand = FALSE)
                          trx$fragments <- counts
                          trx$fpkm <- round(.rpkm(C = counts, S = width(trx), N = tdsObj@transcriptsNormalization), 3)
                          trx <- GenomicRanges::sort(x = trx, decreasing = FALSE, ignore.strand = TRUE)

                          message("Done!")
                  }

                  tdsObj@transcripts <- trx
                  assign(objName, tdsObj, envir = parent.frame())

          }
)


####------------------- annotateTranscripts -------------------####
#' @rdname annotateTranscripts-methods
setMethod("annotateTranscripts",
          signature = signature("TranscriptionDataSet", "GRanges"),
          definition = function(object, annot, min.overlap = 0.3){

                  if (length(object@transcripts) == 0){
                          stop("[ERROR] 'object' doesn't contain transcripts information. Please, call detectTranscripts function on the object of class 'TranscriptionDataSet'.")
                  }

                  if (is.null(names(annot))){
                          stop("[ERROR] Names are not specified in 'annot'.")
                  }

                  if (is.numeric(min.overlap)){
                          if (min.overlap > 1 | min.overlap <= 0){
                                  stop("[ERROR] 'min.overlap' needs to be in a range [0,1].")
                          }
                  } else {
                          stop("[ERROR] 'min.overlap' needs to be numeric.")
                  }

                  objName <- deparse(substitute(object))
                  trx <- object@transcripts

                  # Look for the overlap betwen transcripts and annotations
                  hits <- findOverlaps(trx, annot)

                  # Determine the width of the overlap (picewise intersection)
                  inter <- pintersect(annot[subjectHits(hits)], trx[queryHits(hits)])

                  # Calculate proportion of the overlap
                  # proportion of transcript overlaping annotation
                  trx.overlap.prop <- width(inter) / width(trx[queryHits(hits)])
                  # proportion of annotation overlaping transcript
                  annot.overlap.prop <- width(inter) / width(annot[subjectHits(hits)])
                  # take average of two proportions
                  overlap.prop <- (trx.overlap.prop + annot.overlap.prop) / 2
                  idx <- which(overlap.prop >= min.overlap)

                  # Select overlaps with the overlap proportion larger than 'min.overlap'
                  hits <- hits[idx]

                  # Split annot by transcripts
                  annot.name.by.trx <- split(names(inter[idx]), queryHits(hits))
                  annot.name.list <- CharacterList(replicate(n = length(trx), list("ND")))
                  annot.name.list[as.numeric(names(annot.name.by.trx))] <- annot.name.by.trx
                  trx$annotation.overlap <- annot.name.list

                  object@transcripts <- trx
                  assign(objName, object, envir = parent.frame())
          }
)


####------------------- getTranscripts -------------------####
#' @rdname getTranscripts-methods
setMethod("getTranscripts",
          signature = signature("TranscriptionDataSet"),
          definition = function(object,
                                min.length,
                                min.fpkm,
                                min.coverage){

                  trx <- object@transcripts

                  # Tests
                  # Filter transcripts by length
                  if (!missing(min.length)){
                          if (is.numeric(min.length)){
                                  if (min.length < 0){
                                          stop("[ERROR] 'min.length' needs to be a positive number.")
                                  }
                          } else {
                                  stop("[ERROR] 'min.length' needs to be numeric.")
                          }

                          trx <- trx[width(trx) > min.length]
                  }


                  # Filter transcripts by fpkm
                  if (!missing(min.fpkm)){

                          if (is.numeric(min.fpkm)){
                                  if (min.fpkm < 0){
                                          stop("[ERROR] 'min.fpkm' needs to be a positive number.")
                                  }
                          } else {
                                  stop("[ERROR] 'min.fpkm' needs to be numeric.")
                          }

                          if (!("fpkm" %in% names(values(trx)))){
                                  stop("[ERROR] fpkm values are not available. Please, rerun detectTranscripts function with 'estimate.params'=TRUE.")
                          }

                          trx <- trx[trx$fpkm > min.fpkm]

                  }

                  # Filter transcripts by coverage (proportion of bases covered by fragments)
                  if (!missing(min.coverage)){
                          if (is.numeric(min.coverage)){
                                  if (min.coverage < 0 | min.coverage > 1){
                                          stop("[ERROR] 'min.coverage' needs to be a range [0, 1].")
                                  }
                          } else {
                                  stop("[ERROR] 'min.coverage' needs to be numeric.")
                          }

                          if (!("coverage" %in% names(values(trx)))){
                                  stop("[ERROR] coverage information is not available. Please, rerun detectTranscripts function with 'estimate.params'=TRUE.")
                          }

                          trx <- trx[trx$coverage > min.coverage]

                  }

                  return(trx)
          }
)


####------------------- getTestedGapDistances -------------------####
#' @rdname getTestedGapDistances-methods
setMethod("getTestedGapDistances",
          signature = signature("TranscriptionDataSet"),
          definition = function(object){
                          return(object@gapDistanceTest)
          }
)


####------------------- plotErrorRate -------------------####
#' @rdname plotErrorRate-methods
setMethod("plotErrorRate",
          signature = signature("TranscriptionDataSet"),
          definition = function(object,
                                color = c("#1B9E77", "#D95F02", "#7570B3"),
                                xlab = "Gap distance (kb)",
                                ylab = "Error rate",
                                ...){

                  if (length(object@gapDistanceTest) == 0){
                          stop("[ERROR] 'gapDistanceTest' slot is empty.
                               Make sure that the 'object' was initialized with the constructTDS function call.
                               Please, call the estimateGapDistance function to test multiple gap distances.")
                  }

                  if (length(color) != 3){
                          stop("[ERROR] 'color' needs to contain three colors.")
                  }


                  df <- object@gapDistanceTest

                  graphics::plot(df$sum.two.errors, type = "l", col = color[1], ylim = c(0, 1), xaxt = "n",
                                 ylab = ylab, xlab = xlab, ...)
                  thicks <- which(df$gap.distance %% 1000 == 0)
                  labels <- df$gap.distance[thicks]/1000
                  graphics::axis(side = 1, at = thicks, labels = labels)
                  graphics::lines(df$error.dissected, col = color[2], ...)
                  graphics::lines(df$error.merged, col = color[3], ...)
                  graphics::abline(v = which.min(df$sum.two.errors), col = "grey40", lwd = 1.5, lty = "dashed")
                  graphics::legend("topright",
                                   legend = c("Sum of two errors", "Dissected error", "Merged error"),
                                   col = color, lwd = 2, cex = 0.8)
          }
)


####------------------- exportCoverage -------------------####
#' @rdname exportCoverage-methods
setMethod("exportCoverage",
          signature = signature("TranscriptionDataSet"),
          definition = function(object,
                                file,
                                type,
                                strand,
                                color,
                                filter.by.coverage.cutoff = FALSE,
                                coverage.cutoff = NULL,
                                rpm = FALSE,
                                total.reads){

                  if (length(object@coveragePlus) == 0 | length(object@coverageMinus) == 0){
                          stop("[ERROR] 'object' doesn't contain coverage information.
                               Please, use constructTDS to initialize an object of class 'TranscriptionDataSet'.")
                  }

                  if (missing(type)){
                          type <- "bedGraph"
                  } else {
                          type <- match.arg(arg = type, choices = c("bigWig", "bedGraph"))
                  }

                  if (missing(color)){
                          color <- c(0L, 0L, 255L)
                  } else {
                          if (is.numeric(color)){
                                  color <- as.integer(color)
                          } else {
                                  stop("[ERROR] 'color' needs to be an object o class
                                       'integer' representing the track color (as from col2rgb).")
                          }
                  }

                  if (strand %in% "+"){
                          cov <- object@coveragePlus
                  } else if (strand %in% "-"){
                          cov <- object@coverageMinus * (-1)
                  } else {
                          stop("[ERROR] Wrong 'strand'. Either '+' or '-'.")
                  }

                  if (!is.logical(rpm)){
                          stop("[ERROR] 'rpm' needs to be boolean.")
                  }

                  if (!missing(total.reads)){
                          if (!is.numeric(total.reads)){
                                  stop("[ERROR] 'total.reads' needs to be numeric.")
                          }
                  } else {
                          total.reads <- length(object@fragments)
                  }

                  # Filter by coverage cutoff
                  if (is.logical(filter.by.coverage.cutoff)){
                          if (filter.by.coverage.cutoff){
                                  if (is.null(coverage.cutoff)){
                                          if (length(object@coverageCutoff) == 1){
                                                  coverage.cutoff <- object@coverageCutoff
                                          } else {
                                                  stop("[ERROR] 'coverageCutoff' slot is empty.
                                                       Please, specify 'coverage.cutoff' manually or call estimateBackground function on the 'object'.")
                                          }
                                  } else {
                                          if (is.numeric(coverage.cutoff)){
                                                  if (coverage.cutoff <= 0) {
                                                          stop("[ERROR] 'coverage.cutoff' needs to be above 0.")
                                                  }
                                          } else {
                                                  stop("[ERROR] 'coverage.cutoff' needs to be numeric.")
                                          }
                                  }
                                  cov[abs(cov) < coverage.cutoff] <- 0
                          } else {
                                  if (!is.null(coverage.cutoff)){
                                          message("[WARNING] 'filter.by.coverage.cutoff' is set to FALSE. Coverage profiles will not be filtered by 'coverage.cutoff'.")
                                  }
                          }
                  } else{
                          stop("[ERROR] 'filter.by.coverage.cutoff' needs to be boolean.")
                  }

                  # Normalize
                  if (rpm){
                          cov.rpm <- (10^6 * cov) / total.reads
                          names(cov.rpm) <- names(cov)
                          cov <- cov.rpm
                  }

                  if (type %in% "bigWig"){
                          rtracklayer::export(object = cov, con = file, format = "bigWig")
                  } else if (type %in% "bedGraph"){
                          track.line <- new("GraphTrackLine",
                                            name = file,
                                            color = color,
                                            altColor = color,
                                            autoScale = TRUE,
                                            alwaysZero = TRUE,
                                            type="bedGraph")

                          rtracklayer::export(object = cov, con = file, trackLine = track.line, format = "bedGraph")
                  }
          }
)


####------------------- transcriptsToBed -------------------####
#' @rdname transcriptsToBed-methods
setMethod("transcriptsToBed",
          signature = signature("GRanges"),
          definition = function(object, file, strand.color = c("blue", "red")){

                  if (length(strand.color) != 2){
                          stop("[ERROR] 'strand.color' needs to be a character vector of length two.")
                  }

                  chrom <- seqnames(object)
                  chromStart <- start(object)
                  chromEnd <- end(object)
                  # name <- paste0("transcript_", 1 : length(object))
                  ids <- object$id
                  score <- rep(0, length(object))
                  strand <- strand(object)
                  thichStart <- chromStart
                  thickEnd <- chromEnd
                  col.plus <- paste(grDevices::col2rgb(col = strand.color[1]), collapse = ",")
                  col.minus <- paste(grDevices::col2rgb(col = strand.color[2]), collapse = ",")
                  itemRgb <- ifelse(strand(object) == "+", col.plus, col.minus)

                  df <- data.frame(chrom, chromStart, chromEnd, ids, score, strand, thichStart, thickEnd, itemRgb)
                  utils::write.table(x = df, file = file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
          }
)

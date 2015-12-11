####------------------- show -------------------####
#' Show method for ChipDataSet and TranscriptionDataSet objects
#'
#' Show method for objects of class \code{\link{ChipDataSet}} and \code{\link{TranscriptionDataSet}}
#'
#' @name show
#' @rdname show
#' @aliases show show,TranscriptionDataSet-method
#' @param object A \code{\link{TranscriptionDataSet}} or
#'     \code{\link{ChipDataSet}} object.
#' @return
#'     Displays an overview of the TranscriptionDataSet or ChipDataSet object.
#' @examples
#'     ### Load TranscriptionDataSet object
#'     data(tds)
#'
#'     ### View a short summary of the object
#'     tds
#' @export
setMethod("show",
          signature = signature("TranscriptionDataSet"),
          definition = function(object){
                  cat("S4 object of class", class(object), "\n")

                  cat("=======================================\n")

                  cat("fragments:", length(object@fragments), "\n")

                  cat("fragmentSize:", object@fragmentSize, "\n")

                  cat("coveragePlus:", names(object@coveragePlus), "\n")

                  cat("coverageMinus:", names(object@coverageMinus), "\n")

                  if (length(object@coverageCutoffFdr) == 0){
                          cat("coverageCutoff(): ", object@coverageCutoff, "\n", sep = "")
                  } else {
                          cat("coverageCutoff(Params: fdr=", object@coverageCutoffFdr, "): ",
                              object@coverageCutoff, "\n", sep = "")
                  }

                  if (length(object@gapDistanceTest) > 0){
                          idx <- which.min(object@gapDistanceTest$sum.two.errors)
                          gap.dist <- object@gapDistanceTest$gap.distance[idx]
                          sum.two.errors <- object@gapDistanceTest$sum.two.errors[idx]
                  } else {
                          gap.dist <- numeric()
                          sum.two.errors <- numeric()
                  }

                  if (length(object@gapDistanceTestCovCutoff) == 0){
                          cat("gapDistance(): ", gap.dist, "\n")

                  } else{
                          cat("gapDistance(Params: coverage.cutoff=", object@gapDistanceTestCovCutoff, "): ",
                              gap.dist, " [error rate=", round(sum.two.errors, 3), "]", "\n", sep = "")

                  }

                  if (length(object@transcriptsCovCutoff) == 0 | length(object@transcriptsGapDist) == 0){
                          cat("transcripts():", length(object@transcripts), "\n")

                  } else{
                          cat("transcripts(Params: coverage.cutoff=", object@transcriptsCovCutoff, "; ",
                              "gap.dist=", object@transcriptsGapDist, "): ", sep = "",
                              length(object@transcripts), " transcripts detected", "\n")
                  }
          }
)

#' @rdname show
#' @aliases show show,ChipDataSet-method
setMethod("show",
          signature = signature("ChipDataSet"),
          definition = function(object){
                  cat("S4 object of class", class(object), "\n")

                  cat("=======================================\n")

                  cat("Peaks:", length(object@peaks), "\n")

                  cat("Chromosome distribution:", as.character(S4Vectors::runValue(seqnames(object@peaks))), "\n")

                  if (length(object@features) > 0){
                          features <- names(object@features)
                  } else {
                          features <- character()
                  }
                  cat("Features: ", features, "\n")

                  if (!is.null(object@tssOverlapPrediction$predicted.tssOverlap)){
                          peaks.gene.id <- which(object@tssOverlapPrediction$predicted.tssOverlap$predicted.tssOverlap == "yes")
                          gene.overlap <- length(peaks.gene.id)
                  } else {
                          gene.overlap <- numeric()
                  }
                  cat("Gene associated peaks (predicted): ", gene.overlap, "\n")

                  if (!is.null(object@strandPrediction$predicted.strand)){
                          peak.strand <- table(object@strandPrediction$predicted.strand)
                          peak.strand <- paste0(peak.strand, "(", names(peak.strand), ")")
                          peaks.strand.gene <- table(object@strandPrediction$predicted.strand[peaks.gene.id])
                          peaks.strand.gene <- paste0(peaks.strand.gene, "(", names(peaks.strand.gene), ")")
                  } else {
                          peak.strand <- character()
                          peaks.strand.gene <- character()
                  }
                  # cat("Peaks (all) by strand: ", peak.strand, "\n")
                  cat("Gene associated peaks by strand: ", peaks.strand.gene, "\n")
          }
)

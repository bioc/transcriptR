#' TranscriptionDataSet
#'
#' The \code{TranscriptionDataSet} is a S4 class to store input values,
#' intermediate calculations and results of the transcripts detection and
#' quantification analysis.
#'
#' @name TranscriptionDataSet-class
#' @slot bamFile \code{Character}. Path to a BAM file.
#' @slot fragments \code{\link[GenomicRanges]{GRanges}}. Sequencing reads
#'     extended to the fragment size.
#' @slot fragmentSize \code{Numeric}. Fragment size in base pairs (bp).
#' @slot region \code{\link[GenomicRanges]{GRanges}}. Genomic region(s) to
#'     extract reads from.
#' @slot coveragePlus \code{\link[IRanges]{RleList}}. Fragment coverage profile
#'     for the forward DNA strand.
#' @slot coverageMinus \code{\link[IRanges]{RleList}}. Fragment coverage profile
#'     for the reverse DNA strand.
#' @slot coverageCutoff \code{Numeric}. Background coverage cutoff value.
#' @slot coverageCutoffFdr \code{Numeric}. False Discovery Rate (FDR) used to
#'     estimate background coverage cutoff.
#' @slot gapDistanceTest \code{Data.Frame}. Tested gap distances and
#'     corresponding error rates.
#' @slot gapDistanceTestCovCutoff \code{Numeric}. Coverage cutoff value used for
#'     the gap distance estimation.
#' @slot transcripts \code{\link[GenomicRanges]{GRanges}}. Identified transcripts.
#' @slot transcriptsCovCutoff \code{Numeric}. Coverage cutoff value used for the
#'     transcripts detection.
#' @slot transcriptsGapDist \code{Numeric}. Gap distance value used for the
#'     transcripts detection.
#' @slot transcriptsNormalization \code{Numeric}. Total number of reads used for
#'     normalization when calculating FPKM.
#'
#' @seealso \code{\link{constructTDS}}
#'
#' @author Armen R. Karapetyan
#'
#' @import GenomicRanges
#' @import IRanges
#' @import methods
#' @importFrom BiocGenerics Reduce
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom GenomeInfoDb seqlevelsInUse
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicFeatures exons
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom GenomicFeatures promoters
#' @importFrom GenomicFeatures transcripts
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors runValue
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom caret confusionMatrix
#' @importFrom caret createDataPartition
#' @importFrom caret train
#' @importFrom caret trainControl
#' @importFrom caret twoClassSummary
#' @importFrom chipseq peakCutoff
#' @importFrom grDevices col2rgb
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom pROC plot.roc
#' @importFrom pROC roc
#' @importFrom reshape2 melt
#' @importFrom rtracklayer export
#' @importFrom stats coef
#' @importFrom stats complete.cases
#' @importFrom stats predict
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @import e1071
#'
#' @exportClass TranscriptionDataSet

TranscriptionDataSet <- setClass(Class = "TranscriptionDataSet",

                                 slots = c(bamFile = "character",
                                           fragments = "GRanges",
                                           fragmentSize = "numeric",
                                           region = "GRanges",
                                           coveragePlus = "RleList",
                                           coverageMinus = "RleList",
                                           coverageCutoff = "numeric",
                                           coverageCutoffFdr = "numeric",
                                           gapDistanceTest = "data.frame",
                                           gapDistanceTestCovCutoff = "numeric",
                                           transcripts = "GRanges",
                                           transcriptsCovCutoff = "numeric",
                                           transcriptsGapDist = "numeric",
                                           transcriptsNormalization = "numeric"),

                                 prototype = c(bamFile = character(),
                                               fragments = GenomicRanges::GRanges(),
                                               fragmentSize = numeric(),
                                               region = GenomicRanges::GRanges(),
                                               coveragePlus = IRanges::RleList(),
                                               coverageMinus = IRanges::RleList(),
                                               coverageCutoff = numeric(),
                                               coverageCutoffFdr = numeric(),
                                               gapDistanceTest = data.frame(),
                                               gapDistanceTestCovCutoff = numeric(),
                                               transcripts = GenomicRanges::GRanges(),
                                               transcriptsCovCutoff = numeric(),
                                               transcriptsGapDist = numeric(),
                                               transcriptsNormalization = numeric())
)


# documents data sets that come with the package

#' Example of \code{TranscriptionDataSet} object.
#'
#' \code{tds} is an object of \code{\link{TranscriptionDataSet}} class, containing
#' nuclear RNA-seq data for human chromosome 15 (chr15:63261757-84081194),
#' profiled in prostate cancer LNCaP cells.
#'
#' @name tds
#' @docType data
#' @format \code{\link{TranscriptionDataSet}} object
#' @return \code{\link{TranscriptionDataSet}} object
"tds"

#' Reference annotation (knownGene from UCSC)
#'
#' \code{annot} is an object of \code{\link{GRanges}} class, containing
#' genomic coordinates of the genes located on human chromosome 15 (chr15:63261757-84081194).
#'
#' @name annot
#' @docType data
#' @format \code{\link{GRanges}} object
#' @return \code{\link{GRanges}} object
"annot"

#' ChipDataSet
#'
#' The \code{ChipDataSet} is a S4 class to store input values, intermediate
#' calculations and results of ChIP-seq peaks analysis.
#'
#' @name ChipDataSet-class
#' @slot peaks \code{\link[GenomicRanges]{GRanges}}. ChIP-seq peaks.
#' @slot region \code{\link[GenomicRanges]{GRanges}}. Genomic region(s) to extract peaks from.
#' @slot genomicAnnotation \code{Data.Frame}. Genomic distribution of the peaks
#'     at distinct genomic features (TSSs, exons, introns, intergenic regions).
#' @slot features \code{Data.Frame}. Estimated characteristics (features) of the
#'     peaks.
#' @slot tssOverlapPrediction \code{List}. Prediction of the gene associated
#'     peaks. The following elements are reported:
#'     \itemize{
#'         \item 'predicted.tssOverlap' - predicted class (yes - gene associated;
#'             no - background) and probability of a ChIP-seq peak being
#'             classified as gene associated.
#'         \item 'confusionMatrix' - cross-tabulation of observed and predicted
#'             classes with associated statistics.
#'         \item 'logFitSummary' - statistical significance of the predictors
#'             used in the logistic regression.
#'         \item 'roc' - results of the receiver operating characteristic
#'             analysis.
#'         }
#' @slot strandPrediction \code{List}. Prediction of the peak strandedness. The
#'     following elements are reported:
#'     \itemize{
#'         \item 'predicted.strand' - predicted ChIP-seq peak strand.
#'         \item 'probability.cutoff' - probability cutoff for q2.
#'         \item 'results.plus' - intermediate calculations for the forward DNA
#'             strand.
#'         \item 'results.minus' - intermediate calculations for the reverse DNA
#'             strand.
#'         }
#'
#' @seealso \code{\link{constructCDS}}
#'     \code{\link{predictTssOverlap}}
#'     \code{\link{predictStrand}}
#'
#' @author Armen R. Karapetyan
#'
#' @exportClass ChipDataSet
ChipDataSet <- setClass(Class = "ChipDataSet",

                        slots = c(peaks = "GRanges",
                                  region = "GRanges",
                                  genomicAnnotation = "data.frame",
                                  features = "data.frame",
                                  tssOverlapPrediction = "list",
                                  strandPrediction = "list"),

                        prototype = list(peaks = GenomicRanges::GRanges(),
                                         region = GenomicRanges::GRanges(),
                                         genomicAnnotation = data.frame(),
                                         features = data.frame(),
                                         tssOverlapPrediction = list(predicted.tssOverlap = NULL,
                                                                     confusionMatrix = NULL,
                                                                     logFitSummary = NULL,
                                                                     roc = NULL),
                                         strandPrediction = list(predicted.strand = NULL,
                                                                 probability.cutoff = NULL,
                                                                 results.plus = NULL,
                                                                 results.minus = NULL))
)


# documents data sets that come with the package

#' Example of \code{ChipDataSet} object.
#'
#' \code{cds} is an object of \code{\link{ChipDataSet}} class, containing
#' H3K4me3 active histone mark ChIP-seq peaks from human chromosome 15
#' (chr15:63261757-84081194) profiled in prostate cancer LNCaP cells.
#'
#' @name cds
#' @docType data
#' @format  \code{\link{ChipDataSet}} object
#' @return \code{\link{ChipDataSet}} object
"cds"

####------------------- predictTssOverlap -------------------####
#' predictTssOverlap
#'
#' The function classifies ChIP-seq peaks on gene associated and background using
#' classification model based on a logistic regression.
#'
#' @name predictTssOverlap
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#' @param feature \code{Character}. Specify feature(s) that are be used as
#'     predictors in the classification model. By default, all the features
#'     stored in the \code{\link{ChipDataSet}} object will be used.
#' @param p \code{Numeric}. The proportion of data that is used to train the
#'     model. Default: 0.8.
#'
#' @details
#'     In order to improve the accuracy of the model the data is internally
#'     partitioned into a training and testing data sets. A repeated 10-Fold
#'     Cross-Validation is used to calculate performance measures on the
#'     training data set and to prevent over-fitting.
#'
#'     The model fit and validation is internally acomplised by the functions
#'     implemented in the \href{http://topepo.github.io/caret/index.html}{caret}
#'     package.
#'
#' @return
#'     The slot \code{tssOverlapPrediction} of the provided
#'     \code{TranscriptionDataSet} object will be updated by the the following
#'     elements: 'predicted.tssOverlap', 'confusionMatrix', 'logFitSummary' and
#'     'roc'.
#'
#' @seealso
#'     \code{\link{ChipDataSet}}
#'     \code{\link{constructCDS}}
#'     \href{http://topepo.github.io/caret/index.html}{caret}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' predictTssOverlap(object = cds, feature = "pileup", p = 0.75)
#'
#' ### View a short summary of the gene associated peaks prediction
#' cds
#'
#' ### View peaks and associated prediction
#' getPeaks(cds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname predictTssOverlap-methods
setGeneric(name = "predictTssOverlap",
           def = function(object, feature, p = 0.8)
                   standardGeneric("predictTssOverlap")
)


####------------------- predictStrand -------------------####
#' predictStrand
#'
#' The function evaluates transcription initiation within a peak region by
#' comparing RNA-seq read densities upstream and downstream of an empirically
#' determined transcription start sites. Putative transcription of both forward
#' and reverse genomic strands is tested and the results are stored with each
#' ChIP-seq peak.
#'
#' @name predictStrand
#' @docType methods
#'
#' @param cdsObj A \code{\link{ChipDataSet}} object.
#' @param tdsObj A \code{\link{TranscriptionDataSet}} object.
#' @param coverage.cutoff \code{Numeric}. A cutoff value to discard regions with
#'     the low fragments coverage, representing expression noise. By default,
#'     the value stored in the \code{coverageCutoff} slot of the supplied
#'     \code{TranscriptionDataSet} object is used. The optimal cutoff value can
#'     be calculated by \code{\link{estimateBackground}} function call.
#' @param quant.cutoff \code{Numeric}. A cutoff value for the cumulative
#'     distribution of the RNA-seq signal along the ChIP-seq peak region. Must
#'     be in a range (0, 1). For the details, see step 1 in the "Details"
#'     section below. Default: 0.1.
#' @param win.size \code{Numeric}. The size of the q1 and q2 regions
#'     flanking transcription start position at the 5' and 3', respectively.
#'     For the details, see step 2 in the "Details"  section below.
#'     Default: 2500.
#' @param prob.cutoff \code{Numeric}. A cutoff value for the probability of reads
#'     to be sampled from the q2 flanking region. If not supplied, the value
#'     estimated from the data will be used. Must be in a range (0, 1). For the
#'     details, see step 6 in the "Details" section below.
#'
#' @details
#'     RNA-seq data is incorporated to find direct evidence of active
#'     transcription from every putatively gene associated peak. In order to do
#'     this, we determine the 'strandedness' of the ChIP-seq peaks, using strand
#'      specific RNA-seq data. The following assumptions are made in order to
#'      retrieve the peak 'strandedness':
#'     \itemize{
#'         \item The putatively gene associated ChIP-seq peaks are commonly
#'             associated with transcription initiation.
#'
#'         \item This transcription initiation occurs within the ChIP peak region.
#'
#'         \item When a ChIP peak is associated with a transcription initiation
#'             event, we expect to see a strand-specific increase in RNA-seq
#'             fragment count downstream the transcription initiation site.
#'     }
#'
#'    Each peak in the data set is tested for association with transcription
#'    initiation on both strands of DNA. Steps 1-5 are performed for both
#'    forward and reverse DNA strand separately and step 6 combines the data
#'    from both strands. If the peak is identified as associated with the
#'    transcription on both strands, than it is considered to be a bidirectional.
#'
#'    ChIP peak 'strandedness' prediction steps:
#'     \enumerate{
#'         \item Identify a location within the ChIP-seq peak near the
#'         transcription start site. This is accomplished by calculating the
#'         cumulative distribution of RNA-seq fragments within a peak region.
#'         The position is determined where 100\% - 'quant.cutoff' * 100\% of
#'         RNA-seq fragments are located downstream. This approach performs well
#'         on both gene-poor and gene-dense regions where transcripts may overlap.
#'
#'         \item Two equally sized regions are defined (q1 and q2), flanking the
#'         position identified in (1) on both sides. RNA-seq fragments are
#'         counted in each region.
#'
#'         \item ChIP peaks with an RNA-seq fragment coverage below an estimated
#'         threshold are discarded from the analysis.
#'
#'         \item The probability is calculated for RNA-seq fragments to be
#'         sampled from either q1 or q2. Based on the assumptions we stated
#'         above, a ChIP peak that is associated with transcription initiation
#'         should have more reads in q2 (downstream of the transcription start
#'         position) compared to q1, and subsequently, the probability of a
#'         fragment being sampled from q2 would be higher.
#'
#'         \item ChIP-seq peaks are divided into gene associated and background
#'         based on the prediction.
#'
#'         \item Iteratively, the optimal P(q2) threshold is identified, which
#'         balances out the False Discovery Rate (FDR) and False Negative Rate
#'         (FNR). Peaks with the P(q2) exceeding the estimated threshold are
#'         considered to be associated with the transcription initiation event.
#'     }
#'
#' @return The slot \code{strandPrediction} of the provided
#'     \code{\link{ChipDataSet}} object will be updated by the the following
#'     elements: 'predicted.strand', 'probability.cutoff', 'results.plus' and
#'     'results.minus'.
#'
#' @seealso \code{\link{ChipDataSet}} \code{\link{constructCDS}}
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### Classify peaks on gene associated and background
#' predictTssOverlap(object = cds, feature = "pileup", p = 0.75)
#'
#' ### Predict peak 'strand'
#' predictStrand(cdsObj = cds, tdsObj = tds, coverage.cutoff = 5,
#' quant.cutoff = 0.1, win.size = 2500)
#'
#' ### View a short summary of the 'strand' prediction
#' cds
#'
#' ### View 'strand' prediction
#' getPeaks(cds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname predictStrand-methods
setGeneric(name = "predictStrand",
           def = function(cdsObj, tdsObj, coverage.cutoff,
                          quant.cutoff = 0.1, win.size = 2500, prob.cutoff)
                   standardGeneric("predictStrand")
)


####------------------- addFeature -------------------####
#' addFeature
#'
#' Add feature(s) to the \code{\link{ChipDataSet}} object.
#'
#' @name addFeature
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#' @param feature \code{List}. User specified characteristics of the peaks.
#'
#' @return
#'     The user specified feature(s) will be added to the slot \code{features}
#'     of the provided  \code{\link{ChipDataSet}} object.
#'
#' @seealso \code{\link{constructCDS}} \code{\link{predictTssOverlap}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### As an example create some fake data
#' N <- length(getPeaks(cds))
#' fake.data <- rnorm(n = N)
#'
#' addFeature(object = cds, feature = list(fake = fake.data))
#'
#' ### View newly added feature
#' getPeaks(cds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname addFeature-methods
setGeneric(name = "addFeature",
           def = function(object, feature)
                   standardGeneric("addFeature")
)


####------------------- getPeaks -------------------####
#' getPeaks
#'
#' Retrieve ChIP-seq peak information from the \code{\link{ChipDataSet}} object.
#'
#' @name getPeaks
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @seealso
#'     \code{\link{constructCDS}}
#'     \code{\link{predictTssOverlap}}
#'     \code{\link{predictStrand}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' getPeaks(cds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname getPeaks-methods
setGeneric(name = "getPeaks",
           def = function(object)
                   standardGeneric("getPeaks")
)


####------------------- getGenomicAnnot -------------------####
#' getGenomicAnnot
#'
#' Retrieve genomic distribution of ChIP-seq peaks at distinct genomic features
#' (exons, introns, TSSs, intergenic regions)
#'
#' @name getGenomicAnnot
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#'
#' @return
#'     A four column \code{Data.Frame}, storing information about observed and
#'     expected number of peaks at distinct genomic features.
#'
#' @details
#'     A simple quality check of the supplied ChIP-seq peaks can be performed by
#'     investigating their genomic distribution. Ideally, these peaks should
#'     demonstrate substantial enrichment at TSS regions. Enrichment of the
#'     peaks at a given genomic feature (e.g. TSS) is defined as the ratio
#'     between the observed and expected number of peaks. The expected number of
#'     peaks is calculated from the proportion of the genome covered by the
#'     given genomic feature.
#'
#' @seealso \code{\link{constructCDS}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' getGenomicAnnot(cds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname getGenomicAnnot-methods
setGeneric(name = "getGenomicAnnot",
           def = function(object)
                   standardGeneric("getGenomicAnnot")
)


####------------------- getConfusionMatrix -------------------####
#' getConfusionMatrix
#'
#' Retrieve a cross-tabulation of observed and predicted classes (prediction of
#' gene associated peaks) with associated statistics.
#'
#' @name getConfusionMatrix
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#'
#' @return
#'     An object of \code{Confusion Matrix} class. For the details see
#'     \href{http://topepo.github.io/caret/index.html}{caret} package.
#'
#' @seealso
#'     \code{\link{predictTssOverlap}} \code{\link[caret]{confusionMatrix}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' getConfusionMatrix(cds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname getConfusionMatrix-methods
setGeneric(name = "getConfusionMatrix",
           def = function(object)
                   standardGeneric("getConfusionMatrix")
)


####------------------- getPredictorSignificance -------------------####
#' getPredictorSignificance
#'
#' Retrieve significance of each predictor used in the classification model fit
#' (prediction of gene associated peaks).
#'
#' @name getPredictorSignificance
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#'
#' @return A vector of p-values.
#'
#' @seealso \code{\link{predictTssOverlap}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' getPredictorSignificance(cds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname getPredictorSignificance-methods
setGeneric(name = "getPredictorSignificance",
           def = function(object)
                   standardGeneric("getPredictorSignificance")
)


####------------------- getQuadProb -------------------####
#' getQuadProb
#'
#' Retrieve all internal calculations performed by \code{\link{predictStrand}}
#' function.
#'
#' @name getQuadProb
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#' @param strand \code{Character}. The strand to extract calculations for.
#'     One of ["+", "-"].
#'
#' @return
#'     A nine column \code{Data.Frame}, where each row corresponds to
#'     a ChIP-seq peak and each column keeps one of the intermediate calculations:
#'     \itemize{
#'         \item \code{max.cov} - maximum coverage of the RNA-seq fragments
#'             inside the peak region.
#'
#'         \item \code{pass.cov.treshold} - whether the \code{max.cov} exceeds
#'             the \code{coverage.cutoff}, either user defined or estimated from
#'             RNA-seq data by \code{\link{estimateBackground}} function call and
#'             stored in \code{\link{TranscriptionDataSet}} object.
#'
#'         \item \code{q1q2.sepline.coord} - genomic coordinate corresponding to
#'             the transcription start position inside the peak region.
#'
#'         \item \code{q1.coord} - genomic coordinates of q1.
#'
#'         \item \code{q2.coord} - genomic coordinates of q2.
#'
#'         \item \code{q1.count} - total number of fragments in q1.
#'
#'         \item \code{q2.count} - total number of fragments in q2.
#'
#'         \item \code{q1.prob} - probability of a fragment being sampled from
#'             the q1.
#'
#'         \item \code{q2.prob} - probability of a fragment being sampled from
#'            the q2.
#'     }
#'
#' @seealso \code{\link{predictStrand}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' head(getQuadProb(cds, strand = "+"))
#' head(getQuadProb(cds, strand = "-"))
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname getQuadProb-methods
setGeneric(name = "getQuadProb",
           def = function(object, strand)
                   standardGeneric("getQuadProb")
)


####------------------- getProbTreshold -------------------####
#' getProbTreshold
#'
#' Retrieve estimated P(q2) threshold, used to select peaks with a putative
#' transcription initiation event.
#'
#' @name getProbTreshold
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#'
#' @return Estimated P(q2) threshold.
#'
#' @seealso \code{\link{predictStrand}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' getProbTreshold(cds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname getProbTreshold-methods
setGeneric(name = "getProbTreshold",
           def = function(object)
                   standardGeneric("getProbTreshold")
)


####------------------- plotGenomicAnnot -------------------####
#' plotGenomicAnnot
#'
#' Visualize genomic distribution of ChIP-seq peaks.
#'
#' @name plotGenomicAnnot
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#' @param plot.type \code{Character}. One of ("distrib", "enrich").
#'     Default: "enrich".
#' @param xlab \code{Character}. Title of the x-axis.
#' @param ylab \code{Character}. Title of the y-axis.
#' @param color A character vector of length four, specifying colors for
#'     distinct genomic features (TSSs, exons, introns, intergenic regions).
#'
#' @details
#'     Genomic distribution of the peaks can be visualized in two ways, either
#'     by observing the total number of peaks overlapping given genomic feature
#'     or by looking at the enrichment levels.
#'
#' @return \code{\link{ggplot2}} object.
#'
#' @seealso \code{\link{constructCDS}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### Plot the total number of peaks overlapping distinct genomic features
#' plotGenomicAnnot(object = cds, plot.type = "distrib")
#'
#' ### Plot enrichment of the peaks at a given genomic feature (e.g. TSS)
#' plotGenomicAnnot(object = cds, plot.type = "enrich")
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @import ggplot2
#' @rdname plotGenomicAnnot-methods
setGeneric(name = "plotGenomicAnnot",
           def = function(object, plot.type = c("distrib", "enrich"), xlab, ylab, color)
                   standardGeneric("plotGenomicAnnot")
)


####------------------- plotFeatures -------------------####
#' plotFeatures
#'
#' Visualize the relations between predictors and response variable
#' ('tssOverlap').
#'
#' @name plotFeatures
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#' @param plot.type One of ["box", "density"]. Default: "box"
#' @param feature Feature to plot. By default, all the features are plotted.
#' @param ncol \code{Numeric}. Arrange individual plots in columns. By default,
#'     the number of columns correspond to the number of features used for
#'     plotting.
#' @param xlab \code{Character}. Title of the x-axis
#' @param ylab \code{Character}. Title of the y-axis.
#' @param color A character vector of length two. Default: ["#E41A1C","#377EB8"].
#' @param alpha Color transparency. In a range [0, 1]. Default: 1.
#'
#' @details
#'    In order to discriminate between functional or gene associated peaks and
#'    non-functional or background peaks, each peak in the data set is
#'    characterized by several features. Moreover, the user might supply her/his
#'    own list of features with the \code{\link{addFeature}}. Prior to fitting
#'    the logistic model, the relations between predictors and response variable
#'    (tssOverlap) can be explored with \code{plotFeatures}. Based on the plots,
#'    poor predictors can be excluded from the analysis to improve the model
#'    fit.
#'
#' @return \code{\link{ggplot2}} object.
#'
#' @seealso \code{\link{constructCDS}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### The data can be plotted in two ways
#' ### As a boxplot
#' plotFeatures(object = cds, plot.type = "box")
#'
#' ### Or as a density plot
#' plotFeatures(object = cds, plot.type = "density")
#'
#' ### Additionally, only the subset of features can be shown
#' plotFeatures(object = cds, plot.type = "box", feature = c("pileup", "length"))
#'
#' ### The position of the graphs on the plot, can be adjusted by 'ncol' argument
#' plotFeatures(object = cds, plot.type = "box", ncol = 2)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname plotFeatures-methods
setGeneric(name = "plotFeatures",
           def = function(object,
                          plot.type = c("box", "density"),
                          feature,
                          ncol,
                          xlab,
                          ylab,
                          color = c("#E41A1C","#377EB8"),
                          alpha = 1)
                   standardGeneric("plotFeatures")
)


####------------------- plotROC -------------------####
#' plotROC
#'
#' Visualize the performance of the classification model fit (prediction of the
#' gene associated peaks).
#'
#' @name plotROC
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#' @param ... Further arguments passed to plot.
#'
#' @details
#'     The \code{plotROC} is a simple wrapper for the \code{plot} function
#'     implemented in \code{\link{pROC}} package.
#'
#'     The \href{https://en.wikipedia.org/wiki/Receiver_operating_characteristic}{ROC}
#'     curve is created by plotting the true positive rate (sensitivity) against
#'     the false positive rate (1 - specificity) at various threshold settings.
#'     The closer the curve follows the left-hand border and then the top border
#'     of the ROC space, the more accurate the test. The area under the curve
#'     (AUC) is a measure of accuracy.
#'
#'
#' @return ROC plot.
#'
#' @seealso \code{\link{predictTssOverlap}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### Classify peaks on gene associated and background
#' predictTssOverlap(object = cds, feature = "pileup", p = 0.75)
#'
#' ### Visualize the performance of the model fit
#' plotROC(object = cds, col = "red3", grid = TRUE, auc.polygon = TRUE)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname plotROC-methods
setGeneric(name = "plotROC",
           def = function(object, ...)
                   standardGeneric("plotROC")
)


####------------------- peaksToBed -------------------####
#' peaksToBed
#'
#' A convenient way to explore output of the predictions made on the ChIP peaks
#' is to visualize them in the
#' \href{https://genome.ucsc.edu/}{UCSC genome browser}. The peaksToBed function
#' returns a file in
#' \href{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}{BED} format, which
#' can be uploaded directly to the genome browser. To improve the visual
#' perception, peaks are color-coded by the predicted strand.
#'
#' @name peaksToBed
#' @docType methods
#'
#' @param object A \code{\link{ChipDataSet}} object.
#' @param file \code{Character}. A file name.
#' @param strand.pred.color \code{Character}. A vector of length four, specifying
#'     colors of the predicted peak strand. The order of colors corresponds
#'     to 1) plus strand, 2) minus strand, 3) bideractional, 4) strand not
#'     determined. Default: c("blue", "red", "green4", "black").
#' @param gene.associated.peaks \code{Logical}. Whether to return gene
#'     associated peaks (based on the prediction) only. Default: TRUE
#'
#' @return A file in the BED format.
#'
#' @seealso
#'     \code{\link{constructCDS}}
#'     \code{\link{predictTssOverlap}}
#'     \code{\link{predictStrand}}
#'     \href{https://genome.ucsc.edu/}{UCSC genome browser}
#'     \href{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}{BED}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Classify peaks on gene associated and background
#' predictTssOverlap(object = cds, feature = "pileup", p = 0.75)
#'
#' ### Predict peak 'strandedness'
#' predictStrand(cdsObj = cds, tdsObj = tds, coverage.cutoff = 5,
#' quant.cutoff = 0.1, win.size = 2500)
#'
#' # peaksToBed(object = cds, file = "peaks.bed", gene.associated.peaks = TRUE)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname peaksToBed-methods
setGeneric(name = "peaksToBed",
           def = function(object, file,
                          strand.pred.color = c("blue", "red", "green4", "black"),
                          gene.associated.peaks = TRUE)
                   standardGeneric("peaksToBed")
)

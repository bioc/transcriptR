####------------------- estimateBackground -------------------####
#' estimateBackground
#'
#' Gene expression is a stochastic process, which often results in substantial
#' expression noise. To obtain a putative set of transcribed regions, it is
#' necessary to identify those regions that are expressed significantly above
#' the background level. Using a Poisson-based approach for estimating
#' the noise distribution from the frequency of the transcribed regions with the
#' low fragments coverage, \code{\link{estimateBackground}} function returns a
#' coverage cutoff value for a specific
#' \href{https://en.wikipedia.org/wiki/False_discovery_rate}{False Discovery Rate (FDR)}.
#'
#' @name estimateBackground
#' @docType methods

#' @param object A \code{\link{TranscriptionDataSet}} object.
#' @param fdr.cutoff \code{Numeric}. False Discovery Rate cutoff value.
#'     Default: 0.05.
#'
#' @return The slots \code{coverageCutoffFdr} and \code{coverageCutoff}
#'     of the provided \code{TranscriptionDataSet} object will be updated by
#'     the FDR cutoff value used in the calculations and by the corresponding
#'     estimated coverage cutoff value, respectively.
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Estimate coverage cutoff at different FDR levels
#' estimateBackground(object = tds, fdr = 0.01)
#'
#' @seealso \code{\link{constructTDS}} \code{\link[chipseq]{peakCutoff}}
#'
#' @author Armen R. Karapetyan
#'
#' @export estimateBackground
#' @rdname estimateBackground-methods
setGeneric(name = "estimateBackground",
           def = function(object, fdr.cutoff = 0.05)
                   standardGeneric("estimateBackground")
)


####------------------- estimateGapDistance -------------------####
#' estimateGapDistance
#'
#' The ultimate goal of \code{transcriptR} is to identify continuous regions
#' of transcription. However, in some areas of the genome it is not possible
#' to detect transcription, because of the presence of the low mappability
#' regions and (high copy number) repeats. Sequencing reads can not be uniquely
#' mapped to these positions, leading to the formation of gaps in otherwise
#' continuous coverage profiles and segmentation of transcribed regions into
#' multiple smaller fragments. The gap distance describes the maximum allowed
#' distance between adjacent fragments to be merged into one transcript. To
#' choose the optimal value for the gap distance, the detected transcripts
#' should largely be in agreement with available reference annotations.
#' To accomplish this, the function is build on the methodology proposed by
#' \href{http://www.sciencedirect.com/science/article/pii/S009286741100376X}{Hah et al. (Cell, 2011)}.
#' In brief, the two types of erros are defined:
#' \itemize{
#'     \item \code{dissected} error - the ratio of annotations that is segmented
#'         into two or more fragments.
#'     \item \code{merged} error - the ratio of non-overlapping annotations that
#'         merged by mistake in the experimental data.
#' }
#' There is an interdependence between two types of errors. Increasing the gap
#' distance decreases the \code{dissected} error, by detecting fewer, but longer
#' transcripts, while the \code{merged} error will increase as more detected
#' transcripts will span multiple annotations. The gap distance with the lowest
#' sum of two error types is chosen as the optimal value.
#'
#' @name estimateGapDistance
#' @docType methods

#' @param object A \code{\link{TranscriptionDataSet}} object.
#' @param annot \code{\link[GenomicRanges]{GRanges}}. Reference annotations.
#' @param coverage.cutoff \code{Numeric}. A cutoff value to discard regions with
#'     the low fragments coverage, representing expression noise. By default,
#'     the value stored in the \code{coverageCutoff} slot of the supplied
#'     \code{TranscriptionDataSet} object is used. The optimal cutoff value can
#'     be calculated by \code{\link{estimateBackground}} function call.
#' @param filter.annot \code{Logical}. Whether to filter out lowly expressed
#'     annotations, before estimating error rates. Default: TRUE.
#' @param fpkm.quantile \code{Numeric}. A number in a range (0, 1). A cutoff value
#'     used for filtering lowly expressed annotations. The value corresponds to
#'     the FPKM quantile estimated for the supplied annotations. Default: 0.25.
#' @param gap.dist.range A numeric vector specifying a range of gap distances
#'     to test. By default, the range is from 0 to 10000 with a step of 100.
#'
#' @return The slot \code{gapDistanceTest} of the provided
#'     \code{TranscriptionDataSet} object will be updated by the
#'     \code{data.frame}, containing estimated error rates for each
#'     tested gap distance (see \code{\link{getTestedGapDistances}}, for the
#'     details).
#'
#' @references
#'     Hah N, Danko CG, Core L, Waterfall JJ, Siepel A, Lis JT, Kraus WL.
#'     A rapid, extensive, and transient transcriptional response to estrogen
#'     signaling in breast cancer cells. Cell. 2011.
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Load reference annotations (knownGene from UCSC)
#' data(annot)
#'
#' ### Estimate gap distance minimazing error rate
#' ### Define the range of gap distances to test
#' gdr <- seq(from = 0, to = 10000, by = 1000)
#'
#' estimateGapDistance(object = tds, annot = annot, coverage.cutoff = 5,
#' filter.annot = FALSE, gap.dist.range = gdr)
#'
#' ### View estimated gap distance
#' tds
#'
#' @seealso \code{\link{constructTDS}} \code{\link{plotErrorRate}}
#'     \code{\link{getTestedGapDistances}}
#'
#' @author Armen R. Karapetyan
#'
#' @export estimateGapDistance
#' @rdname estimateGapDistance-methods
setGeneric(name = "estimateGapDistance",
           def = function(object,
                          annot,
                          coverage.cutoff,
                          filter.annot = TRUE,
                          fpkm.quantile = 0.25,
                          gap.dist.range = seq(from = 0, to = 10000, by = 100))
                   standardGeneric("estimateGapDistance")
)


####------------------- detectTranscripts -------------------####
#' detectTranscripts
#'
#' The function dissects transcribed regions (transcripts) genome-wide and
#' performs expression level quantification.
#'
#' @name detectTranscripts
#' @docType methods
#'
#' @param object A \code{\link{TranscriptionDataSet}} object.
#' @param estimate.params \code{Logical}. Whether to estimate expression level
#'     and coverage density of the detected transcripts. Default: TRUE.
#' @param coverage.cutoff \code{Numeric}. A cutoff value to discard regions with
#'     the low fragments coverage, representing expression noise. By default,
#'     the value stored in the \code{coverageCutoff} slot of the supplied
#'     \code{TranscriptionDataSet} object is used. The optimal cutoff value can
#'     be calculated by \code{\link{estimateBackground}} function call.
#' @param gap.dist \code{Numeric}. Maximum allowed distance between transcribed
#'     regions to be merged into the one transcript. By default, the value
#'     stored in the \code{gapDistanceTest} slot of the supplied
#'     \code{TranscriptionDataSet} object is used. The optimal gap distance can
#'     be calculated by \code{\link{estimateGapDistance}} function call.
#' @param total.reads \code{Numeric}. Total number of reads used for the
#'     normalization, when calculating FPKM. By default, the total number of
#'     reads stored in the provided
#'     \code{\link{TranscriptionDataSet}} object is used.
#' @param annot \code{\link[GenomicRanges]{GRanges}}. Reference annotations.
#' @param combine.by.annot \code{Logical}. Whether to combine transcripts
#'     overlapping the same reference annotation. Default: FALSE.
#'
#' @details
#'     The function uses two parameters to identify transcribed regions:
#'     \code{coverage.cutoff} and \code{gap.dist} as calculated by the
#'     \code{\link{estimateBackground}} and \code{\link{estimateGapDistance}},
#'     respectively and stored in the \code{\link{TranscriptionDataSet}} object.
#'     Alternatively, the user may specify his/her own values to be passed to
#'     the function. By increasing the \code{gap.dist}, fewer transcripts of
#'     longer size will be identified, and an increase in the \code{coverage.cutoff}
#'     will result in fewer transcripts of shorter size (a typical transcript tends
#'     to have a lower fragments coverage at the 3' end, and thus, the
#'     \code{coverage.cutoff} value will have an impact on the resulting length of
#'     the detected transcript).
#'
#'     If \code{estimate.params} is set TRUE, the following metrics are estimated for
#'     each transcript:
#'         \itemize{
#'             \item \code{length} - transcript length (in base pairs).
#'             \item \code{bases.covered} - the number of bases covered by the
#'                   sequencing fragments.
#'             \item \code{coverage} - the proportion of transcript length covered by
#'                   fragments. Value in the range (0, 1].
#'             \item \code{fragments} - total number of fragments per transcript.
#'             \item \code{fpkm} - Fragments Per Kilobase of transcript per Million
#'                   mapped reads.
#'             }
#'
#'     The \code{coverage} is a measure of how densely the transcript is covered by
#'     the sequencing fragments. Modestly/highly expressed transcripts will have
#'     a value close to 1, whereas lowly expressed transcripts will have a
#'     value close to 0, indicating the sparse distribution of sequencing
#'     fragments along the transcript body.
#'
#' @return The slot \code{transcripts} of the provided
#'     \code{TranscriptionDataSet} object will be updated by the
#'     \code{\link[GenomicRanges]{GRanges}} object, containing detected
#'     transcripts and, if estimated, corresponding expression levels.
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Load reference annotations (knownGene from UCSC)
#' data(annot)
#'
#' ### Detect transcripts
#' detectTranscripts(object = tds, coverage.cutoff = 5, gap.dist = 4000,
#' estimate.params = TRUE, combine.by.annot = FALSE, annot = annot)
#'
#' ### View detected transcripts
#' getTranscripts(tds)
#'
#' @author Armen R. Karapetyan
#'
#' @seealso \code{\link{constructTDS}}
#'
#' @export detectTranscripts
#' @rdname detectTranscripts-methods
setGeneric(name = "detectTranscripts",
           def = function(object,
                          coverage.cutoff,
                          gap.dist,
                          estimate.params = TRUE,
                          total.reads,
                          combine.by.annot = FALSE,
                          annot)
                   standardGeneric("detectTranscripts"))


####------------------- annotateTranscripts -------------------####
#' annotateTranscripts
#'
#' Annotate detected transcripts by the available reference annotations based
#' on genomic overlap.
#'
#' @name annotateTranscripts
#'
#' @param object A \code{\link{TranscriptionDataSet}} object.
#' @param annot \code{\link[GenomicRanges]{GRanges}}. Reference annotations.
#' @param min.overlap \code{Numeric}. A minimal proportion of the overlap between
#'     transcript and annotation. A value in the range (0, 1]. Default: 0.3.
#'
#' @return  An extra column \code{annotation.overlap} will be added to the metadata
#'     portion of the \code{\link[GenomicRanges]{GRanges}} object which is
#'     stored in the \code{transcripts} slot of the provided
#'     \code{\link{TranscriptionDataSet}} object.
#'
#' @details
#'    Genomic overlap between transcript and annotation is calculated as the
#'    mean of two proportions: 1) proportion of the transcript length overlapping
#'    annotation; 2) proportion of the annotation length overlapping transcript.
#'    This approach levels off differences in length between transcript
#'    and annotation and, thus better suitable for cases in which the length of
#'    either transcript or annotation is much longer than of compared element.
#'
#'    If there is an overlap between transcript and annotation, the ID of the
#'    associated annotation will be linked to the transcript.
#'
#' @seealso \code{\link{detectTranscripts}}
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Load reference annotations (knownGene from UCSC)
#' data(annot)
#'
#' ### Detect transcripts
#' detectTranscripts(object = tds, coverage.cutoff = 5, gap.dist = 4000,
#' estimate.params = TRUE, combine.by.annot = FALSE, annot = annot)
#'
#' ### Annotate detected transcripts
#' annotateTranscripts(object = tds, annot = annot)
#'
#' ### View detected transcripts and associated annotations
#' getTranscripts(tds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname annotateTranscripts-methods
setGeneric(name = "annotateTranscripts",
           def = function(object, annot, min.overlap = 0.3)
                   standardGeneric("annotateTranscripts")
)


####------------------- getTranscipts -------------------####
#' getTranscripts
#'
#' Retrieve transcripts information from the \code{\link{TranscriptionDataSet}}
#' object.
#'
#' @name getTranscripts
#'
#' @param object A \code{\link{TranscriptionDataSet}} object.
#' @param min.length \code{Numeric}. A minimum length (in bp) of the reported
#'     transcripts.
#' @param min.fpkm \code{Numeric}. A minimum FPKM of the reported transcripts.
#' @param min.coverage \code{Numeric}. A minimum coverage of the reported
#'     transcripts.
#'
#' @details
#'     The coverage is a measure of how densely the transcript is covered by the
#'     sequencing fragments. Modestly/highly expressed transcripts will have a
#'     value close to 1, whereas lowly expressed transcripts will have a value
#'     close to 0, indicating the sparse distribution of sequencing fragments
#'     along the transcript body.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @seealso \code{\link{detectTranscripts}}
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### View detected transcripts
#' getTranscripts(tds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname getTranscripts-methods

setGeneric(name = "getTranscripts",
           def = function(object,
                          min.length,
                          min.fpkm,
                          min.coverage)
                   standardGeneric("getTranscripts")
)


####------------------- getTestedGapDistances -------------------####
#' getTestedGapDistances
#'
#' Retrieve a \code{data.frame}, containing a range of tested gap distances
#' and estimated error rates.
#'
#' @name getTestedGapDistances
#'
#' @param object A \code{\link{TranscriptionDataSet}} object.
#'
#' @return A \code{data.frame} containing estimated error rates (\code{dissected},
#'     \code{merged} and \code{sum of two errors}) and corresponding gap distances.
#'
#' @seealso \code{\link{estimateGapDistance}}
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' head(getTestedGapDistances(tds))
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname getTestedGapDistances-methods

setGeneric(name = "getTestedGapDistances",
           def = function(object)
                   standardGeneric("getTestedGapDistances")
)


####------------------- breakTranscriptsByPeaks -------------------####
#' breakTranscriptsByPeaks
#'
#' The function divides closely spaced transcripts into individually transcribed
#' units using the detected active transcription start sites.
#'
#' @name breakTranscriptsByPeaks
#'
#' @param tdsObj A \code{\link{TranscriptionDataSet}} object.
#' @param cdsObj A \code{\link{ChipDataSet}} object.
#' @param estimate.params \code{Logical}. Whether to estimate expression level
#'     and coverage density of the newly detected transcripts. Default: TRUE.
#'
#' @return The slot \code{transcripts} of the provided
#'     \code{TranscriptionDataSet} object will be updated by the
#'     \code{\link[GenomicRanges]{GRanges}} object, containing transcripts and,
#'     if estimated, corresponding expression levels.
#'
#' @details
#'     One of the challenges for primary transcript detection concerns the
#'     simultaneous transcription of closely spaced genes, which needs to be
#'     properly divided into individually transcribed units. \code{transcriptR}
#'     combines RNA-seq data with ChIP-seq data of histone modifications that
#'     mark active Transcription Start Sites (TSSs), such as, H3K4me3 or H3K9/14Ac
#'     to overcome this challenge. The advantage of this approach over the use of,
#'     for example, gene annotations is that this approach is data driven and
#'     therefore able to deal also with novel and case specific events. Furthermore,
#'     the integration of ChIP- and RNA-seq data allows the identification all
#'     known and novel active transcription start sites within a given sample.
#'     Transcription initiation within a peak region is investigated by comparing
#'     RNA-seq read densities upstream and downstream of empirically determined TSSs.
#'     Closely spaced transcripts are divided into individually transcribed units
#'     using the detected active TSSs.
#'
#' @seealso \code{\link{predictStrand}}
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### Load reference annotations (knownGene from UCSC)
#' data(annot)
#'
#' ### Detect transcripts
#' detectTranscripts(object = tds, coverage.cutoff = 5, gap.dist = 4000,
#' estimate.params = TRUE, combine.by.annot = FALSE, annot = annot)
#'
#' ### Classify peaks on gene associated and background
#' predictTssOverlap(object = cds, feature = "pileup", p = 0.75)
#'
#' ### Predict peak 'strand'
#' predictStrand(cdsObj = cds, tdsObj = tds, coverage.cutoff = 5,
#' quant.cutoff = 0.1, win.size = 2500)
#'
#' ### If `estimate.params = TRUE`, FPKM and coverage density will be re-calculated
#' breakTranscriptsByPeaks(tdsObj = tds, cdsObj = cds, estimate.params = TRUE)
#'
#' ### View detected transcripts
#' getTranscripts(tds)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname breakTranscriptsByPeaks-methods
setGeneric(name = "breakTranscriptsByPeaks",
           def = function(tdsObj, cdsObj, estimate.params = TRUE)
                   standardGeneric("breakTranscriptsByPeaks"))


####------------------- plotErrorRate -------------------####
#' plotErrorRate
#'
#' A simple helper function that plot results of
#' \code{\link{estimateGapDistance}} function call.
#'
#' @name plotErrorRate
#' @param object A \code{\link{TranscriptionDataSet}} object.
#' @param color \code{Character}. Colors to be used to plot estimated
#'     errors.
#' @param xlab \code{Character}. Label of the x-axis.
#' @param ylab \code{Character}. Lable of the y-axis.
#' @param ... Further arguments passed to plot.
#' @return plot
#'
#' @details
#'     The tested gap distances are plotted on the x-axis and corresponding
#'     error rates on the y-axis. Three curved lines depict the two error
#'     types calculated by \code{\link{estimateGapDistance}} and the sum of
#'     both errors. The vertical dashed line depicts the gap distance with
#'     the smallest sum of two errors.
#'
#' @seealso \code{\link{estimateGapDistance}}
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Load reference annotations (knownGene from UCSC)
#' data(annot)
#'
#' ### Estimate gap distance minimazing error rate
#' ### Define the range of gap distances to test
#' gdr <- seq(from = 0, to = 10000, by = 1000)
#'
#' estimateGapDistance(object = tds, annot = annot, coverage.cutoff = 5,
#' filter.annot = FALSE, gap.dist.range = gdr)
#'
#' plotErrorRate(object = tds, lwd = 2)
#'
#' @author Armen R. Karapetyan
#'
#' @rdname plotErrorRate-methods
#' @export
setGeneric(name = "plotErrorRate",
           def = function(object,
                          color = c("#1B9E77", "#D95F02", "#7570B3"),
                          xlab = "Gap distance (kb)",
                          ylab = "Error rate",
                          ...)
                   standardGeneric("plotErrorRate"))


####------------------- exportCoverage -------------------####
#' exportCoverage
#'
#' RNA-seq coverage profiles for both forward and reverse DNA strand can
#' be visualized separately in the
#' \href{https://genome.ucsc.edu/}{UCSC genome browser} using \code{exportCoverage}.
#' This function can generate tracks in
#' \href{https://genome.ucsc.edu/goldenPath/help/bigWig.html}{BigWig} and
#' \href{http://genome.ucsc.edu/goldenpath/help/bedgraph.html}{bedGraph} formats,
#' which can be uploaded to the genome browser as custom tracks.
#'
#' @name exportCoverage
#' @docType methods

#' @param object A \code{\link{TranscriptionDataSet}} object.
#' @param file \code{Character}. A file name.
#' @param strand \code{Character}. The strand to create a track for.
#'     One of ["+", "-"].
#' @param type \code{Character}. Track type, either "bigWig" or "bedGraph".
#'     Default: "bedGraph".
#' @param color Object of class "integer" representing the track color
#'     (as from col2rgb). Only works with tracks of type "bedGraph". Default:
#'     c(0L, 0L, 255L).
#' @param filter.by.coverage.cutoff \code{Logical}. Whether to discard regions
#'     with low fragment coverage, representing expression noise from the
#'     resulting track. Default: FALSE.
#' @param coverage.cutoff \code{Numeric}. A cutoff value to discard regions with
#'     the low fragment coverage, representing expression noise. By default,
#'     the value stored in the \code{coverageCutoff} slot of the supplied object
#'     is used. The optimal cutoff value can be calculated by
#'     \code{\link{estimateBackground}} function call. Default: NULL.
#' @param rpm \code{Logical}. Whether to perform normalization ('Reads Per
#'     Million'). Default: FALSE.
#' @param total.reads \code{Numeric}. Total number of reads used for
#'     normalization. By default, the total number of
#'     reads stored in the provided \code{\link{TranscriptionDataSet}} object is
#'     used.
#'
#' @return A file in either
#'     \href{https://genome.ucsc.edu/goldenPath/help/bigWig.html}{BigWig} or
#'     \href{http://genome.ucsc.edu/goldenpath/help/bedgraph.html}{bedGraph}
#'     format.
#'
#' @details
#'     There is an option to filter coverage profiles by the coverage cutoff
#'     value, either estimated for a specific FDR via
#'     \code{\link{estimateBackground}} or a user specified value. By default, the
#'     coverage cutoff value stored in the \code{\link{TranscriptionDataSet}} object
#'     is used. In order to make an informed decision about a proper FDR level, it
#'     is useful to explore the output at different FDR levels and determine the
#'     optimal cutoff value.
#'
#' @seealso
#'     \code{\link{estimateBackground}}
#'     \href{https://genome.ucsc.edu/}{UCSC genome browser}
#'     \href{https://genome.ucsc.edu/goldenPath/help/bigWig.html}{BigWig}
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Look at the coverage profile of the regions expressed above the background level
#' # exportCoverage(object = tds, file = "plus.bg", type = "bedGraph", strand = "+",
#' # filter.by.coverage.cutoff = TRUE, coverage.cutoff = 3, rpm = FALSE)
#'
#' ### Or check the raw coverage (all expressed regions)
#' # exportCoverage(object = tds, file = "plus_raw.bg", type = "bedGraph",
#' # strand = "+", filter.by.coverage.cutoff = FALSE, rpm = FALSE)
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname exportCoverage-methods
setGeneric(name = "exportCoverage",
           def = function(object,
                          file,
                          type,
                          strand,
                          color,
                          filter.by.coverage.cutoff = FALSE,
                          coverage.cutoff = NULL,
                          rpm = FALSE,
                          total.reads)
                   standardGeneric("exportCoverage")
)


####------------------- transcriptsToBed -------------------####
#' transcriptsToBed
#'
#' A convenient graphical way to explore the identified transcripts
#' is to visualize them in the
#' \href{https://genome.ucsc.edu/}{UCSC genome browser}.
#' The \code{transcriptsToBed} function returns a file in
#' \href{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}{BED}
#' format, which can be directly uploaded to the genome browser.
#' To improve the visual perception, transcripts are color-coded by
#' DNA strand orientation.
#'
#' @name transcriptsToBed
#' @docType methods
#'
#' @param object A \code{\link{GRanges}} object.
#' @param file \code{Character}. A file name.
#' @param strand.color A character vector of length two, specifying color for
#'     each DNA strand. Default: c("blue", "red").
#'
#' @return A file in the BED format.
#'
#' @seealso
#'     \code{\link{estimateBackground}}
#'     \href{https://genome.ucsc.edu/}{UCSC genome browser}
#'     \href{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}{BED}
#'
#' @examples
#' ### Load TranscriptionDataSet object
#' data(tds)
#'
#' ### Load reference annotations (knownGene from UCSC)
#' data(annot)
#'
#' ### Detect transcripts
#' detectTranscripts(object = tds, coverage.cutoff = 5, gap.dist = 4000,
#' estimate.params = TRUE, combine.by.annot = FALSE, annot = annot)
#'
#' ### View detected transcripts
#' trx <- getTranscripts(tds)
#'
#' ### Export to BED
#' # transcriptsToBed(object = trx, file = "transcripts.bed",
#' # strand.color = c("blue", "red"))
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @rdname transcriptsToBed-methods
setGeneric(name = "transcriptsToBed",
           def = function(object, file, strand.color = c("blue", "red"))
                   standardGeneric("transcriptsToBed")
)

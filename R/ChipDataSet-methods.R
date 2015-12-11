####------------------- constructCDS -------------------####
.readPeaks <- function(chr, start, end, name, reduce = FALSE, gapwidth = 0){
        # Function to convert peaks to GRanges
        # Peaks located close to each other can be merged

        # convert to GRanges
        gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

        # merge peaks
        if (reduce){
                gr <- reduce(gr, min.gapwidth = gapwidth, with.revmap = TRUE)
                if (!missing(name)){
                        gr$merged.peaks <- as(lapply(gr$revmap, function(i) name[i]), "CharacterList")
                } else {
                        gr$merged.peaks <- gr$revmap
                }
                gr$revmap <- NULL
        }

        return(gr)
}

evaluatePeaks <- function(peak, TxDb, tss.region = c(-1000, 1000)){

        # function to evaluate peaks enrichment at the specific genomic features, e.g. exons, introns
        # ideally peaks must be enriched at the gene TSS regions
        if (is.numeric(tss.region)){
                if (length(tss.region) > 2){
                        stop("[ERROR] 'tss.region' too many values provided.")
                }
        } else {
                stop("[ERROR] 'tss.region' needs to be numeric vector.")
        }

        # Extract genomic regions from TxDb
        message("[INFO] Preparing genomic annotations from TxDb...", appendLF = FALSE)

        transcripts.gr <- sort.GenomicRanges(GenomicFeatures::transcripts(TxDb), ignore.strand = TRUE)
        # promoter region
        tss.gr <- GenomicFeatures::promoters(transcripts.gr, upstream = abs(tss.region[1]), downstream = abs(tss.region[2]))
        tss.gr <- reduce(tss.gr, ignore.strand = TRUE)
        # exons
        exons.gr <- reduce(GenomicFeatures::exons(TxDb), ignore.strand = TRUE)
        exons.gr <- setdiff(exons.gr, tss.gr, ignore.strand = TRUE)
        # introns
        introns.gr <- reduce(unlist(GenomicFeatures::intronsByTranscript(TxDb)), ignore.strand = TRUE)
        introns.gr <- setdiff(introns.gr, tss.gr, ignore.strand = TRUE)
        introns.gr <- setdiff(introns.gr, exons.gr, ignore.strand = TRUE)

        # intergenic regions
        intergen.gr <- gaps(c(tss.gr, exons.gr, introns.gr))
        intergen.gr <- intergen.gr[strand(intergen.gr) == "*"]

        message("Done!")


        message("[INFO] Estimating distribution of the peaks...", appendLF = FALSE)
        # Calculate observed peak's distribution
        # Use peak center to search for overlap peaks and genomic features
        centre <- start(peak) + (width(peak)/2)
        peak.centre <- GRanges(seqnames = seqnames(peak), ranges = IRanges(start = centre, width = 1))

        regions <- list(tss = tss.gr, exons = exons.gr, introns = introns.gr, intergenic = intergen.gr)
        observed <- rep(NA, 4)
        names(observed) <- names(regions)
        for (i in seq_along(regions)) {
                hits <- findOverlaps(peak.centre, reduce(regions[[i]], ignore.strand = TRUE), ignore.strand = TRUE)
                observed[i] <- length(unique(queryHits(hits)))
        }

        # Calculate expected peak's distribution
        tss.len <- as.numeric(sum(width(tss.gr)))
        exons.len <- as.numeric(sum(width(exons.gr)))
        introns.len <- as.numeric(sum(width(introns.gr)))
        intergen.len <- as.numeric(sum(width(intergen.gr)))
        total.len <- as.numeric(sum(c(tss.len, exons.len, introns.len, intergen.len)))

        tss.prop <- round(length(peak) * (tss.len/total.len), 0)
        exons.prop <- round(length(peak) * (exons.len/total.len), 0)
        introns.prop <- round(length(peak) * (introns.len/total.len), 0)
        intergen.prop <- round(length(peak) * (intergen.len/total.len), 0)

        expected <- c(tss.prop, exons.prop, introns.prop, intergen.prop)

        # Prepare output
        region <- c(paste0("TSS [", paste(tss.region, collapse = ":"), "]"), "Exons", "Introns", "Intergenic")
        distrib <- data.frame(region = region,
                              observed = observed,
                              expected = expected,
                              ratio = observed/expected)
        rownames(distrib) <- NULL
        distrib$region <- factor(distrib$region, levels = distrib$region)
        message("Done!")

        return(distrib)
}

#' constructCDS
#'
#' The function constructs an object of class \code{\link{ChipDataSet}}, which
#' is a container for holding processed sequencing data and the results of
#' all downstream analyses. All the slots of the created object are filled
#' during the workflow by applying specific functions to the object directly.
#'
#' @name constructCDS
#' @param peaks A path to a file with peaks. The file needs to have at least 3
#'     columns (tab-separated): chromosome, start (peak), end (peak). The 4th
#'     column - name (peak id) is optional.
#' @param reads A path to a BAM file with sequencing reads.
#' @param region \code{\link[GenomicRanges]{GRanges}}. Genomic region(s) to
#'     extract reads from. If not supplied, all the reads from a BAM file are
#'     extracted.
#' @param TxDb \code{\link[GenomicFeatures]{TxDb}} object.
#' @param tssOf \code{Character}. Extract Transcription Start Site (TSS)
#'     regions from either "gene" or "transcript" annotations. Default: "gene".
#' @param tss.region A numeric vector of length two, which specifies the size of
#'     TSS region. Default: -2kb to 2kb.
#' @param reduce.peaks \code{Logical}. Whether to merge neighboring peaks.
#'     Default: FALSE.
#' @param gapwidth \code{Numeric}. A minimum distance (in bp) between peaks to
#'     merge. Default: 1000.
#' @param fragment.size \code{Numeric}. Extend read length to the fragment size.
#' @param unique \code{Logical}. Whether to remove duplicated reads (based on
#'     the genomic coordinates). Default: FALSE.
#' @param swap.strand \code{Logical}. Whether to reverse the strand of the read.
#'     Default: FALSE.
#' @param param \code{\link[Rsamtools]{ScanBamParam}} object influencing what
#'     fields and which records (reads) are imported from the Bam file.
#'     Default: NULL.
#'
#' @details
#'     The function \code{constructCDS} initializes a
#'     \code{\link{ChipDataSet}} object, by providing the paths to the input
#'     files and information relevant to the ChIP-seq library preparation
#'     procedure. During the object construction the following steps are
#'     executed:
#'     \itemize{
#'         \item The peak information is converted into the object of
#'             \code{\link[GenomicRanges]{GRanges}} class.
#'         \item The genomic distribution of the peaks is evaluated (exonic,
#'             intronic, intergenic, TSSs).
#'         \item Each peak in the data set is functionally characterized:
#'             \itemize{
#'                 \item \code{length} - the length of a peak (in base pairs).
#'                 \item \code{fragments} - total number of fragments overlapping
#'                     a peak region.
#'                 \item \code{density} - number of fragments per base pair of
#'                     the peak length.
#'                 \item \code{pileup} - highest fragment pileup in each peak
#'                     region.
#'                 \item \code{tssOverlap} - overlap (binary, yes/no) of the
#'                     peak with the annotated TSS region.
#'                     }
#'     The estimated features are used to predict which of the peaks are gene
#'     associated in the analysis downstream.
#'     }
#'
#'     As many peak-calling algorithms tend to divide broader peaks into the
#'     several narrower closely spaced peaks, it is advised to merge these
#'     end-to-end peaks to decrease the number of false positives and prevent
#'     unnecessary truncation of transcripts in the downstream analysis.
#'
#' @return An object of class \code{\link{ChipDataSet}}.
#'
#' @seealso \code{\link{ChipDataSet}} \code{\link{predictTssOverlap}}
#'
#' @examples
#' ### Load ChipDataSet object
#' data(cds)
#'
#' ### View a short summary of the object
#' cds
#'
#' @author Armen R. Karapetyan
#'
#' @export
#' @importClassesFrom GenomicFeatures TxDb
constructCDS <- function(peaks, reads, region, TxDb,
                         tssOf = c("gene", "transcript"),
                         tss.region = c(-2000, 2000),
                         reduce.peaks = FALSE, gapwidth = 1000,
                         fragment.size, unique = TRUE,
                         swap.strand = FALSE, param = NULL){

        # Tests
        if (!file.exists(peaks)){
                stop("[ERROR] 'peaks' file does not exists.")
        }

        if (!file.exists(reads)){
                stop("[ERROR] 'reads' file does not exists.")
        }

        if (!missing(region)){
                if (!is(region, "GRanges")){stop("[ERROR] 'region' needs to be of class GRanges.")}
        }

        if (!is(TxDb, "TxDb")){
                stop("[ERROR] 'TxDb' needs to be of class 'TxDb'.")
        }

        if (is.numeric(tss.region)){
                if (length(tss.region) > 2){
                        stop("[ERROR] 'tss.region' too many values provided.")
                }
        } else {
                stop("[ERROR] 'tss.region' needs to be a numeric vector of length two.")
        }

        if (!is.logical(reduce.peaks)){
                stop("[ERROR] 'reduce.peaks' needs to be a logical.")
        }

        if (!is.numeric(gapwidth)){
                stop("[ERROR] 'gapwidth' needs to be a numeric.")
        }

        if (missing(tssOf)){
                tssOf <- "gene"
        } else {
                tssOf <- match.arg(arg = tssOf, choices = c("gene", "transcript"))
        }

        message("[INFO] Extracting peaks information...", appendLF = FALSE)
        msg <- FALSE
        peaks.df <- utils::read.table(peaks, header = FALSE, sep = "\t")
        if (base::ncol(peaks.df) < 3){
                stop("[ERROR] 'peaks' file must contain at least three columns: chromosome, peak start, peak end.")
        }
        chr <- peaks.df[, 1]
        start <- peaks.df[, 2]
        end <- peaks.df[, 3]
        if (base::ncol(peaks.df) > 3){
                name <- peaks.df[, 4]
                if (!any(duplicated(name))){
                        peaks.gr <- .readPeaks(chr = chr, start = start,
                                               end = end, name = name,
                                               reduce = reduce.peaks,
                                               gapwidth = gapwidth)
                } else {
                        msg <- TRUE
                        peaks.gr <- .readPeaks(chr = chr, start = start,
                                               end = end, reduce = reduce.peaks,
                                               gapwidth = gapwidth)
                }
        } else {
                peaks.gr <- .readPeaks(chr = chr, start = start,
                                       end = end, reduce = reduce.peaks,
                                       gapwidth = gapwidth)
        }
        message("Done!")

        if (msg){message("[WARNING] Duplicated names are detected. Therefore, they will not be used.")}
        if (!missing(region)){
                peaks.gr <- subsetByOverlaps(peaks.gr, region)
                peaks.gr <- GenomeInfoDb::keepSeqlevels(x = peaks.gr, value = GenomeInfoDb::seqlevelsInUse(peaks.gr))
        } else {
                region <- GRanges()
        }
        peaks.gr$id <- paste0("peak_", seq_along(peaks.gr))

        # Evaluate genomic distriution of the peaks
        suppressWarnings(genomic.dist <- evaluatePeaks(peaks.gr, TxDb, tss.region))

        message("[INFO] Extracting fragments...", appendLF = FALSE)
        if (!missing(fragment.size)){
                # extend peak region by the fragment size on each side
                peak.region <- GRanges(seqnames = seqnames(peaks.gr),
                                  ranges = IRanges(start(peaks.gr) - fragment.size,
                                                   end = end(peaks.gr) + fragment.size))
                fragments <- .readBam(file = reads, region = peak.region,
                                      fragment.size = fragment.size,
                                      unique = unique, paired.end = FALSE,
                                      swap.strand = swap.strand)
        } else {
                peak.region <- peaks.gr
                fragments <- .readBam(file = reads, region = peak.region,
                                      unique = unique, paired.end = FALSE,
                                      swap.strand = swap.strand)
        }
        message("Done!")

        # Create coverage profile (non-strand specific)
        message("[INFO] Creating coverage profile...", appendLF = FALSE)
        cov <- coverage(fragments)
        message("Done!")

        # Estimate features
        message("[INFO] Estimating features: ")

        message("       # peak length")
        length <- width(peaks.gr)

        message("       # total number of fragments per peak")
        fragments <- countOverlaps(peaks.gr, fragments)

        message("       # density (fragments/bp)")
        density <- round(fragments/length, 3)

        message("       # pileup (peak height)")
#         chrs <- names(cov)
#         chrs <- chrs[chrs %in% levels(seqnames(peaks.gr))]
        chrs <- GenomeInfoDb::seqlevelsInUse(peaks.gr)
        pileup <- lapply(chrs,
                         function(chr){
                                 peak.sub <- peaks.gr[seqnames(peaks.gr) == chr]
                                 views <- Views(cov[[chr]],
                                                start = start(peak.sub),
                                                end = end(peak.sub))
                                 viewMaxs(views)
                         }
        )
        pileup <- unlist(pileup)

        # Combine into the data frame
        df.features <- data.frame(length, fragments, density, pileup)

        # Mark peaks overlaping with the tss regions
        message("       # overlap with tss regions")

        if (tssOf %in% "gene"){
                annot <- GenomicFeatures::genes(TxDb)
        } else if (tssOf %in% "transcript"){
                annot <- GenomicFeatures::transcripts(TxDb)
        } else {
                stop("[ERROR] Wrong 'tssOf', needs to be either 'gene' or 'transcript'." )
        }

        suppressWarnings(tss.gr <- GenomicFeatures::promoters(annot, upstream = abs(tss.region[1]), downstream = abs(tss.region[2])))
        tss.gr <- trim(tss.gr)
        hits <- findOverlaps(reduce(tss.gr, ignore.strand = FALSE), peaks.gr, ignore.strand = TRUE)
        tssOverlap <- rep("no", length(peaks.gr))
        tssOverlap[unique(subjectHits(hits))] <- "yes"
        tssOverlap <- factor(tssOverlap, levels  = c("yes", "no"))
        df.features <- cbind(tssOverlap, df.features)

        return(new("ChipDataSet",
                   peaks = peaks.gr,
                   region = region,
                   genomicAnnotation = genomic.dist,
                   features = df.features))
}


####------------------- predictTssOverlap -------------------####
#' @rdname predictTssOverlap-methods
setMethod("predictTssOverlap",
          signature = signature("ChipDataSet"),
          definition = function(object, feature, p = 0.8){

                  objName <- deparse(substitute(object))
                  df <- object@features

                  if (!missing(feature)){
                          if (length(which(colnames(df) %in% feature)) == length(feature)){
                                  col.select <- which(colnames(df) %in% feature)
                                  df <- df[, c(1, col.select)]  # first column contain tssOverlap information
                          } else {
                                  stop("[ERROR] Wrong 'feature'.")
                          }
                  }

#                   # All features used in modeling must be numeric
#                   if (!all(apply(df[, -ncol(df), drop = FALSE], MARGIN = 2, is.numeric))){
#                           stop("[ERROR] Not numeric features selected.")
#                   }

                  message("[INFO] Partitioning data...", appendLF = FALSE)
                  # Partition data into the training and test data sets
                  trainIndex <- caret::createDataPartition(y = df[, 1],
                                                           p = p,
                                                           list = FALSE)

                  # Create training and test data sets
                  training <- df[trainIndex, ]
                  testing <- df[-trainIndex, ]
                  message("Done!")

                  # Train the model
                  message("[INFO] Training model: ")
                  # k-fold cross validation
                  cntr <- caret::trainControl(method = "repeatedcv",
                                              number = 10,
                                              repeats = 5,
                                              summaryFunction = caret::twoClassSummary,
                                              classProbs = TRUE)

                  # Fit logistic regression to the trainind data set
                  suppressMessages(logFit <- caret::train(x = training[, -1, drop = FALSE],
                                                          y = training[, 1],
                                                          method = "glm",
                                                          metric = "ROC",
                                                          trControl = cntr))

                  # Evaluate model on the training data set
                  trainPred.probs <- stats::predict(logFit, training, type = "prob")
                  trainPred <- stats::predict(logFit, training)
                  trainConfMat <- caret::confusionMatrix(trainPred, training[, 1])
                  trainRocCurve <- pROC::roc(response = training[, 1],
                                             predictor = trainPred.probs$yes,
                                             levels = levels( training[, 1]))

                  message("       Accuracy - ", round(trainConfMat$overall[1], 4))
                  message("       Area under the curve (AUC) - ", round(trainRocCurve$auc, 4))


                  # Test model on the testing data set
                  message("[INFO] Testing model: ")
                  testPred.probs <- stats::predict(logFit, testing, type = "prob")
                  testPred <- stats::predict(logFit, testing)
                  testConfMat <- caret::confusionMatrix(testPred, testing[, 1])
                  testRocCurve <- pROC::roc(response = testing[, 1],
                                            predictor = testPred.probs$yes,
                                            levels = levels( testing[, 1]))

                  message("       Accuracy - ", round(testConfMat$overall[1], 4))
                  message("       Area under the curve (AUC) - ", round(testRocCurve$auc, 4))

                  # Apply model to the whole data set
                  message("[INFO] Classifying peaks: ")
                  allPred.prob <- stats::predict(logFit, df, type = "prob")
                  allPred <- stats::predict(logFit, df)
                  allConfMat <- caret::confusionMatrix(allPred, df[, 1])
                  allRocCurve <- pROC::roc(response = df[, 1],
                                           predictor = allPred.prob$yes,
                                           levels = levels(df[, 1]))

                  message("       Accuracy - ", round(allConfMat$overall[1], 4))
                  message("       Area under the curve (AUC) - ", round(allRocCurve$auc, 4))

                  # Arrange output
                  pred <- data.frame(predicted.tssOverlap.prob = allPred.prob$yes, predicted.tssOverlap = allPred)

                  object@tssOverlapPrediction$predicted.tssOverlap <- pred
                  object@tssOverlapPrediction$confusionMatrix <- allConfMat
                  object@tssOverlapPrediction$logFitSummary <- stats::coef(summary(logFit))
                  object@tssOverlapPrediction$roc <- allRocCurve

                  assign(objName, object, envir = parent.frame())
          }
)


####------------------- addFeature -------------------####
#' @rdname addFeature-methods
setMethod("addFeature",
          signature = signature("ChipDataSet"),
          definition = function(object, feature){

                  if (length(object@features) != 0){
                          objName <- deparse(substitute(object))
                          object.features <- object@features
                          object.peaks <- object@peaks
                  }

                  if (!is.list(feature)){
                          stop("[ERROR] 'feature' needs to be a list.")
                  }
                  if (length(names(feature)) == 0){
                          stop("[ERROR] 'feature' needs to be a named list.")
                  }
                  if (length(names(feature)) != length(feature)){
                          stop("[ERROR] Not all the elemets of 'feature' list are named.")
                  }
                  if (!all(unlist(lapply(feature, function(x) length(x) == length(object.peaks) )))){
                          stop("[ERROR] Length of each element in the 'feature' list needs to be equal to the number of peaks stored in 'object'.")
                  }
                  if (any(duplicated(names(feature)))){
                          stop("[ERROR] Duplicated names in the supplied 'feature' list.")
                  }
                  if (any(names(feature) %in% names(object.features))){
                          stop("[ERROR] Repeated names detected.")
                  }

                  feature.df <- as.data.frame(feature)
                  object@features <- cbind(object.features, feature.df)
                  assign(objName, object, envir = parent.frame())
          }
)


####------------------- predictStrand -------------------####
# Plot coverege and cummulative coverage
.plotPeakCoverage <- function(cdsObj, tdsObj, id, quant.cutoff = 0.1, strand){


        peak <- cdsObj@peaks[id]
        frags.per.peak <- subsetByOverlaps(query = tdsObj@fragments, subject = peak)
        starts.per.peak <- resize(frags.per.peak, width = 1, fix = "start")

        if (strand %in% "+"){
                if (length(frags.per.peak[strand(frags.per.peak) == strand]) > 0){
                        cov.frags <- tdsObj@coveragePlus[[as.character(peak@seqnames)]]
                        cov.starts <- coverage(starts.per.peak[strand(starts.per.peak) == "+"])[[as.character(peak@seqnames)]]
                } else {
                        stop("[ERROR] No fragments available. Nothing to plot.")
                }

        } else if (strand %in% "-"){
                if (length(frags.per.peak[strand(frags.per.peak) == strand]) > 0){
                        cov.frags <- tdsObj@coverageMinus[[as.character(peak@seqnames)]]
                        cov.starts <- coverage(starts.per.peak[strand(starts.per.peak) == "-"])[[as.character(peak@seqnames)]]
                } else {
                        stop("[ERROR] No fragments available. Nothing to plot.")
                }
        } else {
                stop("[ERROR] Wrong 'strand.'")
        }

        #fragments coverage
        view.frags <- Views(subject = cov.frags, start = start(peak), end = end(peak))[[1]]
        # starts coverage
        view.starts <- Views(subject = cov.starts, start = start(peak), end = end(peak))[[1]]

        graphics::par(las = 1, mar = c(5.1, 4.1, 4.1, 4.1))

        if (strand %in% "+"){
                cumsum.starts <-  cumsum(view.starts)
                cum.starts <- cumsum.starts/max(cumsum.starts)
                coord <- which(cum.starts >= quant.cutoff)[1]
                if (length(S4Vectors::runValue(cum.starts)) == 1){
                        if ( is.na(S4Vectors::runValue(cum.starts)) ){
                                stop("[ERROR] No fragments starts in this regions.")
                        }
                }

                graphics::plot(view.frags, type = "l", lwd = 2, col = "blue", xlab = "Distance (bp)", ylab = "")
                graphics::mtext("RNA-seq coverage", side = 2, line = 3, cex.lab = 1,las = 3, col = "blue")
                graphics::par(new = TRUE)
                graphics::plot(cum.starts,
                     yaxt = "n", ylab = "", xlab = "", col = "red",
                     type = "l", lty = "dashed", lwd = 2)
                graphics::axis(side = 4, at = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1))
                graphics::mtext("Cumulative coverage", side = 4, line = 3, cex.lab = 1,las = 3, col = "red")
                graphics::abline(h = quant.cutoff, v = coord, lty = "dashed", lwd = 1.5)

        } else if (strand %in% "-"){
                cumsum.starts <-  cumsum(rev(view.starts))
                cum.starts <- rev(cumsum.starts/max(cumsum.starts))
                coord <- which(cum.starts <= quant.cutoff)[1]
                if (length(S4Vectors::runValue(cum.starts)) == 1){
                        if ( is.na(S4Vectors::runValue(cum.starts)) ){
                                stop("[ERROR] No fragments starts in this regions.")
                        }
                }

                graphics::par(mar = c(5.1, 4.1, 4.1, 4.1))
                graphics::plot((-1)*view.frags, type = "l", lwd = 2, col = "blue", xlab = "Distance (bp)", ylab = "")
                graphics::mtext("RNA-seq coverage", side = 2, line = 3, cex.lab = 1,las = 3, col = "blue")
                graphics::par(new = TRUE)
                graphics::plot(cum.starts*-1,
                     yaxt = "n", ylab = "", xlab = "", col = "red",
                     type = "l", lty = "dashed", lwd = 2)
                graphics::axis(side = 4, at = (-1)*seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1))
                graphics::mtext("Cumulative coverage", side = 4, line = 3, cex.lab = 1,las = 3, col = "red")
                graphics::abline(h = (-1)*quant.cutoff, v = coord, lty = "dashed", lwd = 1.5)
        }
}

# Estimate transcription start site inside the peak (based on RNA-seq cumulative coverage)
.q1q2SeparatingLine <- function(peaks, cov.frags, cov.starts, cov.cutoff, cum.cov.cutoff, strand){


        chrs <- GenomeInfoDb::seqlevelsInUse(peaks)
        views.frags.rvl <- RleViewsList(lapply(names(cov.frags[chrs]),
                                               function(x){
                                                       Views(subject = cov.frags[[x]],
                                                             start = start(peaks[seqnames(peaks) == x]),
                                                             end = end(peaks[seqnames(peaks) == x]))
                                               }))

        views.starts.rvl <- RleViewsList(lapply(names(cov.starts[chrs]),
                                            function(x){
                                                    Views(subject = cov.starts[[x]],
                                                          start = start(peaks[seqnames(peaks) == x]),
                                                          end = end(peaks[seqnames(peaks) == x]))
                                            }))

        .q1q1Pos <- function(x, strand){
                if (strand %in% "+"){cumsum.x <-  cumsum(x)}
                else if (strand %in% "-"){cumsum.x <-  cumsum(rev(x))}
                cumsum.x <- cumsum.x/max(cumsum.x)
                which(cumsum.x > cum.cov.cutoff)[1]
        }
        q1q2.pos <- viewApply(views.starts.rvl, function(x) .q1q1Pos(x, strand))

        # transcription start
        if (strand %in% "+"){sepline.coord <- start(peaks) + unlist(q1q2.pos)}
        else if (strand %in% "-"){sepline.coord <- end(peaks) - unlist(q1q2.pos)}

        # maximal RNA coverage per peak
        max.cov <- unlist(viewMaxs(views.frags.rvl))

        # peak passing the coverage treshold
        pass.cov.cutoff <- rep("no", length(peaks))
        pass.cov.cutoff[which(max.cov > cov.cutoff)] <- "yes"

        df <- data.frame(max.cov, pass.cov.treshold = pass.cov.cutoff, q1q2.sepline.coord = sepline.coord)
        return(df)
}

# Calculate probability of the read being sampled from q1 or q2
.calculateProbs <- function(peaks, frags, trans.start.pos, win.size, strand){

        df <- data.frame(q1.coord = rep(NA, length(peaks)),
                         q2.coord = rep(NA, length(peaks)),
                         q1.count = rep(NA, length(peaks)),
                         q2.count = rep(NA, length(peaks)))

        # Select peaks with transcription start infromation available
        idx <- !is.na(trans.start.pos)

        # Create quadrants (transcription start coordinate) +- win.size
        if (strand %in% "+"){
                q1 <- GRanges(seqnames = seqnames(peaks[idx]), strand = strand,
                              ranges = IRanges(start = trans.start.pos[idx] - win.size - 1,
                                               end = trans.start.pos[idx] - 1))

                q2 <- GRanges(seqnames = seqnames(peaks[idx]), strand = strand,
                              ranges = IRanges(start = trans.start.pos[idx],
                                               end = trans.start.pos[idx] + win.size))
        } else if (strand %in% "-"){
                q1 <- GRanges(seqnames = seqnames(peaks[idx]), strand = strand,
                              ranges = IRanges(start = trans.start.pos[idx] - win.size,
                                               end = trans.start.pos[idx]))

                q2 <- GRanges(seqnames = seqnames(peaks[idx]), strand = strand,
                              ranges = IRanges(start = trans.start.pos[idx] + 1,
                                               end = trans.start.pos[idx] + win.size + 1))
        }

#
#         q1 <- GRanges(seqnames = seqnames(peaks[idx]), strand = strand,
#                       ranges = IRanges(start = trans.start.pos[idx] - win.size - 1,
#                                        end = trans.start.pos[idx]))
#
#         q2 <- GRanges(seqnames = seqnames(peaks[idx]), strand = strand,
#                       ranges = IRanges(start = trans.start.pos[idx],
#                                        end = trans.start.pos[idx] + win.size - 1))
#
#         # Create quadrants (transcription start coordinate) +- win.size
#         q1 <- GRanges(seqnames = seqnames(peaks[idx]), strand = strand,
#                       ranges = IRanges(start = trans.start.pos[idx] - win.size - 1,
#                                        end = trans.start.pos[idx] - 1))
#
#         q2 <- GRanges(seqnames = seqnames(peaks[idx]), strand = strand,
#                       ranges = IRanges(start = trans.start.pos[idx],
#                                        end = trans.start.pos[idx] + win.size))


        if (strand %in% "+"){
                df$q1.coord[idx] <- paste0(as.character(seqnames(q1)), ":", start(q1), "-", end(q1))
                df$q2.coord[idx] <- paste0(as.character(seqnames(q1)), ":", start(q2), "-", end(q2))
        } else {
                df$q1.coord[idx] <- paste0(as.character(seqnames(q1)), ":", start(q2), "-", end(q2))
                df$q2.coord[idx] <- paste0(as.character(seqnames(q1)), ":", start(q1), "-", end(q1))
        }

        # Count reads in quadrants
        q1.count <- countOverlaps(query = q1, subject = frags[strand(frags) == strand])
        q2.count <- countOverlaps(query = q2, subject = frags[strand(frags) == strand])

        # calculate probabilities
        if (strand %in% "+"){
                df$q1.count[idx] <- q1.count
                df$q2.count[idx] <- q2.count
        } else {
                df$q1.count[idx] <- q2.count
                df$q2.count[idx] <- q1.count
        }

        df$q1.prob <- df$q1.count/rowSums(df[, 3:4])
        df$q2.prob <- 1 - df$q1.prob
        return(df)
}

# Estimate probability cutoff
.estimateProbCutoff <- function(df){
        # df - data.frame, 2 columns: 1) values, 2) classes
        bins <- seq(from = 0, to = 1, by = 0.01)
        fdr <- rep(NA, length(bins))
        fnr <- rep(NA, length(bins))
        for (i in seq_along(bins)){
                a <- table(df[ df[, 1] >= bins[i], 2])
                b <- table(df[, 2]) - a
                confmat <- rbind(a, b)
                tp <- confmat[1,1]
                tn <- confmat[2,2]
                fp <- confmat[1,2]
                fn <- confmat[2,1]
                fdr[i] <- fp / (tp + fp)
                fnr[i] <- fn / (tp + fn)
        }
        prob.cutoff <- bins[which(fnr >= fdr)[1]]
        out <- list(probability.cutoff = prob.cutoff, fdr = fdr, fnr = fnr)
}

#' @rdname predictStrand-methods
setMethod("predictStrand",
          signature = signature("ChipDataSet", "TranscriptionDataSet"),
          definition = function(cdsObj, tdsObj, coverage.cutoff, quant.cutoff = 0.1, win.size = 2500, prob.cutoff){

                  if (is.null(cdsObj@tssOverlapPrediction$predicted.tssOverlap)){
                          stop("[ERROR] 'tssOverlapPrediction' slot is empty.")
                  }

                  if (missing(coverage.cutoff)){
                          if (length(tdsObj@coverageCutoff) == 1){
                                  coverage.cutoff <- tdsObj@coverageCutoff
                          } else {
                                  stop("[ERROR] 'coverageCutoff' slot is empty.
                                       Please, specify 'coverage.cutoff' manually or call estimateBackground function on the object of class TranscriptionDataSet.")
                          }
                  } else {
                          if (is.numeric(coverage.cutoff)){
                                  if (coverage.cutoff <= 0){
                                          stop("[ERROR] 'coverage.cutoff' needs to be above 1.")
                                  }
                          } else {
                                  stop("[ERROR] 'coverage.cutoff' needs to be numeric.")
                          }
                  }

                  if (!missing(quant.cutoff)){
                          if (is.numeric(quant.cutoff)){
                                  if (quant.cutoff <= 0 | quant.cutoff >= 1){
                                          stop("[ERROR] 'quant.cutoff' needs to be in a range (0, 1).")
                                  }
                          } else {
                                  stop("[ERROR] 'quant.cutoff' needs to be numeric.")
                          }
                  }

                  if (!missing(prob.cutoff)){
                          if (is.numeric(prob.cutoff)){
                                  if (prob.cutoff <= 0 | prob.cutoff >= 1){
                                          stop("[ERROR] 'prob.cutoff' needs to be in a range (0, 1).")
                                  }
                          } else {
                                  stop("[ERROR] 'prob.cutoff' needs to be numeric.")
                          }
                  }


                  objName <- deparse(substitute(cdsObj))
                  peaks <- cdsObj@peaks

                  # Prepare coverage profiles for fragments starts
                  message("[INFO] Creating coverage profiles for the starts of the fragments...", appendLF = FALSE)

                  starts <- resize(tdsObj@fragments, width = 1, fix = "start")
                  cov.s.pls <- coverage(starts[strand(starts) == "+"])
                  cov.s.mns <- coverage(starts[strand(starts) == "-"])

                  message("Done!")

                  # Estimate transcription start inside peaks (based on cummulative coverage profile)
                  message("[IFNO] Estimating transcription start position...", appendLF = FALSE)

                  sepline.plus <- .q1q2SeparatingLine(peaks = peaks,
                                                      cov.frags = tdsObj@coveragePlus,
                                                      cov.starts = cov.s.pls,
                                                      cov.cutoff = coverage.cutoff,
                                                      cum.cov.cutoff = quant.cutoff,
                                                      strand = "+")

                  sepline.minus <- .q1q2SeparatingLine(peaks = peaks,
                                                       cov.frags = tdsObj@coverageMinus,
                                                       cov.starts = cov.s.mns,
                                                       cov.cutoff = coverage.cutoff,
                                                       cum.cov.cutoff = quant.cutoff,
                                                       strand = "-")

                  # check if there are peak which have enough fragments to pass the coverage treshold
                  # but there are no fragments starts within the peak region to estimate q1q2.sepline.coord
                  peaks.nostarts.plus <- which(sepline.plus$pass.cov.treshold == "yes" &
                                                       is.na(sepline.plus$q1q2.sepline.coord ))
                  peaks.nostarts.minus <- which(sepline.minus$pass.cov.treshold == "yes" &
                                                        is.na(sepline.minus$q1q2.sepline.coord ))
                  if (length(peaks.nostarts.plus) > 0){
                          sepline.plus[peaks.nostarts.plus, "pass.cov.treshold"] <- "no"
                  }
                  if (length(peaks.nostarts.minus) > 0){
                          sepline.minus[peaks.nostarts.minus, "pass.cov.treshold"] <- "no"
                  }

                  message("Done!")

                  # Calculate probabilities for reads being sampled from either q1 or q2 quadrants
                  # quadrand q1: transcription start coordinate - windows size
                  # quadrand q2: transcription start coordinate + windows size
                  message("[INFO] Calculating probabilities...", appendLF = FALSE)

                  probs.plus <- .calculateProbs(peaks = peaks,
                                                frags = starts,
                                                trans.start.pos = sepline.plus$q1q2.sepline.coord,
                                                win.size = win.size, strand = "+")

                  probs.minus <- .calculateProbs(peaks = peaks,
                                                 frags = starts,
                                                 trans.start.pos = sepline.minus$q1q2.sepline.coord,
                                                 win.size = win.size, strand = "-")

                  message("Done!")

                  idx.plus <- which(sepline.plus$pass.cov.treshold == "yes" & !is.na(sepline.plus$q1q2.sepline.coord))
                  idx.minus <- which(sepline.minus$pass.cov.treshold == "yes" & !is.na(sepline.minus$q1q2.sepline.coord))
                  pred.yes <- which(cdsObj@tssOverlapPrediction$predicted.tssOverlap$predicted.tssOverlap == "yes")
                  pred.no <- which(cdsObj@tssOverlapPrediction$predicted.tssOverlap$predicted.tssOverlap == "no")

                  # Estimate P(q2) cutoff
                  message("[INFO] Estimating probability cutoff value...", appendLF = FALSE)

                  x <- c(probs.plus$q2.prob[ intersect(idx.plus, pred.yes) ],
                         probs.minus$q2.prob[ intersect(idx.minus, pred.yes) ])
                  y <- c(probs.plus$q2.prob[ intersect(idx.plus, pred.no) ],
                         probs.minus$q2.prob[ intersect(idx.minus, pred.no) ])

                  x <- data.frame(prob = x, class = "gene")
                  y <- data.frame(prob = y, class = "background")
                  df <- rbind(x, y)
                  df <- df[stats::complete.cases(df), ]
                  prob.cutoff.est <- .estimateProbCutoff(df)
                  prob.cutoff.q2 <- prob.cutoff.est$probability.cutoff

                  message("Done!")

                  message("[INFO] Finalizing output...", appendLF = FALSE)

                  if (!missing(prob.cutoff)){
                          prob.cutoff.q2 <- prob.cutoff
                  }

                  pred.plus <- which(sepline.plus$pass.cov.treshold == "yes" &
                                             probs.plus$q2.prob > prob.cutoff.q2)
                  pred.minus <- which(sepline.minus$pass.cov.treshold == "yes" &
                                              probs.minus$q2.prob > prob.cutoff.q2)
                  pred.bi <- intersect(pred.plus, pred.minus)
                  pred.plus <- setdiff(pred.plus, pred.bi)
                  pred.minus <- setdiff(pred.minus, pred.bi)

                  pred.strand <- rep(".", length(peaks))
                  pred.strand[pred.plus] <- "+"
                  pred.strand[pred.minus] <- "-"
                  pred.strand[pred.bi] <- "bi"
                  pred.strand <- factor(pred.strand, levels = c("+", "-", "bi", "."))
                  message("Done!")

                  cdsObj@strandPrediction$predicted.strand <- pred.strand
                  cdsObj@strandPrediction$probability.cutoff <- prob.cutoff.est$probability.cutoff
                  cdsObj@strandPrediction$results.plus <- cbind(sepline.plus, probs.plus)
                  cdsObj@strandPrediction$results.minus <- cbind(sepline.minus, probs.minus)

                  assign(objName, cdsObj, envir = parent.frame())

          }
)


####------------------- getPeaks -------------------####
#' @rdname getPeaks-methods
setMethod("getPeaks",
          signature = signature("ChipDataSet"),
          definition = function(object){

                  if (length(object@peaks) != 0){
                          peaks <- object@peaks
                  } else {
                          stop("[ERROR] 'peaks' slot is empty.")
                  }

                  if (length(object@features) != 0){
                          features <- object@features
                          values(peaks) <- cbind(values(peaks), S4Vectors::DataFrame(features))
                  }

                  if (!is.null(object@tssOverlapPrediction$predicted.tssOverlap)){
                          predicted.tssOverlap <- object@tssOverlapPrediction$predicted.tssOverlap
                          values(peaks) <- cbind(values(peaks), S4Vectors::DataFrame(predicted.tssOverlap))
                  }

                  if (!is.null(object@strandPrediction$predicted.strand)){
                          predicted.strand <- object@strandPrediction$predicted.strand
                          values(peaks) <- cbind(values(peaks), S4Vectors::DataFrame(predicted.strand))
                  }
                  return(peaks)
          }
)


####------------------- getGenomicAnnot -------------------####
#' @rdname getGenomicAnnot-methods
setMethod("getGenomicAnnot",
          signature = signature("ChipDataSet"),
          definition = function(object){

                          return(object@genomicAnnotation)
          }
)


####------------------- getConfusionMatrix -------------------####
#' @rdname getConfusionMatrix-methods
setMethod("getConfusionMatrix",
          signature = signature("ChipDataSet"),
          definition = function(object){

                  return(object@tssOverlapPrediction$confusionMatrix)
          }
)


####------------------- getPredictorSignificance -------------------####
#' @rdname getPredictorSignificance-methods
setMethod("getPredictorSignificance",
          signature = signature("ChipDataSet"),
          definition = function(object){

                  # -1 to remove intercept pvalue
                  return(object@tssOverlapPrediction$logFitSummary[-1, 4])
          }
)


####------------------- getQuadProb -------------------####
#' @rdname getQuadProb-methods
setMethod("getQuadProb",
          signature = signature("ChipDataSet"),
          definition = function(object, strand){

                  if (strand %in% "+"){
                          return(object@strandPrediction$results.plus)
                  } else if (strand %in% "-"){
                          return(object@strandPrediction$results.minus)
                  } else {
                          stop("[ERROR] Wrong strand. Either '+' or '-'.")
                  }
          }
)


####------------------- getProbTreshold -------------------####
#' @rdname getProbTreshold-methods
setMethod("getProbTreshold",
          signature = signature("ChipDataSet"),
          definition = function(object){

                  return(object@strandPrediction$probability.cutoff)
          }
)


####------------------- plotGenomicAnnot -------------------####
if(getRversion() >= "2.15.1")  utils::globalVariables(c("tssOverlap", "value", "variable", "region", "ratio"))

#'@rdname plotGenomicAnnot-methods
setMethod("plotGenomicAnnot",
          signature = signature("ChipDataSet"),
          definition = function(object, plot.type = c("distrib", "enrich"), xlab, ylab, color){


                  if (length(object@genomicAnnotation) != 0){
                          df <- object@genomicAnnotation
                  } else {
                          stop("[ERROR] 'genomicAnnotation' slot is empty.")
                  }

                  if (missing(plot.type)){
                          plot.type <- "enrich"
                  } else {
                          plot.type <- match.arg(plot.type,  c("distrib", "enrich"))
                  }

                  if (missing(color)){
                          vals <- scale_fill_brewer(palette = "Dark2")
                  } else {
                          if (length(color) < 4 | length(color) > 4){
                                  stop("[ERROR] Please, specify 4 colors.")
                          } else {
                                  vals <- scale_fill_manual(values = color)
                          }
                  }

                  if (plot.type %in% "enrich"){
                          if (missing(xlab)){xlab <- "Genomic regions"}
                          if (missing(ylab)){ylab <- "Enrichment (observed/expected)"}

                          ggplot(data = df, aes(x = region, y = ratio, fill = region)) +
                                  geom_bar(stat = "identity") +
                                  vals +
                                  xlab(label = xlab) +
                                  ylab(label = ylab) +
                                  theme_bw() +
                                  theme(legend.position = "none",
                                        text = element_text(size = 14),
                                        axis.title.x = element_text(vjust = 0),
                                        axis.title.y = element_text(vjust = 1.5))
                  } else {
                          suppressMessages(df.mt <- reshape2::melt(df[, 1:3]))
                          if (missing(xlab)){xlab <- " "}
                          if (missing(ylab)){ylab <- "Number of peaks"}

                          ggplot(data = df.mt, aes(x = variable, y = value, fill = region)) +
                                  geom_bar(stat = "identity") +
                                  vals +
                                  xlab(xlab) +
                                  ylab(ylab) +
                                  theme_bw() +
                                  theme(text = element_text(size = 14),
                                        axis.title.x = element_text(vjust = 0),
                                        axis.title.y = element_text(vjust = 1.5))
                  }
          }
)


####------------------- plotFeatures -------------------####
#' @rdname plotFeatures-methods
setMethod("plotFeatures",
          signature = signature("ChipDataSet"),
          definition = function(object,
                                plot.type = c("box", "density"),
                                feature,
                                ncol,
                                xlab,
                                ylab,
                                color = c("#E41A1C","#377EB8"),
                                alpha = 1){

                  if (length(object@features) != 0){
                         df <- object@features
                  } else {
                          stop("[ERROR] 'features' slot is empty.")
                  }

                  # Select features if specified
                  if (!missing(feature)){
                          feature <- sapply(feature, function(x) match.arg(x, colnames(df)[-1]))
                          feature.id <- which(colnames(df) %in% feature)
                          feature.id <- unique(c(1, feature.id))
                          df <- df[, feature.id]
                  }

                  # Discart factor and character variables from the analysis
                  col.class <- unlist(lapply(df, class))
                  idx <- which(col.class %in% c("factor", "character"))[-1]
                  if (length(idx) > 0){
                          message("[INFO] ", colnames(df)[idx], " not numeric and will not be used for plotting.")
                          df <- df[, -idx]
                  }

                  if (base::ncol(df) == 0){
                          stop("[ERROR] No data to plot.")
                  }

                  # Melt data
                  suppressMessages(df.mt <- reshape2::melt(df, id.vars = colnames(df)[1]))

                  # Plot results
                  if (missing(plot.type)){
                          plot.type <- "box"
                  } else {
                          plot.type <- match.arg(plot.type, c("box", "density"))
                  }

                  if (missing(xlab)){
                          if (plot.type %in% "box") {
                                  xlab <- "Peak overlap with the TSS region\n (based on reference annotations)"
                          } else {
                                  xlab <- "Value"
                          }
                  }

                  if (missing(ylab)){
                          if (plot.type %in% "box") {
                                  ylab <- "Value"
                          } else {
                                  ylab <- "Density"
                          }
                  }

                  if (!missing(ncol)){
                          if (is.numeric(ncol)){
                                  if(ncol < 1){
                                          stop("[ERROR] Wrong number of columns.")
                                  }
                          } else {
                                  stop("[ERROR] 'ncol' needs to be numeric.")
                          }
                  } else {
                          ncol <- dim(df)[2] - 1
                  }

                  if (plot.type %in% "box"){
                          p <- ggplot(data = df.mt, aes(x = tssOverlap, y = value, fill = tssOverlap)) +
                                  geom_boxplot(alpha = alpha) +
                                  facet_wrap(~ variable, ncol = ncol, scales = "free") +
                                  scale_fill_manual(values = color) +
                                  xlab(label = xlab) +
                                  ylab(label = ylab) +
                                  theme_bw() +
                                  theme(legend.position = "none",
                                        text = element_text(size = 14),
                                        axis.title.x = element_text(vjust = -0.75),
                                        panel.margin = unit(1, "lines"))
                  }

                  if (plot.type %in% "density"){
                          p <- ggplot(data = df.mt, aes(x = value, colour = tssOverlap)) +
                                  geom_density(alpha = alpha, size = 1) +
                                  facet_wrap(~ variable , ncol = ncol, scales = "free") +
                                  scale_color_manual(values = color, name = "Peak overlap\nwith the TSS region") +
                                  xlab(label = xlab) +
                                  ylab(label = ylab) +
                                  theme_bw() +
                                  theme(text = element_text(size = 14),
                                        axis.title.x = element_text(vjust = -0.75),
                                        panel.margin = unit(1, "lines"))
                  }
                  p
          }
)


####------------------- plotROC -------------------####
#' @rdname plotROC-methods
setMethod("plotROC",
          signature = signature("ChipDataSet"),
          definition = function(object, ...){

                  if (!is.null(object@tssOverlapPrediction$roc)){
                          roc <- object@tssOverlapPrediction$roc
                  } else {
                          stop("[ERROR] 'tssOverlapPrediction' slot is empty.")
                  }

                  pROC::plot.roc(roc, ...)
                  graphics::legend("bottomright", bty = "n",
                                   legend = paste0("AUC: ", round(roc$auc, 4)),
                                   text.col = "red4", cex = 1.2)
          }
)


####------------------- peaksToBed -------------------####
#' @rdname peaksToBed-methods
setMethod("peaksToBed",
          signature = signature("ChipDataSet"),
          definition = function(object, file,
                                strand.pred.color = c("blue", "red", "green4", "black"),
                                gene.associated.peaks = TRUE){

                  if (length(object@peaks) != 0){
                          peaks <- object@peaks
                  } else {
                          stop("[ERROR] 'peaks' slot is empty.")
                  }

                  if (!is.null(object@strandPrediction$predicted.strand)){
                          peaks$predicted.strand <- object@strandPrediction$predicted.strand
                  } else {
                          stop("[ERROR] 'strandPrediction' slot is empty.")
                  }

                  if (gene.associated.peaks){
                          if (!is.null(object@tssOverlapPrediction$predicted.tssOverlap)){
                                  predicted.peaks <- object@tssOverlapPrediction$predicted.tssOverlap$predicted.tssOverlap
                                  peaks <- peaks[ which(predicted.peaks == "yes") ]
                          } else {
                                  stop("[ERROR] 'tssOverlapPrediction' slot is empty.")
                          }
                  }

                  # Transform colors to rgb
                  strand.pred.color.rgb <- sapply(strand.pred.color, grDevices::col2rgb)

                  chrom <- seqnames(peaks)
                  chromStart <- start(peaks)
                  chromEnd <- end(peaks)
                  # id <- paste0("peak_", seq_along(peaks))
                  id <- peaks$id
                  score <- rep(0, length(peaks))
                  strand <- rep(".", length(peaks))
                  thichStart <- chromStart
                  thickEnd <- chromEnd
                  itemRgb <- ifelse(test = peaks$predicted.strand == "+",
                                    paste(strand.pred.color.rgb[, 1], collapse = ","),
                                    ifelse(test = peaks$predicted.strand == "-",
                                           paste(strand.pred.color.rgb[, 2], collapse = ","),
                                           ifelse(test = peaks$predicted.strand == "bi",
                                                  paste(strand.pred.color.rgb[, 3], collapse = ","),
                                                  paste(strand.pred.color.rgb[, 4], collapse = ","))
                                    )
                  )
                  df <- data.frame(chrom, chromStart, chromEnd, id, score, strand, thichStart, thickEnd, itemRgb)
                  utils::write.table(x = df, file = file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
          }
)

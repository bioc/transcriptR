context("TranscriptionDataSet methods")

test_that("estimateBackground works", {
        data(tds)
        # No output is produced
        expect_silent(estimateBackground(object = tds, fdr.cutoff = 0.1))
        # Check the result of the function call with different fdr settings
        estimateBackground(object = tds, fdr.cutoff = 0.1)
        expect_equal(tds@coverageCutoff, 3.561)
        estimateBackground(object = tds, fdr.cutoff = 0.05)
        expect_equal(tds@coverageCutoff, 4.063)
        estimateBackground(object = tds, fdr.cutoff = 0.01)
        expect_equal(tds@coverageCutoff, 5.327)
})

test_that("estimateGapDistance works", {
        data(tds)
        annot <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
        # Check the estimated gap distance
        estimateGapDistance(object = tds, annot = annot, coverage.cutoff = 5.327,
                            filter.annot = FALSE)
        gd.test <- getTestedGapDistances(tds)[which.min(getTestedGapDistances(tds)$sum.two.errors), ]
        gd.min <- gd.test[, "gap.distance"]
        expect_is(gd.test, class = "data.frame")
        expect_equal(gd.min, 6900)

        # Check the effect of gene filtering on the resulting gap distance
        estimateGapDistance(object = tds, annot = annot, coverage.cutoff = 5.327,
                            filter.annot = TRUE, fpkm.quantile = 0.5)
        gd.test <- getTestedGapDistances(tds)[which.min(getTestedGapDistances(tds)$sum.two.errors), ]
        gd.min <- gd.test[, "gap.distance"]
        expect_equal(gd.min, 3200)
})

test_that("detectTranscripts works", {
        data(tds)
        annot <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)

        detectTranscripts(object = tds, coverage.cutoff = 5.327, gap.dist = 6900,
                          estimate.params = FALSE)
        expect_is(getTranscripts(tds), "GRanges")
        expect_equal(length(getTranscripts(tds)), 341)
        expect_equal(tds@transcriptsCovCutoff, 5.327)
        expect_equal(tds@transcriptsGapDist, 6900)
        expect_equal(tds@transcriptsNormalization, 1120345)
        exp.names <- c("id", "length")
        expect_named(object = GenomicRanges::mcols(getTranscripts(tds)), expected = exp.names)

        detectTranscripts(object = tds, coverage.cutoff = 5.327, gap.dist = 6900,
                          estimate.params = FALSE, combine.by.annot = TRUE, annot = annot)
        expect_equal(length(getTranscripts(tds)), 274)

        detectTranscripts(object = tds, coverage.cutoff = 5.327, gap.dist = 6900,
                          estimate.params = TRUE)
        exp.names <- c("id", "length", "bases.covered", "coverage", "fragments","fpkm")
        expect_named(object = GenomicRanges::mcols(getTranscripts(tds)), expected = exp.names)
        fpkm.high <- getTranscripts(tds)$fpkm

        detectTranscripts(object = tds, coverage.cutoff = 5.327, gap.dist = 6900,
                          estimate.params = TRUE, total.reads = length(tds@fragments) * 10)
        fpkm.low <- getTranscripts(tds)$fpkm
        expect_gt(sum(fpkm.high), sum(fpkm.low))
        expect_equal(tds@transcriptsNormalization, 1120345 * 10)
})

test_that("annotateTranscripts works", {
        data(tds)
        annot <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)

        detectTranscripts(object = tds, coverage.cutoff = 5.327, gap.dist = 6900,
                          estimate.params = FALSE)
        annotateTranscripts(object = tds , min.overlap = 0.5, annot = annot)
        exp.names <- c("id", "length", "annotation.overlap")
        expect_named(GenomicRanges::mcols(getTranscripts(tds)), exp.names)

        detectTranscripts(object = tds, coverage.cutoff = 5.327, gap.dist = 6900,
                          estimate.params = TRUE)
        annotateTranscripts(object = tds , min.overlap = 0.5, annot = annot)
        exp.names <- c("id", "length", "bases.covered", "coverage", "fragments","fpkm", "annotation.overlap")
        expect_named(GenomicRanges::mcols(getTranscripts(tds)), exp.names)
})

test_that("getTranscipts works", {
        data(tds)
        detectTranscripts(object = tds, coverage.cutoff = 5.327, gap.dist = 6900,
                          estimate.params = FALSE)
        expect_is(getTranscripts(tds), "GRanges")
})

test_that("getTestedGapDistances works", {
        data(tds)
        annot <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
        # Check the estimated gap distance
        estimateGapDistance(object = tds, annot = annot,
                            coverage.cutoff = 5.327, filter.annot = FALSE)
        expect_is(getTestedGapDistances(object = tds), "data.frame")
})

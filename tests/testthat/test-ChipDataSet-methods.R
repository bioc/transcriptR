context("ChipDataSet methods")

test_that("cds slots", {
        data(cds)

        # Check slots of the new object
        expect_is(cds@peaks, "GRanges")
        expect_is(cds@region, "GRanges")
        expect_is(cds@genomicAnnotation, "data.frame")
        expect_is(cds@features, "data.frame")
        expect_is(cds@tssOverlapPrediction, "list")
        expect_is(cds@strandPrediction, "list")

        expect_equal(length(cds@tssOverlapPrediction), 4)
        expect_equal(length(cds@strandPrediction), 4)

        expect_null(cds@tssOverlapPrediction$predicted.tssOverlap)
        expect_null(cds@tssOverlapPrediction$confusionMatrix)
        expect_null(cds@tssOverlapPrediction$logFitSummary)
        expect_null(cds@tssOverlapPrediction$roc)

        expect_null(cds@strandPrediction$predicted.strand)
        expect_null(cds@strandPrediction$probability.cutoff)
        expect_null(cds@strandPrediction$results.plus)
        expect_null(cds@strandPrediction$results.minus)
})

test_that("Predictions work", {
        data(cds)
        data(tds)

        # Check the slots of the object after function call
        predictTssOverlap(object = cds, feature = "pileup")
        predictStrand(cdsObj = cds, tdsObj = tds, coverage.cutoff = 5.327)

        expect_is(cds@tssOverlapPrediction$predicted.tssOverlap, "data.frame")
        expect_is(cds@tssOverlapPrediction$confusionMatrix, "confusionMatrix")
        expect_is(cds@tssOverlapPrediction$logFitSummary, "matrix")
        expect_is(cds@tssOverlapPrediction$roc, "roc")

        expect_is(cds@strandPrediction$predicted.strand, "factor")
        expect_is(cds@strandPrediction$probability.cutoff, "numeric")
        expect_is(cds@strandPrediction$results.plus, "data.frame")
        expect_is(cds@strandPrediction$results.minus, "data.frame")

        # Check the results of the strand prediction
        # expect_equal(as.numeric(table(cds@strandPrediction$predicted.strand)), c(57,44,42,240))
        # expect_equal(cds@strandPrediction$probability.cutoff, 0.64)
        exp.names <- c("max.cov", "pass.cov.treshold","q1q2.sepline.coord","q1.coord", "q2.coord",
                       "q1.count","q2.count","q1.prob","q2.prob")
        expect_named(cds@strandPrediction$results.plus, exp.names)
        expect_named(cds@strandPrediction$results.minus, exp.names)
})


test_that("addFeature works", {
        data(cds)
        addFeature(object = cds, feature = list(fake = runif(n = length(cds@peaks))))
        exp.names <- c("tssOverlap","length","fragments","density","pileup","fake")
        expect_named(cds@features, exp.names)
})

test_that("getters work", {
        data(cds)
        data(tds)

        expect_is(getPeaks(cds), "GRanges")
        expect_is(getGenomicAnnot(cds), "data.frame")
        expect_null(getConfusionMatrix(cds))
        expect_null(getPredictorSignificance(cds))
        expect_null(getQuadProb(cds, strand = "+"))
        expect_null(getQuadProb(cds, strand = "-"))
        expect_null(getProbTreshold(cds))

        predictTssOverlap(object = cds, feature = "pileup")
        predictStrand(cdsObj = cds, tdsObj = tds, coverage.cutoff = 5.327)

        expect_is(getConfusionMatrix(cds), "confusionMatrix")
        expect_is(getPredictorSignificance(cds), "numeric")
        expect_is(getQuadProb(object = cds, strand = "+"), "data.frame")
        expect_is(getQuadProb(object = cds, strand = "-"), "data.frame")
        expect_is(getProbTreshold(cds), "numeric")

        expect_error(getQuadProb(cds, strand = FALSE))
        expect_error(getQuadProb(cds, strand = 1))
        expect_error(getQuadProb(cds, strand = NA))
        expect_error(getQuadProb(cds, strand = NaN))
        expect_error(getQuadProb(cds, strand = "."))

        expect_equal(getProbTreshold(object = cds), 0.64)
})

context("Testing objects construction")

test_that("TranscriptionDataSet object construction", {
        bam.rna <- system.file("extdata", "rna_reads.bam", package = "transcriptR")
        tds <- constructTDS(file = bam.rna, fragment.size = 250, swap.strand = TRUE)
        expect_is(tds, "TranscriptionDataSet")
})

test_that("constructTDS error is expected", {
        bam.rna <- system.file("extdata", "rna_reads.bam", package = "transcriptR")
        expect_error(constructTDS(file = "file"))
        expect_error(constructTDS(file = TRUE))

        expect_error(constructTDS(file = bam.rna, fragment.size = "1"))
        expect_warning(constructTDS(file = bam.rna, fragment.size = 1000),
                       regexp = "'fragment.size' seems to be rather large")
        expect_error(constructTDS(file = bam.rna, fragment.size = FALSE))

        expect_error(constructTDS(file = bam.rna, swap.strand = "TRUE"))
        expect_error(constructTDS(file = bam.rna, swap.strand = 10))

        expect_error(constructTDS(file = bam.rna, paired.end = "TRUE"))
        expect_error(constructTDS(file = bam.rna, paired.end = 1))
})

test_that("ChipDataSet object constraction", {
        bam.chip <- system.file("extdata", "chip_reads.bam", package = "transcriptR")
        peaks <- system.file("extdata", "peaks.txt", package = "transcriptR")
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
        cds <- constructCDS(peaks = peaks, reads = bam.chip, TxDb = txdb)
        expect_is(cds, "ChipDataSet")
        expect_is(cds@features, "data.frame")
        expect_is(cds@genomicAnnotation, "data.frame")
        expect_named(cds@genomicAnnotation, c("region", "observed",
                                              "expected", "ratio"))
        expect_named(cds@features, c("tssOverlap", "length",
                                     "fragments", "density", "pileup"))
})

test_that("constructCDS error is expected", {
        bam.chip <- system.file("extdata", "chip_reads.bam", package = "transcriptR")
        peaks <- system.file("extdata", "peaks.txt", package = "transcriptR")
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

        expect_error(constructCDS(peaks = "peaks", reads = bam.chip, TxDb = txdb))
        expect_error(constructCDS(peaks = 1, reads = bam.chip, TxDb = txdb))
        expect_error(constructCDS(peaks = FALSE, reads = bam.chip, TxDb = txdb))

        expect_error(constructCDS(peaks = peaks, reads = "bam.chip", TxDb = "txdb"))
        expect_error(constructCDS(peaks = peaks, reads = 1, TxDb = FALSE))
        expect_error(constructCDS(peaks = peaks, reads = FALSE, TxDb = txdb))
})

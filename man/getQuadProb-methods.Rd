% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChipDataSet-generics.R, R/ChipDataSet-methods.R
\docType{methods}
\name{getQuadProb}
\alias{getQuadProb}
\alias{getQuadProb,ChipDataSet-method}
\title{getQuadProb}
\usage{
getQuadProb(object, strand)

\S4method{getQuadProb}{ChipDataSet}(object, strand)
}
\arguments{
\item{object}{A \code{\link{ChipDataSet}} object.}

\item{strand}{\code{Character}. The strand to extract calculations for.
One of ["+", "-"].}
}
\value{
A nine column \code{Data.Frame}, where each row corresponds to
    a ChIP-seq peak and each column keeps one of the intermediate calculations:
    \itemize{
        \item \code{max.cov} - maximum coverage of the RNA-seq fragments
            inside the peak region.

        \item \code{pass.cov.treshold} - whether the \code{max.cov} exceeds
            the \code{coverage.cutoff}, either user defined or estimated from
            RNA-seq data by \code{\link{estimateBackground}} function call and
            stored in \code{\link{TranscriptionDataSet}} object.

        \item \code{q1q2.sepline.coord} - genomic coordinate corresponding to
            the transcription start position inside the peak region.

        \item \code{q1.coord} - genomic coordinates of q1.

        \item \code{q2.coord} - genomic coordinates of q2.

        \item \code{q1.count} - total number of fragments in q1.

        \item \code{q2.count} - total number of fragments in q2.

        \item \code{q1.prob} - probability of a fragment being sampled from
            the q1.

        \item \code{q2.prob} - probability of a fragment being sampled from
           the q2.
    }
}
\description{
Retrieve all internal calculations performed by \code{\link{predictStrand}}
function.
}
\examples{
### Load ChipDataSet object
data(cds)

### Load TranscriptionDataSet object
data(tds)

head(getQuadProb(cds, strand = "+"))
head(getQuadProb(cds, strand = "-"))

}
\author{
Armen R. Karapetyan
}
\seealso{
\code{\link{predictStrand}}
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChipDataSet-generics.R, R/ChipDataSet-methods.R
\docType{methods}
\name{plotGenomicAnnot}
\alias{plotGenomicAnnot}
\alias{plotGenomicAnnot,ChipDataSet-method}
\title{plotGenomicAnnot}
\usage{
plotGenomicAnnot(object, plot.type = c("distrib", "enrich"), xlab, ylab,
  color)

\S4method{plotGenomicAnnot}{ChipDataSet}(object, plot.type = c("distrib",
  "enrich"), xlab, ylab, color)
}
\arguments{
\item{object}{A \code{\link{ChipDataSet}} object.}

\item{plot.type}{\code{Character}. One of ("distrib", "enrich").
Default: "enrich".}

\item{xlab}{\code{Character}. Title of the x-axis.}

\item{ylab}{\code{Character}. Title of the y-axis.}

\item{color}{A character vector of length four, specifying colors for
distinct genomic features (TSSs, exons, introns, intergenic regions).}
}
\value{
ggplot2 object.
}
\description{
Visualize genomic distribution of ChIP-seq peaks.
}
\details{
Genomic distribution of the peaks can be visualized in two ways, either
    by observing the total number of peaks overlapping given genomic feature
    or by looking at the enrichment levels.
}
\examples{
### Load ChipDataSet object
data(cds)

### Plot the total number of peaks overlapping distinct genomic features
plotGenomicAnnot(object = cds, plot.type = "distrib")

### Plot enrichment of the peaks at a given genomic feature (e.g. TSS)
plotGenomicAnnot(object = cds, plot.type = "enrich")

}
\author{
Armen R. Karapetyan
}
\seealso{
\code{\link{constructCDS}}
}


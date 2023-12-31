% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show.R
\docType{methods}
\name{show}
\alias{show}
\alias{show,ChipDataSet-method}
\alias{show,TranscriptionDataSet-method}
\title{Show method for ChipDataSet and TranscriptionDataSet objects}
\usage{
\S4method{show}{TranscriptionDataSet}(object)

\S4method{show}{ChipDataSet}(object)
}
\arguments{
\item{object}{A \code{\link{TranscriptionDataSet}} or
\code{\link{ChipDataSet}} object.}
}
\value{
Displays an overview of the TranscriptionDataSet or ChipDataSet object.
}
\description{
Show method for objects of class \code{\link{ChipDataSet}} and \code{\link{TranscriptionDataSet}}
}
\examples{
    ### Load TranscriptionDataSet object
    data(tds)

    ### View a short summary of the object
    tds
}


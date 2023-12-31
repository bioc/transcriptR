% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChipDataSet-generics.R, R/ChipDataSet-methods.R
\docType{methods}
\name{plotFeatures}
\alias{plotFeatures}
\alias{plotFeatures,ChipDataSet-method}
\title{plotFeatures}
\usage{
plotFeatures(object, plot.type = c("box", "density"), feature, ncol, xlab,
  ylab, color = c("#E41A1C", "#377EB8"), alpha = 1)

\S4method{plotFeatures}{ChipDataSet}(object, plot.type = c("box", "density"),
  feature, ncol, xlab, ylab, color = c("#E41A1C", "#377EB8"), alpha = 1)
}
\arguments{
\item{object}{A \code{\link{ChipDataSet}} object.}

\item{plot.type}{One of ["box", "density"]. Default: "box"}

\item{feature}{Feature to plot. By default, all the features are plotted.}

\item{ncol}{\code{Numeric}. Arrange individual plots in columns. By default,
the number of columns correspond to the number of features used for
plotting.}

\item{xlab}{\code{Character}. Title of the x-axis}

\item{ylab}{\code{Character}. Title of the y-axis.}

\item{color}{A character vector of length two. Default: ["#E41A1C","#377EB8"].}

\item{alpha}{Color transparency. In a range [0, 1]. Default: 1.}
}
\value{
ggplot2 object.
}
\description{
Visualize the relations between predictors and response variable
('tssOverlap').
}
\details{
In order to discriminate between functional or gene associated peaks and
   non-functional or background peaks, each peak in the data set is
   characterized by several features. Moreover, the user might supply her/his
   own list of features with the \code{\link{addFeature}}. Prior to fitting
   the logistic model, the relations between predictors and response variable
   (tssOverlap) can be explored with \code{plotFeatures}. Based on the plots,
   poor predictors can be excluded from the analysis to improve the model
   fit.
}
\examples{
### Load ChipDataSet object
data(cds)

### The data can be plotted in two ways
### As a boxplot
plotFeatures(object = cds, plot.type = "box")

### Or as a density plot
plotFeatures(object = cds, plot.type = "density")

### Additionally, only the subset of features can be shown
plotFeatures(object = cds, plot.type = "box", feature = c("pileup", "length"))

### The position of the graphs on the plot, can be adjusted by 'ncol' argument
plotFeatures(object = cds, plot.type = "box", ncol = 2)

}
\author{
Armen R. Karapetyan
}
\seealso{
\code{\link{constructCDS}}
}


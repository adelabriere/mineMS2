% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms2Lib.R
\docType{methods}
\name{calculateCoverage,ms2Lib-method}
\alias{calculateCoverage,ms2Lib-method}
\title{Calculate the coverage of all the patterns in the dataset.}
\usage{
\S4method{calculateCoverage}{ms2Lib}(x)
}
\arguments{
\item{x}{The ms2Lib to bo computed.}
}
\value{
The m2l object with all the coverage calculated.
}
\description{
Calculate the coverage, the total intensity covered by the patterns on the mass graphs.
 This calculation can be qute long.
}
\examples{
#Loading the data
data(m2l)

#Calculate the coverage for an ms2Lib object
m2l <- calculateCoverage(m2l)
}

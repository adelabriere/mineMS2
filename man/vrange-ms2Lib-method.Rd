% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms2Lib.R
\docType{methods}
\name{vrange,ms2Lib-method}
\alias{vrange,ms2Lib-method}
\title{Return The Range Of Iteration Of MS2Lib Object}
\usage{
\S4method{vrange}{ms2Lib}(m2l, type = c("S", "L", "P"), reduced = TRUE,
  as.number = FALSE)
}
\arguments{
\item{m2l}{AN ms2Lib object}

\item{type}{"S","L" or "P"}

\item{reduced}{Used only if "type" is set to "P", shall the filtered pattern set be returned.}

\item{as.number}{If as number is selected integer are returned without the prefix.}
}
\value{
A character vector giving the existing ids in the ms2Lib object,
}
\description{
Return the full range of iteration for different objects for an MS2lib object.
}
\examples{
#Loading the data
data(m2l)

#Range of iterations for spectra
vrange(m2l,"S")
}

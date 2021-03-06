% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexPlot.R
\name{plotPatterns}
\alias{plotPatterns}
\title{Plot patterns and occurences of an ms2Lib object}
\usage{
plotPatterns(m2l, ids = NULL, components = NULL, occurences = TRUE,
  full = FALSE, byPage = 9, titles = NULL)
}
\arguments{
\item{m2l}{An ms2Lib object}

\item{ids}{The ids to be plotted or NULL if all the ids are supposed to be plot}

\item{components}{An id giving the ocmponent of each spectrum to be plotted in first page.}

\item{occurences}{Shall include the occurences be plotted aswell as the graph.}

\item{full}{Shall the full patterns be plotted (it can take some times)}

\item{byPage}{Maximum number of occurences by page.}

\item{titles}{A vector giving the titles of the MS-MS spectra.}
}
\value{
Nothing
}
\description{
Plot patterns graph and occurences of an ms2Lib object.
}
\examples{
#Loading the data
data(m2l)

#plotting the patterns
plotPatterns(m2l,c("P30","P51"))
}

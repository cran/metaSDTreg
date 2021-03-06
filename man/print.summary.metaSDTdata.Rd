% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaSDTdata.r
\name{print.summary.metaSDTdata}
\alias{print.summary.metaSDTdata}
\title{Print method for a summary.metaSDTdata object}
\usage{
\method{print}{summary.metaSDTdata}(x, ...)
}
\arguments{
\item{x}{A metaSDTdata object.}

\item{...}{For future methods.}
}
\value{
Invisible.
}
\description{
Print method for a summary.metaSDTdata object
}
\examples{
metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
summary(metadata)

}

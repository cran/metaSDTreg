% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaSDTreg.r
\name{print.metaSDTreg}
\alias{print.metaSDTreg}
\title{Print method for metaSDTreg}
\usage{
\method{print}{metaSDTreg}(x, ...)
}
\arguments{
\item{x}{An object of class metaSDTreg.}

\item{...}{For future methods.}
}
\value{
Invisible.
}
\description{
Print method for metaSDTreg
}
\examples{
## Declare simulated data as metaSDTdata
metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')

## Fit model to subset of data
(fit <- metaSDTreg(A ~ signal,
            data=metadata,
            subset = m <= 20))

}

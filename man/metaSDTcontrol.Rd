% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaSDTreg.r
\name{metaSDTcontrol}
\alias{metaSDTcontrol}
\title{Control for metaSDTreg}
\usage{
metaSDTcontrol(
  meth = "BFGS",
  start.vals = NULL,
  hessianMeth = TRUE,
  returnMaxLik = FALSE,
  loglikUse = c("R", "C++"),
  scoreFun = FALSE,
  c.maxLik = list(),
  c.clm = list(maxIter = 1000, maxModIter = 20),
  ...
)
}
\arguments{
\item{meth}{Method used by maxLik to perform maximisation, which must support constrained optimisation, see \code{\link[maxLik]{maxLik}}. Defaults to 'BFGS'.}

\item{start.vals}{Starting values for the maximisation algorithm. Any NA entries will be set to zero. Defaults to NULL, in which case starting values will be obtained from the partial proportional odds model as implemented in the package \pkg{ordinal} using \code{\link[ordinal]{clm}}.}

\item{hessianMeth}{Method to compute the Hessian of the optimisation. Passed to the argument finalHessian in \code{\link[maxLik]{maxLik}}.}

\item{returnMaxLik}{Should the function return the fit object from \code{\link[maxLik]{maxLik}} rather that the \code{\link{metaSDTreg}} object? Defaults to FALSE.}

\item{loglikUse}{For testing purposes. The C++ option is not presently included. Should the procedure use the log-likelihood implemented in R (default), or that in C++? Option will be dropped in future versions.}

\item{scoreFun}{For testing purposes. Should the optimisation procedure use the analytic score function? If FALSE (and 'meth' relies on gradients) numerical methods are used. Presently, only numerical methods are implemented. Defaults to FALSE. Option default may change in future versions.}

\item{c.maxLik}{A list that is passed to the 'MaxControl' object in \code{\link[maxLik]{maxLik}}. The list should contain control arguments to be passed to the maximisation method used by maxLik (the default being maxBFGS). In particular, 'iterlim' is the number of iterations (default 200) and a larger number of 'printlevel' (default 0) prints more information from the maximisation procedure.}

\item{c.clm}{A list that is passed to clm.control containing controls to be used by \code{\link[ordinal]{clm}}. The clm is used to determine starting values. By default, maxIter is set to 1000 and maxModIter is set to 20.}

\item{...}{For future methods.}
}
\value{
A list of class 'metaSDTcontrol' containing the control arguments.
}
\description{
Control function for metaSDTreg.
The function allows passing arguments to the maxLik and clm procedure used for optimisation and starting values, respectively.
}
\examples{
names(metaSDTcontrol())

}

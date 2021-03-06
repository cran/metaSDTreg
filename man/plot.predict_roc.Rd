% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_roc.r
\name{plot.predict_roc}
\alias{plot.predict_roc}
\title{Plot predicted ROC curve}
\usage{
\method{plot}{predict_roc}(x, thr = seq(-10, 10, length = 1000), rline = TRUE, ...)
}
\arguments{
\item{x}{A 'predict_roc' object to plot.}

\item{thr}{The sequence of thresholds parametrising the ROC curve, if this is a function. Default to a length 1000 sequence from -10 to 10.}

\item{rline}{Logical, should the line of random discrimination be added to the plot? Defaults to TRUE.}

\item{...}{Addtional arguments passed to \code{plot}.}
}
\value{
Invisible.
}
\description{
Plot method for objects of class 'predict_roc'.
}
\details{
If neither 'xlab' nor 'ylab' is passed to \code{plot} the function supplies default x- and y-axis labels based on the type of ROC curve.
}
\examples{
## Declare simulated data as metaSDTdata
metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')

## Fit model to subset of data
fit <- metaSDTreg(A ~ signal,
            data=metadata,
            subset = m <= 20)

## Model-predicted signal-specific ROC curve
signalROC <- predict_roc(fit, type = 's')
plot(signalROC, type = 'l')

}

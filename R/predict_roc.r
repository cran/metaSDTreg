##' Generic predict_roc method
##'
##' Predict ROC curve from signal detection theory model.
##'
##' @param object Signal detection theory model.
##' @param ... Further arguments passed to methods.
##'
##' @return An object of class 'predict_roc'.
##'
##' @seealso \code{\link{predict_roc.metaSDTreg}}, \code{\link{predict_roc.metaSDTdata}}.
##' 
##' @export
##' 
predict_roc <- function(object, ...) {
    UseMethod("predict_roc")
}


##' Predicted ROC curve
##'
##' Predict ROC curves from metaSDTreg object.
##'
##' The 'metaSDTreg' object given to the function must have named coefficients with names as they would be if \code{metaSDTreg} is run without user-supplied starting values.
##'
##' A ROC curve is a 2-D curve parametrised by some x given by c(FA(x), HR(x)) where FA is the false alarm rate and HR is the hit rate. For example, for type 1 ROC,
##' \deqn{FA(x) = 1 - pnorm(x - s0*d),} 
##' \deqn{HR(x) = 1 - pnorm(x - s1*d),}
##' where \eqn{d} is the signal sensitivity.
##' 
##' Note that the predicted ROC curve is for a reference individual in the regression, i.e. additional covariates are not entered into the ROC so that reparametrisation of the 'metaSDTreg' model is needed to change predictions.  
##'
##' @param object An object of class \code{metaSDTreg}.
##' @param type The type of ROC curve to predict. A character string, where '1' requests the type 1 ROC curve (the default), 'n' requests the type 2 noise-specific and 's' the type 2 signal-specific ROC curve.
##' @param s0 Numeric, the value of 'signal' to regard as 'noise'. Defaults to 0.
##' @param s1 Numeric, the value of 'signal' to regard as 'signal'. Defaults to 1.
##' @param ... For future methods
##' 
##' @return A function of class 'predict_roc' containing the appropriate ROC curve. This is a function of x which returns c(FA,HR), where FA is the false alarm rate and HR is the hit rate.
##'
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Fit model to subset of data
##' fit <- metaSDTreg(A ~ signal,
##'             data=metadata,
##'             subset = m <= 20)
##'
##' ## Model-predicted signal-specific ROC curve
##' signalROC <- predict_roc(fit, type = 's')
##' 
##' @references Maniscalco, B., & Lau, H. (2014).
##' Signal Detection Theory Analysis of Type 1 and Type 2 Data: Meta-d , Response-Specific Meta-d , and the Unequal Variance SDT Model.
##' In S. M. Fleming, & C. D. Frith (Eds.), The {Cognitive} {Neuroscience} of {Metacognition}
##' (pp. 25 66). : Springer Berlin Heidelberg.
##'
##' @export
##' 
predict_roc.metaSDTreg <- function(object, type = c("1", "n", "s"), s0 = 0, s1 = 1, 
    ...) {
    ## Checks
    if (!inherits(object, "metaSDTreg")) {
        stop("Object not of class 'metaSDTreg'.")
    }
    if (length(intersect(c("d.n", "d", "d.s", "theta"), names(stats::coef(object)))) == 
        0) {
        stop("Coefficients in object have non-standard names.")
    }
    
    t0 <- type[1]
    coefs <- stats::coef(object)
    out <- list()
    
    ## Calculate appropriate ROC; type 1
    if (t0 == "1") {
        out$type1ROC <- function(x) {
            c(FA = 1 - stats::pnorm(x, mean = s0 * coefs["d"]), HR = 1 - stats::pnorm(x, 
                mean = s1 * coefs["d"]))
        }
    } else {
        if (t0 == "n") {
            # type 2, noise
            out$type2ROCn <- function(x) {
                c(FA = truncnorm::ptruncnorm(x, mean = s1 * coefs["d.n"], b = coefs["theta"]), 
                  HR = truncnorm::ptruncnorm(x, mean = s0 * coefs["d.n"], b = coefs["theta"]))
            }
        } else {
            if (t0 == "s") {
                # type 2, signal
                out$type2ROCs <- function(x) {
                  c(FA = 1 - truncnorm::ptruncnorm(x, mean = s0 * coefs["d.s"], a = coefs["theta"]), 
                    HR = 1 - truncnorm::ptruncnorm(x, mean = s1 * coefs["d.s"], a = coefs["theta"]))
                }
            } else {
                stop("Invalid type")
            }
        }
    }
    
    out <- c(out, type = t0)
    class(out) <- "predict_roc"
    
    return(out)
}


##' Observed ROC points
##'
##' The observed points of the ROC curve from a 'metaSDTdata' object.
##'
##' Note that the type 1 ROC points arise by using each criterion in turn to decide between 'signal' and 'noise'. Since this involves also the type 2 thresholds, such a curve is also sometimes referred to as a 'pseudo' ROC curve.
##'
##' @param object A 'metaSDTdata' object from which to calculate observed ROC points.
##' @param type The type of ROC curve to predict. A character string, where '1' requests the type 1 ROC curve, 'n' requests the type 2 noise-specific and 's' the type 2 signal-specific ROC curve.
##' @param s0 Numeric, the value of object$signal to regard as 'noise'. Defaults to 0.
##' @param s1 Numeric, the value of object$signal to regard as 'signal'. Defaults to 1.
##' @param ... For future methods
##' 
##' @return A matrix two-column matrix of class 'predict_roc' with one row of c(FA, HR) per threshold (FA: False Alarm rate, HR: Hit Rate).
##' 
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Observed signal-specific ROC curve
##' signalROC <- predict_roc(metadata, type = 's')
##' 
##' @references Maniscalco, B., & Lau, H. (2014).
##' Signal Detection Theory Analysis of Type 1 and Type 2 Data: Meta-d , Response-Specific Meta-d , and the Unequal Variance SDT Model.
##' In S. M. Fleming, & C. D. Frith (Eds.), The {Cognitive} {Neuroscience} of {Metacognition}
##' (pp. 25 66). : Springer Berlin Heidelberg.
##' 
##' @export
##'
predict_roc.metaSDTdata <- function(object, type = c("1", "n", "s"), s0 = 0, s1 = 1, 
    ...) {
    ## Checks
    if (!inherits(object, "metaSDTdata")) {
        stop("Object not of class 'metaSDTdata'.")
    }
    
    t0 <- type[1]
    df <- object[object$signal %in% c(s0, s1), , drop = FALSE]
    out <- list()
    
    ## type 1 -- dummy matrix of indicators
    if (t0 == "1") {
        d.mat <- stats::model.matrix(~A - 1 + signal, df)
        ## -- cumulative probabilities
        d.mat[, 1:(ncol(d.mat) - 1)] <- t(apply(d.mat[, 1:(ncol(d.mat) - 1)], MARGIN = 1, 
            FUN = cumsum))
        
        cp.mat <- apply(d.mat[, 1:(ncol(d.mat) - 2)], MARGIN = 2, FUN = function(x) {
            prop.table(table(x, d.mat[, ncol(d.mat)]), margin = 2)[2, ]
        })
        
        ## -- matrix of survival probabilities
        p.mat <- t(1 - cp.mat)
        dimnames(p.mat) <- list(c(), c("FA", "HR"))
        
        out$type1ROC <- p.mat
    } else {
        if (t0 == "n") {
            # type 2, noise
            df2 <- df[df$A <= attr(object, "L"), , drop = FALSE]
            df2$A <- factor(df2$A, levels = levels(df2$A)[1:attr(object, "L")])
            
            d.mat <- stats::model.matrix(~A - 1 + signal, df2)
            
            d.mat[, 1:(ncol(d.mat) - 1)] <- t(apply(d.mat[, 1:(ncol(d.mat) - 1)], 
                MARGIN = 1, FUN = cumsum))
            
            cp.mat <- apply(d.mat[, 1:(ncol(d.mat) - 2)], MARGIN = 2, FUN = function(x) {
                prop.table(table(x, d.mat[, "signal"]), margin = 2)[2, ]
            })
            
            p.mat <- t(1 - cp.mat)
            dimnames(p.mat) <- list(c(), c("FA2n", "HR2n"))
            
            out$type2ROCn <- p.mat
        } else {
            if (t0 == "s") {
                # type 2, signal
                df2 <- df[df$A > attr(object, "L"), , drop = FALSE]
                df2$A <- factor(df2$A, levels = levels(df2$A)[(attr(object, "L") + 
                  1):nlevels(df2$A)])
                
                d.mat <- stats::model.matrix(~A - 1 + signal, df2)
                
                d.mat[, 1:(ncol(d.mat) - 1)] <- t(apply(d.mat[, 1:(ncol(d.mat) - 
                  1)], MARGIN = 1, FUN = cumsum))
                
                cp.mat <- apply(d.mat[, 1:(ncol(d.mat) - 2)], MARGIN = 2, FUN = function(x) {
                  prop.table(table(x, d.mat[, "signal"]), margin = 2)[2, ]
                })
                
                p.mat <- t(1 - cp.mat)
                dimnames(p.mat) <- list(c(), c("FA2s", "HR2s"))
                
                out$type2ROCs <- p.mat
            } else {
                stop("Invalid type")
            }
        }
    }
    
    out <- c(out, type = t0)
    class(out) <- "predict_roc"
    
    return(out)
}


##' Extract coordinates for predicted ROC curve
##'
##' @param x A 'predict_roc' object to plot.
##' @param thr The sequence of thresholds parametrising the ROC curve, if this is a function.
##'
##' @keywords internal
##' 
ROCcoords <- function(x, thr) {
    if (inherits(x[[1]], "function")) {
        pp <- vapply(thr,
                     FUN = x[[1]],
                     FUN.VALUE = array(double(), dim = c(2, 1)))
        pp <- t(pp[, 1, ])  # simplify to 2-dim array
    } else {
        pp <- x[[1]]
    }
    
    return(pp)
}



##' Plot predicted ROC curve
##'
##' Plot method for objects of class 'predict_roc'. 
##' 
##' If neither 'xlab' nor 'ylab' is passed to \code{plot} the function supplies default x- and y-axis labels based on the type of ROC curve.
# ##' When plotting a predicted ROC curve, the points (0,0) and (1,1) are always
# plotted.
##' 
##' @param x A 'predict_roc' object to plot.
##' @param thr The sequence of thresholds parametrising the ROC curve, if this is a function. Default to a length 1000 sequence from -10 to 10.
##' @param rline Logical, should the line of random discrimination be added to the plot? Defaults to TRUE.
##' @param ... Addtional arguments passed to \code{plot}.
##'
##' @return Invisible.
##'
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Fit model to subset of data
##' fit <- metaSDTreg(A ~ signal,
##'             data=metadata,
##'             subset = m <= 20)
##'
##' ## Model-predicted signal-specific ROC curve
##' signalROC <- predict_roc(fit, type = 's')
##' plot(signalROC, type = 'l')
##'
##' @export
##' 
plot.predict_roc <- function(x, thr = seq(-10, 10, length = 1000), rline = TRUE, 
    ...) {
    pp <- ROCcoords(x, thr)
    
    dots <- list(...)
    pargs <- c(list(formula = pp[, 2] ~ pp[, 1]), dots)
    
    if (!("xlab" %in% names(dots)) & !("ylab" %in% names(dots)) == TRUE) {
        lab0 <- c("Type 1 ", "Type 2 noise-specific ", "Type 2 signal-specific ")
        lab1 <- paste0(lab0[match(x[["type"]], c("1", "n", "s"))], c("FA", "HR"))
        pargs$xlab <- lab1[1]
        pargs$ylab <- lab1[2]
    }
    
    do.call(graphics::plot, args = pargs)
    if (rline) {
        graphics::abline(0, 1, lty = 2)
    }
    
    invisible()
}

##' Lines of predicted ROC curve
##'
##' Plot lines method for objects of class 'predict_roc'. 
##' 
##' @param x A 'predict_roc' object to plot.
##' @param thr The sequence of thresholds parametrising the ROC curve, if this is a function. Default to a length 1000 sequence from -10 to 10.
##' @param ... Addtional arguments passed to \code{lines}.
##' 
##' @return Invisible.
##'
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Fit model to subset of data
##' fit <- metaSDTreg(A ~ signal,
##'             data=metadata,
##'             subset = m <= 20)
##'
##' ## Plot observed type 1 ROC points
##' plot(predict_roc(metadata, type = '1'), xlim = 0:1, ylim = 0:1)
##'
##' ## Add Model-predicted ROC curve (estimated from subset of data)
##' lines(predict_roc(fit, type = '1'))
##' 
##' @export
##' 
lines.predict_roc <- function(x, thr = seq(-10, 10, length = 1000), ...) {
    pp <- ROCcoords(x, thr)
    
    dots <- list(...)
    pargs <- c(list(x = pp[, 1], y = pp[, 2]), dots)
    
    do.call(graphics::lines, args = pargs)
    invisible()
}

##' Points from predicted ROC curve
##'
##' Plot points method for objects of class 'predict_roc'. 
##' 
##' @param x A 'predict_roc' object to plot.
##' @param thr The sequence of thresholds parametrising the ROC curve, if this is a function. Default to a length 1000 sequence from -10 to 10.
##' @param ... Addtional arguments passed to \code{lines}.
##' 
##' @return Invisible.
##'
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Fit model to subset of data
##' fit <- metaSDTreg(A ~ signal,
##'             data=metadata,
##'             subset = m <= 20)
##'
##' ## Plot observed type 1 ROC points
##' plot(predict_roc(metadata, type = '1'), xlim = 0:1, ylim = 0:1, pch = 'x')
##'
##' ## Add Model-predicted ROC curve (estimated from subset of data)
##' points(predict_roc(fit, type = '1'))
##'
##' @export
##' 
points.predict_roc <- function(x, thr = seq(-10, 10, length = 1000), ...) {
    pp <- ROCcoords(x, thr)
    
    dots <- list(...)
    pargs <- c(list(x = pp[, 1], y = pp[, 2]), dots)
    
    do.call(graphics::points, args = pargs)
    invisible()
}

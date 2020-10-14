##' Simulated metacognition experiment
##'
##' Data is simulated from the meta-SDT model, formula (3)-(4) in Kristensen, Sandberg and Bibby (2020).
##' The true values of the simulation parameters are
##' c(tau.n.2=-0.5, tau.n.1=-0.2, theta=0,
##' tau.s.1=0.8, tau.s.2=1, d.n=0.5, d=1, d.s=0.75).
##'
##' @format A data frame with 25000 obs. of  5 variables:
##' \describe{
##'   \item{id}{Subject id.}
##'   \item{m}{Observation id, within subject.}
##'   \item{S}{Signal indicator (0 or 1).}
##'   \item{resp}{Type 1 response (0 or 1).}
##'   \item{conf}{Type 2 response, ordered (1, 2 or 3).} 
##' }
##'
##' @references
##' Kristensen, S. B., Sandberg, K., & Bibby, B. M. (2020). Regression
##' methods for metacognitive sensitivity. Journal of Mathematical
##' Psychology, 94. <doi:10.1016/j.jmp.2019.102297>.
##' 
"simMetaData"


##' Meta-SDT log-likelihood
##' 
##' Meta SDT log-likelihood function, internal function. Used by \code{metaSDTreg}.
##' 
##' @param parm Parameters of the model, c(intercepts^T,d.n,d,d.s,beta^T)^T, where length(intercepts)=K0-1 and lenght(beta)=ncol(X).
##' @param df Data frame
##' @param A Matrix of responses
##' @param X Matrix of covariates
##' @param X.n Presently not in use
##' @param X.s Presently not in use
##' @param L0 Levels of the type 2 task
##' @param K0 Levels of the ordinal variable A
##'
##' @keywords internal
##' 
metaSDT.loglik <- function(parm, df, A, X, X.n = NULL, X.s = NULL, L0, K0) {
    ## Common mean structure
    if (ncol(X) > 0) {
        common <- X %*% parm[(K0 + 3):length(parm)]
    } else {
        common <- 0
    }
    
    ## Cumulative probabilities
    p1L <- vapply(c(parm[1:(L0 - 1)], Inf), FUN = truncnorm::ptruncnorm, FUN.VALUE = vector("numeric", 
        length = length(df$signal)), b = parm[L0], mean = parm[K0] * df$signal + 
        common)
    
    pL <- stats::pnorm(parm[L0], mean = parm[K0 + 1] * df$signal + common)
    
    pLK <- vapply(c(parm[(L0 + 1):(K0 - 1)], Inf), FUN = truncnorm::ptruncnorm, FUN.VALUE = vector("numeric", 
        length = length(df$signal)), a = parm[L0], mean = parm[K0 + 2] * df$signal + 
        common)
    
    ## Combine into P
    p.mat <- Matrix::Matrix(cbind(pL, 1 - pL, p1L - cbind(0, p1L[, 1:(L0 - 1)]), 
        pLK - cbind(0, pLK[, 1:(L0 - 1)])), sparse = FALSE)
    
    ## Output log-likelihood
    return(sum(A * log(p.mat)))
}


##' Starting values from PPO
##'
##' Used by \command{metaSDTreg} to calculate starting values using the partial proportional odds model as implemented in \pkg{ordinal}. Internal function.
##'
##' The function calls the \code{\link[ordinal]{clm}} function using \code{all.vars(update(formula0, . ~ . - signal))} as the formula argument and using \code{~ signal} as the nominal argument (also see \code{\link{metaSDTreg}} for the meaning of the signal variable).
##' The coefficients of the signal variable are weighted together within the noise- and signal-specific type 2 responses using their inverse standard errors as weights.
##' 
##' @param formula0 See \code{\link{metaSDTreg}}.
##' @param data0 Data frame.
##' @param control0 See \code{\link{metaSDTcontrol}}.
##' @param L0 Number of levels in ordinal type 2 response.
##' @param K0 Number of levels of ordinal response 'A'.
##' @param n.covar0 Number of covariates aside from the 'signal' variable.
##'
##' @return A named vector of starting values.
##'
##' @keywords internal
##' 
starting.vals.PPO <- function(formula0, data0, control0, L0, K0, n.covar0) {
    
    s.mod <- ordinal::clm(stats::update(formula0, A ~ . - signal), nominal = ~signal, 
        link = "probit", data = data0, control = control0$c.clm)
    
    ## Weigh d estimates together using inverse SE as weights
    b.s <- stats::coef(s.mod)
    se.s <- sqrt(diag(stats::vcov(s.mod)))
    w.s <- 1/se.s
    
    if (n.covar0 > 0) {
        s.vals0 <- c(b.s[1:(K0 - 1)], sum(b.s[K0:(K0 + L0 - 2)] * w.s[K0:(K0 + L0 - 
            2)])/sum(w.s[K0:(K0 + L0 - 2)]), b.s[K0 + L0 - 1], sum(b.s[(K0 + L0):(2 * 
            K0 - 2)] * w.s[(K0 + L0):(2 * K0 - 2)])/sum(w.s[(K0 + L0):(2 * K0 - 2)]), 
            b.s[(2 * K0 - 1):(2 * K0 - 1 + n.covar0 - 1)])
    } else {
        s.vals0 <- c(b.s[1:(K0 - 1)], sum(b.s[K0:(K0 + L0 - 2)] * w.s[K0:(K0 + L0 - 
            2)])/sum(w.s[K0:(K0 + L0 - 2)]), b.s[K0 + L0 - 1], sum(b.s[(K0 + L0):(2 * 
            K0 - 2)] * w.s[(K0 + L0):(2 * K0 - 2)])/sum(w.s[(K0 + L0):(2 * K0 - 2)]))
    }
    
    ## Change slope signs and names
    s.vals <- s.vals0 * c(rep(1, K0 - 1), rep(-1, length(s.vals0) - K0 + 1))
    
    if (n.covar0 > 0) {
        names(s.vals) <- c(paste0("tau.n.", (L0 - 1):1), "theta", paste0("tau.s.", 
            1:(L0 - 1)), "d.n", "d", "d.s", names(s.vals)[(K0 + 3):length(names(s.vals))])
    } else {
        names(s.vals) <- c(paste0("tau.n.", (L0 - 1):1), "theta", paste0("tau.s.", 
            1:(L0 - 1)), "d.n", "d", "d.s")
    }
    
    ## Output
    return(s.vals)
    
}


##' Meta-SDT regression
##' 
##' Fit the meta SDT regression model.
##'
##' The input data frame should be of class 'metaSDTdata' as constructed by the function \code{\link{metaSDTdata}}. This will contain a variable, 'signal' which presently must be included on the right-hand side of the formula. The left hand side must be the ordinal variable 'A' also contained in the metaSDTdata object.
##' The 'signal' variable has a special interpretation as its coefficients are the signal sensitivities (type 1 d-prime and the two response specific type 2 sensitivities meta-d-prime). Presently all other variables, as determined by \code{all.vars(update(formula, . ~ . - signal ))}, enter proportionally in the model, i.e. they share cofficients across response scale levels. See Kristensen & al. (2020) for details.
##' This can be written formally as follows (in the notation of Kristensen & al. (2020)): 
##'     \deqn{    V_{\mathcal{N}} = V_{\mathcal{S}} = V ,  }
##' and \eqn{ \beta_{\mathcal{N}} }, \eqn{ \beta_{\mathcal{S}} } and \eqn{ \beta } agree on all but the first entrance which is the signal sensitivities \eqn{ d_{mathcal{N}} }, \eqn{ d_{mathcal{S}} } and \eqn{ d }, respectively.
##' Note that 'formula' specifies the mean of the latent variables and not the threshold model. Accordingly, attemps to remove the intercept in 'formula' will be ignored with a warning.
##' The function fails when the left-hand side of 'formula' is not the ordinal variable 'A' and when 'data' is not a 'metaSDTdata' object. Future versions may be less defensive and simply issue a warning in these cases.
##' Note that constrained optimisation is used to maximise the likelihood under ordinality of the thresholds. The variance-covariance matrix is not adjusted following the constrained estimation, cf. documentation of \code{\link[maxLik]{maxLik}}.
##' When interpreting results from \code{summary.metaSDTreg()} care should be taken regarding the z-testor and associated p-value for threshold parameters, as the distribution of the statistic depends on the null hypothesis which will only be reasonable under special circumstances as in other types of ordinal models.  
##' 
##' @param formula Formula specifying the regression model. Presently, the left-hand side should be the ordinal variable A, while the right-hand side must contain the signal variable from the metaSDTdata object. Note that the variable 'signal' has a special interpretation in the model, see 'Details' below.
##' @param data Data frame to fit the model on. Should be declared as metaSDTdata using the \code{\link{metaSDTdata}} function.
##' @param subset Optional argument specifying a subset of the data to be used. A logical statement, which is evaluated in 'data'. Only rows with TRUE are used in the analysis.
##' @param na.action Method for handling missing values, see na.action. Default is na.fail which stops the function with an error in the presence of missing values among variables entering 'formula'.
##' @param control Optimisation control, see \code{\link{metaSDTcontrol}}.
##'
##'
##' @return An object of class 'metaSDTreg'. This is a list object containing,
##' \itemize{
##'  \item logLik: The log likelihood after optimisation.
##'  \item coefficients: Estimated coefficients.
##'  \item vcov: Variance-covariance matrix of the maxLik object.
##'  \item call: The call issued to metaSDTreg.
##'  \item na.act: The NA action.
##' }
##' @references
##' Kristensen, S. B., Sandberg, K., & Bibby, B. M. (2020). Regression
##' methods for metacognitive sensitivity. Journal of Mathematical
##' Psychology, 94. <doi:10.1016/j.jmp.2019.102297>.
##'
##' @import maxLik
##' @import ordinal
##' @import Matrix
##' @importFrom stats na.fail
##'
##'
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##' 
##' ## Fit function to data of first 20 replications per subject
##' fit_sub <- metaSDTreg(A ~ signal,
##'                       data=metadata,
##'                       subset = m <= 20)
##' summary(fit_sub)
##'
##' \donttest{
##' ## Fit model to estimate thresholds and sensitivities
##' fit <- metaSDTreg(A ~ signal,
##'             data=metadata)
##'
##' ## True values are
##' ## c(tau.n.2=-0.5, tau.n.1=-0.2, theta=0,
##' ##   tau.s.1=0.8, tau.s.2=1, d.n=0.5, d=1, d.s=0.75)
##' coef(fit)
##' }
##'
##' 
##' @export
##' 
metaSDTreg <- function(formula, data, subset, na.action = na.fail, control = metaSDTcontrol()) {
    ## Checks
    if (!inherits(data, "metaSDTdata")) {
        stop("Data frame not declared as meta SDT data.")
    }
    if (is.na(match("L", names(attributes(data))))) {
        stop("Data frame does not have attributes matching those of meta SDT data.")
    }
    if (!inherits(formula, "formula")) {
        stop("'formula' is not a formula.")
    }
    if (all.vars(stats::update(formula, . ~ 0)) != "A") {
        stop("Presently, the outcome must be the ordinal variable A in the meta SDT data.")
    }
    if (!("signal" %in% all.vars(formula))) {
        stop("The model is not meaningful without the signal variable from the meta SDT data.")
    }
    if (is.factor(data[["signal"]])) {
        stop("The model is not meaningful for factor variable signal.")
    }
    if (attr(stats::terms(formula), "intercept") == 0) {
        warning("Attempt to remove intercept will be ignored.")
    }
    if (!inherits(control, "metaSDTcontrol")) {
        stop("Control not of class 'metaSDTcontrol'.")
    }
    
    L <- attr(data, "L")
    K <- 2 * L
    
    ## Subsetting and NA values
    if (missing(subset)) {
        df0 <- data
    } else {
        subset_call <- substitute(subset)
        df0 <- data[eval(subset_call, envir = data, enclos = parent.frame()), , drop = FALSE]
    }
    df0 <- df0[all.vars(formula)]
    data <- do.call(na.action, args = list(df0))
    
    ## Build design matrices X.mat = Matrix::Matrix(
    ## stats::model.matrix(stats::update(formula,.~.-signal), data=data) )
    X.mat <- stats::model.matrix(stats::update(formula, . ~ . - signal), data = data)
    X.mat <- X.mat[, -1, drop = FALSE]  # remove '(Intercept)'
    
    ## X.s.mat= #Presently not used X.n.mat=
    
    n.covar <- ncol(X.mat)
    
    ## Matrix with dummy variables for responses
    data$a <- factor(data$A, ordered = TRUE)
    A.mat <- Matrix::Matrix(stats::model.matrix(~a - 1, data), sparse = TRUE)
    A.mat <- cbind(as.numeric(data$a <= L), as.numeric(data$a > L), A.mat)
    
    
    ## Code constraint matrices for the intercepts (K-2 constraints)
    constr.mats <- list(ineqA = cbind(diag(-1, nrow = K - 2, ncol = K - 1) + diag(1, 
        nrow = K - 1, ncol = K - 1)[-1, ], array(0, dim = c(K - 2, 3)), array(0, 
        dim = c(K - 2, n.covar))), ineqB = rep(0, K - 2))
    
    ## Starting values
    if (is.null(control$start.vals)) {
        tryCatch({
            s.vals <- starting.vals.PPO(formula0 = formula, data0 = data, control0 = control, 
                L0 = L, K0 = K, n.covar0 = n.covar)  # PPO model
        }, error = function(e) {
            stop("Starting values: ", conditionMessage(e), "\n")
        }, warning = function(w) {
            warning(paste0("Starting values: ", conditionMessage(w), "\n"))
        })
    } else {
        s.vals <- control$start.vals
    }
    
    if (!exists("s.vals")) {
        stop("Starting values could not be determined.")
    }
    
    if (sum(is.na(s.vals)) > 0) {
        warning("Some starting values are NA, set to start at zero.")
    }
    s.vals[is.na(s.vals)] <- 0
    
    ## Presently not implemented
    if (control$scoreFun) {
        scorefun <- NULL  # metaSDT.scorefun
    } else {
        scorefun <- NULL
    }
    
    
    ## Maximise log likelihood
    if (control$loglikUse == "C++") {
        warning(paste("C++ implementation of Log-likelihood not implemented in present version.", 
            "Using R implementation."))
        metaSDTloglik <- metaSDT.loglik
        fit <- maxLik::maxLik(metaSDTloglik, grad = scorefun, signal = data$signal, 
            A = A.mat, X = X.mat, L = L, start = s.vals, method = control$meth, finalHessian = control$hessianMeth, 
            constraints = constr.mats, control = control$c.maxlik)
    } else {
        if (control$loglikUse == "R") {
            fit <- maxLik::maxLik(logLik = metaSDT.loglik, grad = scorefun, df = data, 
                A = A.mat, X = X.mat, L0 = L, K0 = K, start = s.vals, method = control$meth, 
                finalHessian = control$hessianMeth, constraints = constr.mats, control = control$c.maxlik)
        } else {
            stop("Invalid log-likelihood type: ", control$loglikUse)
        }
    }
    
    
    ## Check if maxLik object requested
    if (control$returnMaxLik) {
        return(fit)
    }
    
    ## Build output
    out <- list()
    out$logLik <- fit$maximum
    out$coefficients <- fit$estimate
    out$vcov <- stats::vcov(fit)
    out$call <- match.call()
    out$na.act <- na.action
    
    class(out) <- "metaSDTreg"
    
    ## Output
    return(out)
}


##' Control for metaSDTreg
##'
##' 
##' Control function for metaSDTreg.
##' The function allows passing arguments to the maxLik and clm procedure used for optimisation and starting values, respectively.
##' 
##' @param meth Method used by maxLik to perform maximisation, which must support constrained optimisation, see \code{\link[maxLik]{maxLik}}. Defaults to 'BFGS'.
##' @param start.vals Starting values for the maximisation algorithm. Any NA entries will be set to zero. Defaults to NULL, in which case starting values will be obtained from the partial proportional odds model as implemented in the package \pkg{ordinal} using \code{\link[ordinal]{clm}}.
##' @param hessianMeth Method to compute the Hessian of the optimisation. Passed to the argument finalHessian in \code{\link[maxLik]{maxLik}}.
##' @param returnMaxLik Should the function return the fit object from \code{\link[maxLik]{maxLik}} rather that the \code{\link{metaSDTreg}} object? Defaults to FALSE.
##' @param loglikUse For testing purposes. The C++ option is not presently included. Should the procedure use the log-likelihood implemented in R (default), or that in C++? Option will be dropped in future versions.
##' @param scoreFun For testing purposes. Should the optimisation procedure use the analytic score function? If FALSE (and 'meth' relies on gradients) numerical methods are used. Presently, only numerical methods are implemented. Defaults to FALSE. Option default may change in future versions.
##' @param c.maxLik A list that is passed to the 'MaxControl' object in \code{\link[maxLik]{maxLik}}. The list should contain control arguments to be passed to the maximisation method used by maxLik (the default being maxBFGS). In particular, 'iterlim' is the number of iterations (default 200) and a larger number of 'printlevel' (default 0) prints more information from the maximisation procedure.
##' @param c.clm A list that is passed to clm.control containing controls to be used by \code{\link[ordinal]{clm}}. The clm is used to determine starting values. By default, maxIter is set to 1000 and maxModIter is set to 20.
##' @param ... For future methods.
##'
##' @return A list of class 'metaSDTcontrol' containing the control arguments.
##'
##' @examples
##' names(metaSDTcontrol())
##'
##' @export
##' 
metaSDTcontrol <- function(meth = "BFGS", start.vals = NULL, hessianMeth = TRUE, 
    returnMaxLik = FALSE, loglikUse = c("R", "C++"), scoreFun = FALSE, c.maxLik = list(), 
    c.clm = list(maxIter = 1000, maxModIter = 20), ...) {
    
    out <- list(meth = meth, start.vals = start.vals, hessianMeth = hessianMeth, 
                returnMaxLik = returnMaxLik, loglikUse = loglikUse[1], scoreFun = scoreFun, 
                c.maxLik = do.call(maxLik::maxControl, args = c.maxLik),
                c.clm = do.call(ordinal::clm.control, 
                                args = c.clm), ...)
    
    class(out) <- "metaSDTcontrol"
    
    return(out)
}

##' Print method for metaSDTreg
##' 
##' @param x An object of class metaSDTreg.
##' @param ... For future methods.
##'
##' @return Invisible.
##'
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Fit model to subset of data
##' (fit <- metaSDTreg(A ~ signal,
##'             data=metadata,
##'             subset = m <= 20))
##'
##' @method print metaSDTreg
##'
##' @export
##' 
print.metaSDTreg <- function(x, ...) {
    
    print(summary(as.list(x)))
}


##' Coefficients method for metaSDTreg
##' 
##' @param object An object of class metaSDTreg.
##' @param ... For future methods.
##'
##' @return A named vector of parameter estimates.
##'
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Fit model to subset of data
##' fit <- metaSDTreg(A ~ signal,
##'             data=metadata,
##'             subset = m <= 20)
##' coef(fit)
##' 
##' @export
##' 
coef.metaSDTreg <- function(object, ...) {
    
    object$coefficients
}

##' Variance-covariance method for metaSDTreg
##' 
##' @param object An object of class metaSDTreg.
##' @param ... For future methods.
##'
##' @return A matrix of variances and covariances.
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
##' ## Standard errors
##' sqrt(diag(vcov(fit)))
##'
##' 
##' @export
##' 
vcov.metaSDTreg <- function(object, ...) {
    
    object$vcov
}

##' Summary method for metaSDTreg
##' 
##' @param object An object of class metaSDTreg.
##' @param ... For future methods.
##'
##' @return A list with class 'summary.metaSDTreg' containing summaries.
##'
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Fit model to subset of data
##' fit <- metaSDTreg(A ~ signal,
##'             data=metadata,
##'             subset = m <= 20)
##' summary(fit)
##'
##' @method summary metaSDTreg
##' 
##' @export
##' 
summary.metaSDTreg <- function(object, ...) {
    
    se <- sqrt(diag(object$vcov))
    ests <- object$coefficients
    z.vals <- ests/se
    
    TAB <- cbind(Estimate = ests, `Std. Error` = se, `z value` = z.vals, `P(>|z|)` = 2 * 
        (1 - stats::pnorm(abs(z.vals))))
    
    out <- list()
    out$call <- object$call
    out$logLik <- object$logLik
    out$coefficients <- TAB
    
    
    class(out) <- "summary.metaSDTreg"
    return(out)
}

##' Print summary method for metaSDTreg
##' 
##' @param x An object of class summary.metaSDTreg.
##' @param ... For future methods.
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
##' summary(fit)
##'
##' @method print summary.metaSDTreg
##'
##' @export
##' 
print.summary.metaSDTreg <- function(x, ...) {
    
    cat("Call: \n")
    print(x$call)
    cat("\n", paste("Log likelihood:", round(x$logLik, 2)), "\n \n")
    
    stats::printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
    
}

##' Log-likelihood of metaSDTreg
##'
##' Extract the log-likelihood from a metaSDTreg object.
##'
##' @param object Object of class 'metaSDTreg'.
##' @param ... For further methods.
##'
##' @return Numeric, the likelihood at the maximum found by the optimisation procedure.
##' 
##' @examples
##' ## Declare simulated data as metaSDTdata
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' ## Fit model to subset of data
##' fit <- metaSDTreg(A ~ signal,
##'             data=metadata,
##'             subset = m <= 20)
##' logLik(fit)
##' 
##' @export
##'
logLik.metaSDTreg <- function(object, ...) {
    return(object$logLik)
}

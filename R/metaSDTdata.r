##' Construct metaSDTdata
##' 
##' Constructor function for a metaSDTdata object.
##'
##' If type1 or type 2 is not an ordered factor, the function returns a warning. The function constructs a data frame containing variables named c('type1','type2','A','signal') along with any other variable in 'data' that is not given as an argument to the function. Because of this a warning is issued when the names c('type1','type2','A','signal') are present in 'data'.
##'
##' @param data Data frame to be converted. Data should be in long format with one row corresponding to a single response.
##' @param type1 A string naming the variable containing the type 1 response, which should be an ordered factor with two levels where the first level corresponds to 'noise',
##' @param type2 A string naming the variable containing the ordinal type 2 response, which should be an ordered factor.
##' @param signal A string naming the variable containing the signal.
##' 
##' @return A data object of class 'metaSDTdata'. This has attributes 'L', the number of levels in the ordinal type 2 rating, and 'K' which is two times L (the number of levels of the ordinal variable 'A').
##'
##' @examples
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##'
##' summary(metadata)
##' summary.data.frame(metadata)
##' 
##' @export
##' 
metaSDTdata <- function(data, type1, type2, signal) {
    ## Checks
    if (!inherits(data, "data.frame")) {
        stop("data must be a data frame.")
    }
    if (!is.factor(data[[type1]]) | !is.ordered(data[[type1]]) | nlevels(data[[type1]] != 
        2)) {
        warning(paste(type1, "is not an ordered factor with two levels."))
    }
    if (!is.factor(data[[type2]]) | !is.ordered(data[[type2]])) {
        warning(paste(type2, "is not an ordered factor."))
    }
    if (length(intersect(c("type1", "type2", "A", "signal"), names(data))) > 0) {
        warning("Variables 'type1', 'type2', 'signal' and 'A' are overwritten.")
    }
    
    ## Create output data frame containing the essential variables
    out <- subset(data, select = c(type1 = type1, type2 = type2, signal = signal))
    names(out) <- c("type1", "type2", "signal")
    
    ## Add any covariates
    out <- data.frame(out, subset(data, select = setdiff(names(data), c(type1, type2, 
        signal))))
    
    ## Create ordinal variable
    out$type1 <- factor(as.numeric(out$type1) - 1)  # response is 0,1
    out$type2 <- factor(as.numeric(out$type2), ordered = T)  # ordered 1,2,...,L
    out$A <- ifelse(out$type1 == 0, nlevels(out$type2) - as.numeric(out$type2) + 
        1, as.numeric(out$type2) + nlevels(out$type2))  # if resp 0/1
    out$A <- factor(out$A, ordered = T)
    
    ## Set attributes
    attr(out, "L") <- nlevels(out$type2)
    attr(out, "K") <- nlevels(out$type2) * 2
    
    ## Set class and output
    class(out) <- append("metaSDTdata", class(data))
    return(out)
}


##' Print method for metaSDT data
##'
##' @param x The metaSDTdata object.
##' @param ... Additional arguments passed to print.data.frame
##'
##' @return Invisible.
##'
##' @examples 
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##' print(metadata)
##'
##' @method print metaSDTdata
##'
##' @export
##' 
print.metaSDTdata <- function(x, ...) {
    print.data.frame(x)
}


##' Summarise a metaSDTdata object as a cognitive experiment.
##'
##' @param object The metaSDTdata object.
##' @param ... Additional arguments. Presently not used.
##'
##' @return A list with class 'summary.metaSDTdata' containing summaries.
##'
##' @examples 
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##' summary(metadata)
##'
##' @method summary metaSDTdata
##' 
##' @export
##' 
summary.metaSDTdata <- function(object, ...) {
    out <- list()
    out$tab <- table(A = object$A, signal = object$signal)
    out$L <- attr(object, "L")
    out$slevs <- unique(object$signal)
    
    class(out) <- "summary.metaSDTdata"
    return(out)
}


##' Print method for a summary.metaSDTdata object
##' 
##' @param x A metaSDTdata object.
##' @param ... For future methods.
##'
##' @return Invisible.
##'
##' @examples 
##' metadata <- metaSDTdata(simMetaData, type1='resp', type2='conf', signal='S')
##' summary(metadata)
##'
##' @method print summary.metaSDTdata
##'
##' @export
##' 
print.summary.metaSDTdata <- function(x, ...) {
    cat("Data from a metacognitive experiment. \n", paste("\n Type 2 task response has", 
        x$L, "levels"), paste("\n Signal variable has levels: \n"), unique(x$slevs), 
        paste("\n Responses observed: \n"))
    print.table(stats::addmargins(x$tab))
}

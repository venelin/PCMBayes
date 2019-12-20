#' Check if an object is a ParameterPrior object.
#' @param o an R object.
#' @return a logical indicating whether the S3 class of \code{o} inherits from
#' \code{'ParameterPrior'}.
#' @export
is.ParameterPrior <- function(o) { inherits(o, "ParameterPrior") }

#' Create a ParameterPrior object
#' @param d,r character strings denoting the names of a density function and of
#' a random generator function for some distribution, e.g. "dunif" and "runif".
#' Alternatively, these can be function objects, but this is not recommended,
#' because it can result in multiple duplication of these objects in memory.
#' Default values: \code{NULL}.
#' @param p a list with 'static' arguments common for \code{d} and \code{r}.
#' By 'static', it is meant arguments that do not change during a sampling
#' procedure, e.g. the rate of an exponentional prior or the boundaries of a
#' uniform prior. Default value: \code{NULL}
#' @param p.d a list of static arguments to be passed specifically to
#' \code{d} (added to \code{p}). Default value: \code{NULL}.
#' @param p.r a list of static arguments to be passed specifically to
#' \code{r} (added to \code{p}). Default value: \code{NULL}.
#' @param ... additional arguments to be stored as member elements in the
#' returned ParameterPrior object.
#'
#' @return an object of S3 class \code{c("ParameterPrior", "Prior")}.
#' @export
ParameterPrior <- function(
  d = NULL, r = NULL, p = NULL, p.d = NULL, p.r = NULL, ...) {
  structure(
    c(as.list(environment()), list(...)), class = c("ParameterPrior", "Prior"))
}

#' @export
PriorDensity.ParameterPrior <- function(prior, vecParams, log = TRUE, ...) {
  if(is.function(prior$d)) {
    do.call(prior$d, c(list(v = vecParams, log = log), prior$p, prior$p.d))
  } else if(is.character(prior$d)) {
    d <- get(prior$d)
    if(!is.function(d)) {
      stop(paste0("PriorDensity.ParameterPrior:: ", prior$d,
                  "could not be found."))
    } else {
      do.call(d, c(list(vecParams, log = log), prior$p, prior$p.d))
    }
  } else {
    stop(
      paste0("PriorDensity.ParameterPrior:: prior$d should be a function ",
             " or a character string denoting a function."))
  }
}

#' Check if an object is a Prior object.
#' @param o an R object.
#' @return a logical indicating whether the S3 class of \code{o} inherits from
#' \code{'Prior'}.
#' @export
is.Prior <- function(o) { inherits(o, "Prior") }


#' Evaluate a prior distribution density over parameter values.
#'
#' @details This is an S3 generic function. Call
#' \code{methods("PCMPriorDensity")} to see implementing methods.
#'
#' @param prior a Prior object.
#' @param vecParams a numerical vector.
#' @param log logical indicating if a log-prior should be returned
#' (default: TRUE).
#' @param returnVector logical indicating if a vector of length
#' \code{length(vecParams)} should be returned. Default: \code{TRUE}.
#' @param ... additional arguments passed to methods.
#' @return The returned value depends on \code{returnVector}:
#' \itemize{
#' \item If \code{returnVector} is \code{TRUE} (default) it should be a numeric
#' vector of length equal to \code{length(vecParams)}.
#' If the prior distributions for the elements in \code{vecParams},
#' are independent, then each element in the returned vector corresponds to the
#' (log-)prior density of the corresponding element in \code{vecParams}.
#' Otherwise, the first element in the returned vector corresponds to the
#' joint (log-)prior density and the other elements are set to 1 (if
#' \code{log=FALSE}) or to 0 (if \code{log=TRUE}).
#' \item If \code{returnVector} is \code{FALSE}, this should be a single real
#' number corresponding to the joint prior (log-)density of the elements in
#' \code{vecParams}.
#' }
#' @export
PriorDensity <- function(prior, vecParams, log = TRUE, returnVector = TRUE...) {
  if(!is.Prior(prior)) {
    stop("PriorDensity:: prior argument should be a Prior object.")
  }
  UseMethod("PriorDensity", prior)
}

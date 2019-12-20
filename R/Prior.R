#' Check if an object is a Prior object.
#' @param o an R object.
#' @return a logical indicating whether the S3 class of \code{o} inherits from
#' \code{'Prior'}.
#' @export
is.Prior <- function(o) { inherits(o, "Prior") }

#' Add a prior to an object
#'
#' This is an S3 generic with methods depending on the S3 class of the \code{o},
#' argument.
#'
#' @inheritParams PCMBase::PCMAddToListAttribute
#' @param o an object.
#' @param prior a prior object.
#' @return if inplace is TRUE (default) nothing is returned. Otherwise, a
#' modified version of \code{o} is returned.
#'
#' @examples
#' library(PCMBase)
#' model <- PCMBaseTestObjects$model_MixedGaussian_ab
#'
#' AddPrior(
#'   model,
#'   prior = ParameterPrior(
#'     d = "dunif", r = "runif",
#'     p = list(min = PCMParamGetShortVector(PCMParamLowerLimit(model)),
#'              max = PCMParamGetShortVector(PCMParamUpperLimit(model)))))
#'
#' AddPrior(
#'   model, "a$Sigma_x", enclos = "diag(?[,,1])",
#'   prior = ParameterPrior(
#'     d = "dunif", r = "runif",
#'     p = list(min = c(0.2, 0.1, 0.1), max = c(5, 5, 5))))
#'
#' modelPrior <- PCMPrior(model)
#' vec <- PCMParamGetShortVector(model)
#' PriorDensity(modelPrior, vec)
#' @export
AddPrior <- function(
  o, member = "", enclos = "?", prior, inplace = TRUE) {
  UseMethod("AddPrior", o)
}

#' Evaluate a prior distribution density over parameter values.
#'
#' @details This is an S3 generic function. Call
#' \code{methods("PCMPriorDensity")} to see implementing methods.
#'
#' @param prior a PCMPrior object.
#' @param vecParams a numerical vector.
#' @param log logical indicating if a log-prior should be returned
#' (default: TRUE).
#' @param ... additional arguments passed to methods.
#' @return The returned value depends on the S3 class of \code{prior}.
#' \itemize{
#' \item For any class except \code{'PCMPrior'}, it should be a numeric
#' vector of length equal to \code{length(vecParams)}.
#' If the the prior distributions for the elements in \code{vecParams},
#' are independent, then each element in the returned vector corresponds to the
#' (log-)prior density of the corresponding element in \code{vecParams}.
#' Otherwise, the first element in the returned vector corresponds to the
#' joint (log-)prior density and the other elements are set to 1 (if
#' \code{log=FALSE}) or to 0 (if \code{log=TRUE}).
#' \item For class \code{'PCMPrior'}, this should be a single real
#' number corresponding to the joint prior (log-)density of the elements in
#' \code{vecParams}.
#' }
#' @export
PriorDensity <- function(prior, vecParams, log = TRUE, ...) {
  if(!is.Prior(prior)) {
    stop("PriorDensity:: prior argument should be a Prior object.")
  }
  UseMethod("PriorDensity", prior)
}

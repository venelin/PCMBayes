#' Check if an object is a PCMPrior object.
#' @param o an R object.
#' @return a logical indicating whether the S3 class of \code{o} inherits from
#' \code{'PCMPrior'}.
#' @export
is.PCMPrior <- function(o) { inherits(o, "PCMPrior") }

#' Create a PCMPrior object for a PCM model object
#' @param model a PCM object.
#' @param ...  additional arguments to be stored as member elements in the
#' returned PCMPrior object.
#' @return an object of S3 class \code{c("PCMPrior", "CompositePrior", "Prior")}.
#' @export
PCMPrior <- function(model, ...) {
  structure(
    list(priors = PCMCombineListAttribute(model, "prior"), ...),
    class = c("PCMPrior", "CompositePrior", "Prior"))
}

#' Add a prior to a PCM object
#'
#' @inheritParams PCMBase::PCMAddToListAttribute
#' @param o a PCM object.
#' @param prior a prior object.
#' @return if inplace is TRUE (default) nothing is returned. Otherwise, a
#' modified version of \code{o} is returned.
#'
#' @examples
#' library(PCMBase)
#' model <- PCMBaseTestObjects$model_MixedGaussian_ab
#'
#' PCMAddPrior(
#'   model,
#'   prior = ParameterPrior(
#'     d = "dunif", r = "runif",
#'     p = list(min = PCMParamGetShortVector(PCMParamLowerLimit(model)),
#'              max = PCMParamGetShortVector(PCMParamUpperLimit(model)))))
#'
#' PCMAddPrior(
#'   model, "a$Sigma_x", enclos = "diag(?[,,1])",
#'   prior = ParameterPrior(
#'     d = "dunif", r = "runif",
#'     p = list(min = c(0.2, 0.1, 0.1), max = c(5, 5, 5))))
#'
#' modelPrior <- PCMPrior(model)
#' vec <- PCMParamGetShortVector(model)
#' PriorDensity(modelPrior, vec)
#' @export
PCMAddPrior <- function(
  o, member = "", enclos = "?", prior, inplace = TRUE) {

  if(!is.PCM(o)) {
    stop(
      "PCMAddPrior:: The argument o should inherit from S3 class PCM.")
  }
  if(!is.Prior(prior)) {
    stop(
      "PCMAddPrior:: The argument prior should inherit from S3 class Prior.")
  }

  if(inplace) {
    eval(substitute(PCMAddToListAttribute(
      name = "prior", value = prior, object = o, member = member,
      enclos = enclos, inplace = TRUE)),
      parent.frame())
  } else {
    o <- PCMAddToListAttribute(
      name = "prior", value = prior, object = o, member = member,
      enclos = enclos, inplace = FALSE)
    o
  }
}



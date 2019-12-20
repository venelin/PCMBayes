#' Check if an object is a MGPMStatePrior object.
#' @param o an R object.
#' @return a logical indicating whether the S3 class of \code{o} inherits from
#' \code{'MGPMStatePrior'}.
#' @export
is.MGPMStatePrior <- function(o) { inherits(o, "MGPMStatePrior") }

#' Create a MGPMStatePrior object for a MGPMState object
#' @param state a MGPMState object.
#' @param ...  additional arguments to be stored as member elements in the
#' returned MGPMStatePrior object.
#' @return an object of S3 class \code{c("MGPMStatePrior", "Prior")}.
#' @export
MGPMStatePrior <- function(state, ...) {
  structure(
    list(priors = PCMCombineListAttribute(state, "prior"), ...),
    class = c("MGPMStatePrior", "Prior"))
}


#' @export
AddPrior.MGPMState <- function(
  o, member = "", enclos = "?", prior, inplace = TRUE) {

  if(!is.Prior(prior)) {
    stop(
      "AddPrior.MGPMState:: The argument prior should inherit from S3 class Prior.")
  }

  if(!member %in% c())




  if(inplace) {
    eval(substitute(PCMAddToListAttribute(
      name = "prior", value = prior, object = o, member = member,
      enclos = enclos, spec = FALSE, fixed = TRUE, inplace = TRUE)),
      parent.frame())
  } else {
    o <- PCMAddToListAttribute(
      name = "prior", value = prior, object = o, member = member,
      enclos = enclos, spec = FALSE, fixed = TRUE, inplace = FALSE)
    o
  }
}

#' @export
PriorDensity.MGPMStatePrior <- function(prior, vecParams, log = TRUE, ...) {
  vecDens <- rep(NA_real_, length(vecParams))
  if(!is.list(prior$priors)) {
    stop(paste0("PriorDensity.MGPMStatePrior:: the member prior$priors ",
                "should be a list of Prior objects."))
  } else {
    for(pr in prior$priors) {
      vecDens[pr$pos] <- PriorDensity(pr, vecParams[pr$pos], log = log, ...)
    }
    if(log) {
      sum(vecDens)
    } else {
      prod(vecDens)
    }
  }
}

#' Check if an object is a Prior object.
#' @param o an R object.
#' @return a logical indicating whether the S3 class of \code{o} inherits from
#' \code{'Prior'}.
#' @export
is.Prior <- function(o) { inherits(o, "Prior") }

#' Create a ParameterPrior object
#' @param enclos a character string containing the special symbol '?'. This
#' symbol is to be replaced by expressions denoting objects to which this
#' ParameterPrior object will be added. The result of this
#' substitution can be anything but, usually would be a valid R expression.
#' Default: "?". See also \code{\link{AddPrior}}.
#' @param pos an integer denoting an offset position in a parameter vector.
#' Default: \code{NULL}.
#' @param d,r character strings denoting the names of a density function and of
#' a random generator function for some distribution, e.g. "dunif" and "runif".
#' @param p a list with 'static' arguments common for \code{d} and \code{r}.
#' By 'static', it is meant arguments that do not change during a sampling
#' procedure, e.g. the rate of an exponentional prior or the boundaries of a
#' uniform prior.
#' @param p.d a list of static arguments to be passed specifically to
#' \code{d} (added to \code{p}).
#' @param p.r a list of static arguments to be passed specifically to
#' \code{r} (added to \code{p}).
#' @param ... additional arguments to be stored as member elements in the
#' returned ParameterPrior object.
#'
#' @return an object of S3 class \code{c("ParameterPrior", "Prior")}.
#' @export
ParameterPrior <- function(
  enclos = "?", pos = NULL,
  d = NULL, r = NULL, p = NULL, p.d = NULL, p.r = NULL, ...) {
  structure(
    c(as.list(environment()), list(...)), class = c("ParameterPrior", "Prior"))
}

#' Check if an object is a ParameterPrior object.
#' @param o an R object.
#' @return a logical indicating whether the S3 class of \code{o} inherits from
#' \code{'ParameterPrior'}.
#' @export
is.ParameterPrior <- function(o) { inherits(o, "ParameterPrior") }

#' Create a PCMPrior object for a PCM model object
#' @param model a PCM object.
#' @param ...  additional arguments to be stored as member elements in the
#' returned PCMPrior object.
#' @return an object of S3 class \code{c("PCMPrior", "Prior")}.
#' @export
PCMPrior <- function(model, ...) {
  structure(
    list(priors = PCMCombineListAttribute(model, "prior"), ...),
    class = c("PCMPrior", "Prior"))
}

#' Check if an object is a PCMPrior object.
#' @param o an R object.
#' @return a logical indicating whether the S3 class of \code{o} inherits from
#' \code{'PCMPrior'}.
#' @export
is.PCMPrior <- function(o) { inherits(o, "PCMPrior") }

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

#' @export
AddPrior.PCM <- function(
  o, member = "", enclos = "?", prior, inplace = TRUE) {

  if(!is.Prior(prior)) {
    stop(
      "AddPrior:: The argument prior should inherit from S3 class Prior.")
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

#' @export
PriorDensity.PCMPrior <- function(prior, vecParams, log = TRUE, ...) {
  vecDens <- rep(NA_real_, length(vecParams))
  if(!is.list(prior$priors)) {
    stop(paste0("PriorDensity.PCMPrior:: the member prior$priors ",
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



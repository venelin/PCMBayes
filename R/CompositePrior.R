#' @export
PriorDensity.CompositePrior <- function(
  prior, vecParams, log = TRUE, returnVector = TRUE, ...) {
  vecDens <- rep(NA_real_, length(vecParams))
  if(!is.list(prior$priors)) {
    stop(paste0("PriorDensity.CompositePrior:: the member prior$priors ",
                "should be a list of Prior objects."))
  } else {
    for(pr in prior$priors) {
      vecDens[pr$pos] <-
        PriorDensity(pr, vecParams[pr$pos], log = log, returnVector = TRUE, ...)
    }
    dens <- if(log) {
      sum(vecDens)
    } else {
      prod(vecDens)
    }

    if(returnVector) {
      vecDens[1L] <- dens
      if(log) {
        vecDens[-1L] <- 0.0
      } else {
        vecDens[-1L] <- 1.0
      }
      vecDens
    } else {
      dens
    }
  }
}

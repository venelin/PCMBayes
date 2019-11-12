# PCMMCMC <- function(X, tree, prior, propose) {
#
#   MetropolisHastings <- function(numIterations) {
#     for(i in i + seq_len(numIterations)) {
#       proposal <- PCMMove()
#
#       new.pars <- proposal$pars
#       #new.cache <- prop$prop$cache
#       #accept.type <- c(accept.type,paste(prop$decision,prop$move,sep="."))
#       propName <- paste(prop$decision,prop$move,sep=".")
#       acceptN[propName] <- acceptN[propName]+1
#       pr2 <- prior(new.pars,cache)
#       hr <- prop$hr
#       new <- lik.fn(new.pars, cache, dat, model=model)
#       nll <- new$loglik
#       nll <- ifelse(is.na(nll), -Inf, nll)
#       if (runif(1) < exp(nll-oll+pr2-pr1+hr)){
#         oldpar <- new.pars
#         pr1 <- pr2
#         oll <- nll
#         accept[propName] <- accept[propName]+1
#       } else {
#         #accept <- c(accept,0)
#       }
#     }
#   }
#
#   mcmc <- list(
#     X = X,
#     tree = tree,
#     prior = prior,
#     propose = propose,
#     states = list())
#   class(mcmc) <- c("PCMMCMC", class(mcmc))
#   mcmc
# }
#
#


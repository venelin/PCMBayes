#' Context of a MCMC
#' @import PCMBase
#' @inheritParams PCMBase::PCMLik
#' @param model a template used to build MixedGaussian model objects.
#' @return an object of S3 class 'MCMCContext'.
#' @export
MCMCContext <- function(
  X, tree, model, SE = matrix(0, PCMNumTraits(model), PCMTreeNumTips(tree))) {

  ctx <- list(
    k = nrow(X),
    X = X,
    SE = SE,
    model = model,
    treeOriginal = tree,
    tree = PCMTree(tree))
  PCMTreeSetLabels(ctx$tree)
  class(ctx) <- "MCMCContext"
  ctx
}

NamesModelTypes <- function(modelTemplate) {
  names(modelTemplate)[sapply(modelTemplate, is.PCM)]
}

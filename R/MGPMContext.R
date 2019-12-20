#' Context of a MGPM used for MCMC sampling
#' @import PCMBase
#' @inheritParams PCMBase::PCMLik
#' @param model a MGPM template used to build MGPM model objects.
#' @return an object of S3 class 'MGPMContext'.
#' @export
MGPMContext <- function(
  X, tree, model, SE = matrix(0, PCMNumTraits(model), PCMTreeNumTips(tree))) {

  ctx <- list(
    k = nrow(X),
    X = X,
    SE = SE,
    model = model,
    treeOriginal = tree,
    tree = PCMTree(tree))
  PCMTreeSetLabels(ctx$tree)
  class(ctx) <- "MGPMContext"
  ctx
}

NamesModelTypes <- function(modelTemplate) {
  names(modelTemplate)[sapply(modelTemplate, is.PCM)]
}

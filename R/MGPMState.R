#' Context of an MGPM MCMC
#' @import PCMBase
#' @inheritParams PCMBase::PCMLik
#' @param model a MGPM template used to build MGPM model objects.
#' @return a MGPM context object.
#' @export
MGPMContext <- function(
  X, tree, model, SE = matrix(0, PCMNumTraits(model), PCMTreeNumTips(tree))) {

  ctx <- list(
    k = nrow(X),
    X = X,
    SE = SE,
    modelTypes = model[sapply(model, is.PCM)],
    argsMixedGaussian = model[!sapply(model, is.PCM)],
    treeOriginal = tree,
    tree = PCMTree(tree))
  PCMTreeSetLabels(ctx$tree)
  class(ctx) <- "MGPMContext"
  ctx
}

#' (Re-)construct an MGPMState object from a state vector
#'
#' @param s a numerical vector representing an MGPM state. Some of the entries
#' in this vector are interpreted as integers denoting nodes in a tree,
#' regimes or model types (see Details).
#' @param ctx a MGPMContext object.
#' @return a MGPMState object.
#' @seealso \code{\link{MGPMStateVector}}
#' @export
MGPMState <- function(s, ctx) {
  # number of tips
  N <- PCMTreeNumTips(ctx$tree)

  # number of shifts
  K <- as.integer(s[1L])

  # shift nodes
  n <- s[1L + seq_len(K)]

  # shift locations on the branches leading to the shift nodes
  l <- s[1L + K + seq_len(K)]

  # regimes: these are integers among n
  r <- as.integer(s[1L + K + K + seq_len(K)])

  # order n, l and r according to PCMTreePreorder
  ordern <- order(match(n, ctx$tree$edge[ctx$preorderNodes, 2]))

  n <- n[ordern]
  l <- l[ordern]
  r <- r[ordern]

  regimeNames <- unique(c(N+1L, r))

  # Number of different regimes
  R <- length(regimeNames)

  # The node labels in ctx$tree must be 1,2,...,M, with N+1 being the root and
  # 1,...,N being the tips.
  tree <- ctx$tree

  minSingletonBranchLength <- getOption("PCMBayes.min.singleton.length", 0.1)
  # Insert singleton nodes wherever needed
  needASingleton <- (l >= minSingletonBranchLength)

  if(sum(needASingleton) > 0L) {
    # ids of rows in tree$edge and tree$edge.length
    branchesNeedingSingletons <- match(n[needASingleton], tree$edge[, 2L])
    # Offsets of the singletons in root-ward direction
    posSingletons <-
      tree$edge.length[branchesNeedingSingletons] - l[needASingleton]

    if(sum(posSingletons > minSingletonBranchLength) > 0L)
      tree <- PCMTreeInsertSingletons(
        tree,
        n[needASingleton][posSingletons >= minSingletonBranchLength],
        posSingletons[posSingletons >= minSingletonBranchLength])
  }

  # After possible insertion of singleton nodes, the node indices are not
  # anymore the same, so we use the node labels to set the partition
  part.names <- as.character(c(N+1L, n))
  part.regime <- structure(character(length(part.names)), names = part.names)
  part.regime[] <- as.character(c(N+1L, r))
  PCMTreeSetPartRegimes(tree, part.regime, setPartition = TRUE)

  # model types associated with the regimes
  m <- s[1L + K + K + K + seq_len(R)]
  names(m) <- regimeNames

  # Create the MGPM object
  model <- do.call(
    MixedGaussian,
    c(list(k = ctx$k,
           modelTypes = ctx$modelTypes,
           mapping = m),
      ctx$argsMixedGaussian))

  # Load the parameter values into the model (the number of parameters is stored
  # in P)
  P <- PCMParamLoadOrStore(
    model, s, offset = 1L + K + K + K + R, k = ctx$k, R = R,
    load = TRUE)

  structure(list(
    K = K, R = R, P = P, n = n, l = l, r = r, m = m, model = model, tree = tree),
    class = "MGPMState")
}

#' Construct an MGPM state vector.
#'
#' @param K an integer (see Details).
#' @param n an integer vector of length \code{K} denoting shift nodes (see Details).
#' @param l a real vector of length \code{K} (see Details).
#' @param r an integer vector of length \code{K} (see Details).
#' @param m an integer vector of length \code{R}, where \code{R} denotes the
#' total number of regimes (see Details).
#' @param v a real vector of length \code{P}, where \code{P} denotes the total
#' number of numeric parameters of the MGPM model represented by the state
#' (see Details).
#' @return a numerical vector representing the concatenation of
#' \code{K, n, l, r, m, v}, with S3 class set to 'MGPMStateVector'.
#'
#' @details
#' The MGPMState of an MCMC is represented by a numerical vector of the form:
#' \deqn{\vec{s}=(K, n_2,...,n_K, l_2,...,l_{K+1}, r_2,...,r_{K+1}, m_1,...,m_R, v_1,...,v_P)^T}
#' The element of this vector are described as follows:
#' \describe{
#' \item{K: }{number of shifts;}
#' \item{\eqn{(n_2,...,n_{K+1}): }}{shift nodes - the shifts in the model occur
#' at points within the branches leading to the shift nodes in tip-ward
#' direction. The corresponding locations of these points are specified by
#' \eqn{(l_2,...,l_{K+1})}. The nodes \eqn{(n_2,...,n_{K+1})} should be ordered
#' according to \code{\link{PCMTreePreorder}(ctx$tree)}.}
#' \item{\eqn{(l_2,...,l_{K+1}): }}{offset of the shift points measured as
#' distances from the beginnings of the branches leading to shift nodes in
#' tip-ward direction.}
#' \item{\eqn{(r_2,...,r_{K+1}): }}{regime index vector. This is an integer
#' vector with elements among \eqn{(N+1,n_2,...,n_{K+1})}, indicating the regime
#' associated with each part in the tree. The regimes are named as the shift
#' nodes, with N+1 corresponding to the part (and regime) originating at the
#' root. It is possible to have lumped regimes, that is, different parts having
#' the same regime. This regime-lumping must obey the following rules:
#' \enumerate{
#' \item neighbor parts cannot have a lumped regime. Two parts originating at
#' nodes \eqn{n_i} and \eqn{n_j} in the tree are called neighbor parts if they
#' are separated solely by \eqn{n_i} or by \eqn{n_j};
#' \item to resolve the conflict between the shift nodes of the different parts
#' covered by a lumped regime, it is established that the name of a lumped
#' regime must equal the shift-node that appears first, according to
#' \code{\link{PCMTreePreorder}(ctx$tree)}.
#' }}
#' \item{\eqn{(m_1,...,m_R): }}{model type assignment to the unique regimes.
#' This is an integer vector with elements between 1 and M, M denoting the
#' total number of model types possible. Each element corresponds to an
#' element in \code{unique(c(N+1,r_2,...,r_{K+1}))}}
#' \item{\eqn{(v_1,...,v_P):} }{real numbers passed to
#' \code{\link{PCMParamLoadOrStore}}. This is a vectorized form of all model
#' parameters.}
#' }
#' @seealso \code{\link{MGPMState}}
#' @export
MGPMStateVector <- function(K, n, l, r, m, v) {
  structure(
    c(K, n, l, r, m, v), class = c("MGPMStateVector"))
}



# set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")
# tree <- PCMTree(ape::rtree(20))
# PCMTreeSetPartition(tree, c(26, 28))
# PCMTreeGetPartition(tree)
# order <- PCMTreePreorder(tree)
# order(match(c(28, 26), tree$edge[order, 2]))

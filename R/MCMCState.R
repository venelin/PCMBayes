#' Check if an object is of S3 class 'MCMCState'
#' @param o an object.
#' @return a logical.
#' @export
is.MCMCState <- function(o) {
  inherits(o, "MCMCState")
}

#' If necessary, convert an object to an MCMCState object.
#' @param s a \code{\link{MCMCState}} or a \code{\link{MCMCStateVector}}
#' object.
#' @param ctx a \code{\link{MCMCContext}} object. This is needed only if
#'  \code{s} is \code{\link{MCMCStateVector}} object.
#'
#' @return an MCMCState object corresponding to \code{s}. If \code{s} is already
#' a \code{\link{MCMCState}} object, it is returned as is. Otherwise, a
#' \code{\link{MCMCState}} object is constructed from \code{s} and \code{ctx}.
#' @export
as.MCMCState <- function(s, ctx = NULL) {
  if(is.MCMCState(s)) {
    s
  } else if(is.MCMCStateVector(s)) {
    if(!is.MCMCContext(ctx)) {
      stop(paste0(
        "as.MCMCState: ctx must be a MCMCContext object if s is a ",
        "MCMCStateVector."))
    }
    s <- MCMCState(s, ctx)
  } else {
    stop("as.MCMCState: s be a MCMCState or a MCMCStateVector object.")
  }
}

#' Construct a MCMCState object from a state vector
#'
#' @param s a numerical vector representing a MCMC state. Some of the
#'  entries in this vector are interpreted as integers denoting nodes in a
#'  tree, regimes or model types (see \code{\link{MCMCStateVector}}).
#' @param ctx a MCMCContext object.
#' @return a MCMCState object.
#'
#' @seealso \code{\link{MCMCStateVector}}
#' @export
MCMCState <- function(s, ctx) {
  # number of tips
  N <- PCMTreeNumTips(ctx$tree)

  # Initialize default values:
  K <- 0L
  R <- 1L
  P <- NA_integer_
  n <- integer(0L)
  l <- double(0L)
  r <- N+1L
  m <- 1L
  v <- NULL
  model <- NULL
  # The node labels in ctx$tree must be 1,2,...,M, with N+1 being the root and
  # 1,...,N being the tips.
  tree <- ctx$tree


  if(length(s) >= 1L) {
    # number of shifts
    K <- as.integer(s[1L])
    if(K > 0L && !(length(s) >= 1L + K + K + K)) {
      stop(paste0("MCMCState:: The vector s is wrong length, given that K=", K,
                  ". It should have at least ", 1L + K + K + K, " elements."))
    } else {
      # shift nodes
      n <- as.integer(s[1L + seq_len(K)])

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

      m <- rep(1L, R)
      # model types associated with the regimes
      if(length(s) >= 1L + K + K + K + R) {
        m <- as.integer(s[1L + K + K + K + seq_len(R)])
      }
      names(m) <- regimeNames

      # Create the MCMC object
      namesModelTypes <- NamesModelTypes(ctx$model)
      model <- do.call(
        MixedGaussian,
        c(list(k = ctx$k,
               modelTypes = ctx$model[namesModelTypes],
               mapping = m),
          attr(ctx$model, "spec")[
            setdiff(names(attr(ctx$model, "spec")), namesModelTypes)]))

      # Load the parameter values into the model (the number of parameters is
      # stored in P)
      P <- as.integer(attr(model, "p"))
      v <- double(P)
      if(length(s) >= 1L + K + K + K + R + P) {
        v <- s[-seq_len(1L + K + K + K + R)]
      }
      P1 <- PCMParamLoadOrStore(model, v, offset = 0, k = ctx$k, R = R,load = TRUE)
      if(P != P1) {
        stop(paste0(
          "MCMCState:: something is wrong with the attribute p of the model",
          " and the number of parameters returned by PCMParamLoadOrStore."))
      }
    }
  }

  structure(list(
    K = K, R = R, P = P, n = n, l = l, r = r, m = m,
    v = v, model = model, tree = tree),
    class = "MCMCState")
}



#' Check if an object is of S3 class 'MCMCStateVector'
#' @param o an object.
#' @return a logical.
#' @export
is.MCMCStateVector <- function(o) {
  inherits(o, "MCMCStateVector")
}

#' Construct a MCMC state vector.
#'
#' @param K an integer (see Details).
#' @param n an integer vector of length \code{K} denoting shift nodes (see Details).
#' @param l a real vector of length \code{K} (see Details).
#' @param r an integer vector of length \code{K} (see Details).
#' @param m an integer vector of length \code{R}, where \code{R} denotes the
#' total number of regimes (see Details).
#' @param v a real vector of length \code{P}, where \code{P} denotes the total
#' number of numeric parameters of the MCMC model represented by the state
#' (see Details).
#' @return a numerical vector representing the concatenation of
#' \code{K, n, l, r, m, v}, with S3 class set to 'MCMCStateVector'.
#'
#' @details
#' The MCMCState of an MCMC is represented by a numerical vector:
#' \deqn{\vec{s}=(K, n_2,...,n_{K+1}, l_2,...,l_{K+1}, r_2,...,r_{K+1}, m_1,...,m_R, v_1,...,v_P)^T}
#' The element of this vector are described as follows:
#' \describe{
#' \item{\eqn{K}: }{number of shifts;}
#' \item{\eqn{(n_2,...,n_{K+1})^T}: }{shift nodes - the shifts in the model occur
#' at points within the branches leading to the shift nodes in tip-ward
#' direction. The corresponding locations of these points are specified by
#' \eqn{(l_2,...,l_{K+1})^T}. The nodes \eqn{(n_2,...,n_{K+1})^T} should be ordered
#' according to \code{\link{PCMTreePreorder}(ctx$tree)}.}
#' \item{\eqn{(l_2,...,l_{K+1})^T}: }{offset of the shift points measured as
#' distances from the beginnings of the branches leading to shift nodes in
#' tip-ward direction.}
#' \item{\eqn{(r_2,...,r_{K+1})^T}: }{regime index vector. This is an integer
#' vector with elements among \eqn{(N+1,n_2,...,n_{K+1})^T}, indicating the regime
#' associated with each part in the tree. The regimes are named as the shift
#' nodes, with N+1 corresponding to the part (and regime) originating at the
#' root. It is possible to have lumped regimes, that is, different parts
#' of the tree having the same regime. This regime-lumping must obey the
#' following rules:
#' \enumerate{
#' \item neighbor parts cannot have a lumped regime. Two parts originating at
#' nodes \eqn{n_i} and \eqn{n_j} in the tree are called neighbor parts if they
#' are separated solely by \eqn{n_i} or by \eqn{n_j};
#' \item to resolve the conflict between the shift nodes of the different parts
#' covered by a lumped regime, it is established that the name of a lumped
#' regime must equal the shift-node that appears first, according to
#' \code{\link{PCMTreePreorder}(ctx$tree)}.
#' }}
#' \item{\eqn{(m_1,...,m_R)^T}: }{model type assignment to the unique regimes.
#' This is an integer vector with elements between 1 and M, M denoting the
#' total number of model types possible. Each element corresponds to an
#' element in \code{unique(c(N+1,r_2,...,r_{K+1}))}}
#' \item{\eqn{(v_1,...,v_P)^T}: }{real numbers passed to
#' \code{\link{PCMParamLoadOrStore}}. This is a vectorized form of all model
#' parameters.}
#' }
#' @seealso \code{\link{MCMCState}}
#' @export
MCMCStateVector <- function(K, n, l, r, m, v) {
  structure(
    c(as.double(K), as.double(n), as.double(l), as.double(r), as.double(m),
      as.double(v)),
    class = c("MCMCStateVector"))
}


#' @title Indices of different parts of an MCMCState in an MCMCStateVector
#'
#' @description The different parts of an MCMCState are described in
#' \code{\link{MCMCStateVector}}. \code{posK} returns the osition of the
#'  number of shifts (this is always equal to 1). See the section functions
#'  for the others.
#'
#' @param s a \code{\link{MCMCState}} or a \code{\link{MCMCStateVector}}
#' object.
#' @param ctx a \code{\link{MCMCContext}} object. This is needed only if
#'  \code{s} is \code{\link{MCMCStateVector}} object.
#'
#' @return an integer vector.
#'
#' @seealso \code{\link{MCMCState}}, \code{\link{MCMCStateVector}}
#' @export
posK <- function(s, ctx = NULL) {
  1L
}

#' @describeIn posK
#'
#' Positions of the shift nodes;
#'
#' @export
posn <- function(s, ctx = NULL) {
  s <- as.MCMCState(s, ctx)
  1L + seq_len(s$K)
}

#' @describeIn posK
#'
#' Positions of the offsets of the shift-points in tip-ward direction
#' relative to the beginnings of the branches leading to shift nodes;
#'
#' @export
posl <- function(s, ctx = NULL) {
  s <- as.MCMCState(s, ctx)
  1L + s$K + seq_len(s$K)
}

#' @describeIn posK
#'
#' Position of the regime id's for each shift node;
#'
#' @export
posr <- function(s, ctx = NULL) {
  s <- as.MCMCState(s, ctx)
  1L + s$K + s$K + seq_len(s$K)
}

#' @describeIn posK
#'
#' Positions of the model type id's for each regime;
#'
#' @export
posm <- function(s, ctx = NULL) {
  s <- as.MCMCState(s, ctx)
  1L + s$K + s$K + s$K + seq_len(s$R)
}

#' @describeIn posK
#'
#' Positions of the model parameters.
#'
#' @export
posv <- function(s, ctx = NULL) {
  s <- as.MCMCState(s, ctx)
  1L + s$K + s$K + s$K + s$R + seq_len(s$P)
}

#' @export
AddPrior.MCMCState <- function(
  o, member = "", enclos = "?", prior, inplace = TRUE) {

  if(!is.Prior(prior)) {
    stop(
      "AddPrior.MCMCState:: The argument prior should inherit from S3 class Prior.")
  }

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


# set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")
# tree <- PCMTree(ape::rtree(20))
# PCMTreeSetPartition(tree, c(26, 28))
# PCMTreeGetPartition(tree)
# order <- PCMTreePreorder(tree)
# order(match(c(28, 26), tree$edge[order, 2]))


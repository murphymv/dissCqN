

#' @title Networks of Shared Species
#' @description Generate interaction networks comprising only the shared species
#'   across two or more networks.
#' @param net A list of two or more networks to compare, supplied as matrices.
#' @param pairwise Logical, whether to compare networks pairwise (default),
#'   rather than considering species across multiple networks.
#' @param compare.sub Subsets of networks to compare pairwise. These should be
#'   supplied as a list of two sets of network names or indices. If only one set
#'   is supplied, this set is compared to all other networks in \code{net}. If
#'   more than two sets are supplied, only the first two are used.
#' @note If comparing networks pairwise, and subsets are not specified, the
#'   output will contain some network self-comparisons (i.e. redundancy).
#' @return A list of networks of shared species. If comparing pairwise, this
#'   will be of length \emph{a} * \emph{b} * 2 (with \emph{a} and \emph{b} being
#'   the sets of networks to compare), or if considering multiple networks, the
#'   length of the original list.
#' @examples
netShared <- function(net, pairwise = TRUE, compare.sub = NULL) {

  n <- net; cs <- compare.sub

  ## List of networks
  if (!is.list(n)) n <- list(n)
  if (is.null(names(n))) names(n) <- paste0("net", 1:length(n))

  if (length(n) > 1) {

    if (pairwise) {

      ## Networks to compare pairwise
      n2 <- n1 <- n

      ## Compare subsets?
      if (!is.null(cs)) {
        if (!is.list(cs)) cs <- list(cs)
        n1 <- n[cs[[1]]]
        n2 <- if (length(cs) > 1) n[cs[[2]]] else n[!names(n) %in% names(n1)]
      }

      ## Networks of shared species
      ns <- lapply(n1, function(i) {
        lapply(n2, function(j) {
          ss <- Map(function(k, l) k[k %in% l], dimnames(i), dimnames(j))
          lapply(list(i, j), function(k) {
            k[ss[[1]], ss[[2]], drop = FALSE]
          })
        })
      })
      ns <- lapply(rapply(ns, enquote, how = "unlist"), eval)
      setNames(ns, gsub("(.{1}$)", "_\\1", names(ns)))

    } else {

      ## Shared species across all networks
      spp <- unique(unlist(lapply(n, dimnames)))
      ss <- spp[sapply(spp, function(i) {
        all(sapply(n, function(j) i %in% unlist(dimnames(j))))
      })]

      ## Networks of shared species
      lapply(n, function(i) {
        s1 <- rownames(i) %in% ss
        s2 <- colnames(i) %in% ss
        i[s1, s2, drop = FALSE]
      })

    }

  } else n

}


#' @title Assemblage x Species Interaction Matrix
#' @description Generate matrix of assemblages x species interactions from a
#'   set of networks.
#' @param net An interaction network, or list of networks, supplied as matrices.
#' @param shared.spp Logical, whether to use networks of shared species only.
#' @param ... Arguments to \code{\link[dissCqN]{netShared}} (if \code{shared.spp
#'   = TRUE}).
#' @return A matrix with assemblages in rows and species interactions in
#'   columns.
#' @examples
intMat <- function(net, shared.spp = FALSE, ...) {

  n <- net; ss <- shared.spp

  ## List of networks
  if (ss) n <- netShared(n, ...)
  if (!is.list(n)) n <- list(n)
  if (is.null(names(n))) names(n) <- paste0("net", 1:length(n))

  ## List of interactions
  int <- lapply(n, function(i) {
    if (sum(i) > 0) {
      nam <- c(outer(rownames(i), colnames(i), paste, sep = ":"))
      setNames(c(i), nam)
    }
  })
  int.nam <- sort(unique(unlist(lapply(int, names))))

  ## Matrix of interactions
  m <- sapply(int, function(i) {
    if (!is.null(i)) {
      sapply(int.nam, function(j) {
        if (j %in% names(i)) i[[j]] else 0
      })
    } else {
      rep(0, max(1:length(int.nam)))
    }
  })
  m <- if (!is.matrix(m)) matrix(m) else t(m)
  dimnames(m) <- list(names(n), int.nam)
  m[, colSums(m) > 0, drop = FALSE]

}


#' @title Multiple Assemblage Dissimilarity
#' @description Multiple assemblage dissimilarity for orders \emph{q} =
#'   0-\emph{N}.
#' @param mat A matrix with assemblages in rows and species or species
#'   interactions in columns. Alternatively, a list of matrices, which will be
#'   interpreted as interaction networks and used to construct an assemblage x
#'   interaction matrix.
#' @param q Integer, the orders of \code{q} for which to calculate
#'   dissimilarity. Can be any set of integers between 0 and \emph{N} (the
#'   number of assemblages in \code{mat}).
#' @param shared.spp Logical, whether to compare networks of shared species only
#'   (if \code{mat} is a list of networks).
#' @details
#' @return A numeric vector of dissimilarities for the orders of \code{q}.
#' @references Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H., & Chazdon, R. L.
#'   (2008). A Two-Stage Probabilistic Approach to Multiple-Community Similarity
#'   Indices. \emph{Biometrics}, \strong{64}(4), 1178–1186.
#'   \url{https://doi.org/fcvn63}
#'
#'   Diserud, O. H., & Ødegaard, F. (2007). A multiple-site similarity measure.
#'   \emph{Biology Letters}, \strong{3}(1), 20–22. \url{https://doi.org/bwhfx6}
#'
#'   Jost, L., Chao, A., & Chazdon, R. L. (2011). Compositional similarity and
#'   beta diversity. In A. E. Magurran & B. J. McGill (Eds.), \emph{Biological
#'   Diversity: Frontiers in Measurement and Assessment} (pp. 66–84). Oxford
#'   University Press.
#' @examples
dissCqN <- function(mat, q = 0:2, shared.spp = FALSE) {

  m <- mat; ss <- shared.spp

  ## Matrix of assemblages x species/interactions
  m <- if (is.list(m)) intMat(m, ss, pairwise = FALSE) else as.matrix(m)
  if (any(is.na(m))) {
    warning("Rows with missing data (NA) removed.")
    m <- na.omit(m)
  }
  m <- m[rowSums(m) > 0, colSums(m) > 0, drop = FALSE]
  n <- nrow(m)

  ## Orders of dissimilarity
  q <- as.integer(round(q[q >= 0 & q <= max(2, n)]))
  names(q) <- paste0("C", q, "N")

  ## Species' relative abundances
  if (any(q > 0)) {
    m <- sweep(m, 1, rowSums(m), "/")
    m[is.nan(m)] <- 0
  }

  ## CqN dissimilarity
  1 - sapply(q, function(q) {
    if (n > 1) {
      if (q == 0) {
        s <- ncol(m)
        sm <- mean(rowSums(m > 0))
        (n - s / sm) / (n - 1)
      }
      else if (q == 1) {
        1 / log(n) *
          sum(apply(m, 2, function(i) {
            i <- i[i > 0]
            sum(sapply(1:length(i), function(j) {
              i[j] / n * log(1 + sum(i[-j] / i[j]))
            }))
          }))
      }
      else {
        sp <- colSums(m)
        spq <- colSums(m^q)
        (1 / (n^q - n) * sum(sp^q - spq)) /
          (1 / n * sum(spq))
      }
    }
    else NA
  })

}


#' @title Pairwise Dissimilarity Matrix
#' @description Generate pairwise dissimilarity matrix from a set of species
#'   assemblages.
#' @param mat A matrix with assemblages in rows and species or species
#'   interactions in columns. Alternatively, a list of matrices, which will be
#'   interpreted as interaction networks and used to construct an assemblage x
#'   interaction matrix.
#' @param q Integer, the order of \code{q} for which to calculate dissimilarity
#'   (passed to \code{\link[dissCqN]{dissCqN}}). Should be a single integer
#'   value, and if not, only the first value will be used.
#' @param compare.sub Subsets of assemblages to compare pairwise. These should
#'   be supplied as a list of two sets of assemblage names or indices. If only
#'   one set is supplied, this set is compared to all other assemblages in
#'   \code{mat}. If more than two sets are supplied, only the first two are
#'   used. If \code{NULL}, all assemblages are compared.
#' @param shared.spp Logical, whether to compare networks of shared species only
#'   (if \code{mat} is a list of networks).
#' @param parallel The type of parallel processing to use, if any. Can be one of
#'   \code{"snow"}, \code{"multicore"}, or \code{"no"} (for none - the default).
#' @param ncpus Number of system cores to use for parallel processing. If
#'   \code{NULL} (default), all available cores are used.
#' @param cl Optional cluster to use if \code{parallel = "snow"}. If \code{NULL}
#'   (default), a local cluster is created using the specified number of cores.
#' @return A matrix of pairwise assemblage dissimilarities.
#' @examples
dissMat <- function(mat, q = 0, compare.sub = NULL, shared.spp = FALSE,
                    parallel = "no", ncpus = NULL, cl = NULL) {

  m <- mat; cs <- compare.sub; ss <- shared.spp; p <- parallel; nc <- ncpus

  ## Networks?
  net <- is.list(m)
  net2 <- net && ss

  ## Matrix of assemblages x species/interactions
  m <- if (net) intMat(m, ss, compare.sub = cs) else as.matrix(m)

  ## Pairwise comparisons (indices, names)
  s <- 1:nrow(m)
  s1 <- if (net2) s[c(FALSE, TRUE)] else s
  s2 <- if (net2) s[c(TRUE, FALSE)] else s
  n <- if (net) names(mat) else rownames(m)
  if (is.null(n)) n <- paste0(if (net) "network" else "assemblage", s)
  n2 <- n1 <- n

  ## Compare subsets?
  if (!is.null(cs)) {
    if (!is.list(cs)) cs <- list(cs)
    n1 <- cs[[1]]
    if (is.numeric(n1)) n1 <- n[n1]
    n2 <- if (length(cs) > 1) cs[[2]]
    if (is.numeric(n2)) n2 <- n[n2]
    if (is.null(n2)) n2 <- n[!n %in% n1]
    if (!net2) {
      s1 <- which(n %in% n1)
      s2 <- which(n %in% n2)
    }
  }

  ## Dissimilarity matrix
  d <- t(semEff::pSapply(s1, function(i) {
    sapply(s2, function(j) {
      m <- m[c(i, j), , drop = FALSE]
      dissCqN(m, q[1])
    })
  }, p, nc, cl))
  if (nrow(d) != length(s1)) d <- t(d)
  if (net2) {
    ind <- rep(n1, each = length(n2))
    d <- tapply(diag(d), ind, eval, simplify = FALSE)
    d <- do.call(rbind, d)
  }
  dimnames(d) <- list(n1, n2)
  d

}


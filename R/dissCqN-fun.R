

#' @title Networks of Shared Species
#' @description Generate interaction networks comprising only the shared species
#'   across two or more networks.
#' @param net A list of two or more networks to compare, supplied as matrices.
#' @param pairwise Logical, whether to compare networks pairwise (default),
#'   rather than considering species shared across multiple networks.
#' @param compare.sub Subsets of networks to compare pairwise. These should be
#'   supplied as a list of two sets of network names or indices. If only one set
#'   is supplied, this is compared to all other networks in `net`. If more than
#'   two sets are supplied, only the first two are used.
#' @note If comparing networks pairwise, and subsets are not specified, the
#'   output will contain network self-comparisons (redundancy).
#' @return A list of networks of shared species. If comparing pairwise, this
#'   will be of length *n1* * *n2* * 2 (with *n1* and *n2* being the numbers of
#'   networks in each set), or if considering multiple networks, the length of
#'   the original list.
#' @export
netShared <- function(net, pairwise = TRUE, compare.sub = NULL) {

  n <- net; cs <- compare.sub

  # List of networks
  if (!isList(n)) n <- list(n)
  if (is.null(names(n))) names(n) <- paste0("network", 1:length(n))
  nn <- names(n)

  # Networks of shared species
  ns <- if (length(n) > 1) {

    if (pairwise) {

      # Networks to compare pairwise
      n2 <- n1 <- n

      # Compare subsets?
      if (!is.null(cs)) {
        if (!isList(cs)) cs <- list(cs)
        s1 <- cs[[1]]
        n1 <- if (is.character(s1)) n[which(nn %in% s1)] else n[s1]
        s2 <- if (length(cs) == 1) which(!nn %in% names(n1)) else cs[[2]]
        n2 <- if (is.character(s2)) n[which(nn %in% s2)] else n[s2]
      }

      # Networks of shared species
      ns <- lapply(n1, function(i) {
        lapply(n2, function(j) {
          ss <- Map(function(k, l) k[k %in% l], dimnames(i), dimnames(j))
          lapply(list(i, j), function(k) {
            k[ss[[1]], ss[[2]], drop = FALSE]
          })
        })
      })
      ns <- lapply(rapply(ns, enquote, how = "unlist"), eval)
      nsn <- gsub("(.{1}$)", "_\\1", names(ns))
      setNames(ns, nsn)

    } else {

      # Shared species across all networks
      spp <- unique(unlist(lapply(n, dimnames)))
      ss <- spp[sapply(spp, function(i) {
        all(sapply(n, function(j) i %in% unlist(dimnames(j))))
      })]

      # Networks of shared species
      lapply(n, function(i) {
        s1 <- rownames(i) %in% ss
        s2 <- colnames(i) %in% ss
        i[s1, s2, drop = FALSE]
      })

    }

  } else n
  return(ns)

}


#' @title Assemblage x Species Interaction Matrix
#' @description Generate a matrix of assemblages x species interactions from a
#'   set of networks.
#' @param net An interaction network, or list of networks, supplied as matrices.
#' @param shared.spp Logical, whether to use networks of shared species only.
#' @param ... Arguments to [netShared()] (if `shared.spp = TRUE`).
#' @return A matrix with assemblages in rows and species interactions in
#'   columns.
#' @export
intMat <- function(net, shared.spp = FALSE, ...) {

  n <- net

  # List of networks
  if (shared.spp) n <- netShared(n, ...)
  if (!isList(n)) n <- list(n)
  if (is.null(names(n))) names(n) <- paste0("network", 1:length(n))

  # List of interactions
  int <- lapply(n, function(i) {
    if (sum(i) > 0) {
      nam <- c(outer(rownames(i), colnames(i), paste, sep = ":"))
      setNames(c(i), nam)
    }
  })
  int.nam <- sort(unique(unlist(lapply(int, names))))

  # Matrix of interactions
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
  m <- m[, colSums(m) > 0, drop = FALSE]
  return(m)

}


#' @title Multiple Assemblage Dissimilarity
#' @description Multiple assemblage dissimilarity for orders *q* = 0-*N*.
#' @param mat A matrix with assemblages in rows and species or species
#'   interactions in columns. Alternatively, a list of matrices, which will be
#'   interpreted as interaction networks and used to construct an assemblage x
#'   interaction matrix.
#' @param q Integer, the order(s) of *q* for which to calculate dissimilarity.
#'   Can be any set of integers between 0 and *N* (the number of assemblages in
#'   `mat`).
#' @param pairwise Logical, whether to calculate pairwise, rather than multiple
#'   assemblage, dissimilarity.
#' @param compare.sub Subsets of assemblages to compare pairwise. These should
#'   be supplied as a list of two sets of assemblage names or indices. If only
#'   one set is supplied, this is compared to all other assemblages in `mat`. If
#'   more than two sets are supplied, only the first two are used. If `NULL`
#'   (default), all assemblages are compared.
#' @param shared.spp Logical, whether to compare networks of shared species only
#'   (if `mat` is a list of networks).
#' @param parallel The type of parallel processing to use, if any. Can be one of
#'   `"snow"`, `"multicore"`, or `"no"` (for none – the default). Passed to
#'   [pSapply()].
#' @param ncpus Number of system cores to use for parallel processing. If `NULL`
#'   (default), all available cores are used.
#' @param cl Optional cluster to use if `parallel = "snow"`. If `NULL`
#'   (default), a local cluster is created using the specified number of cores.
#' @details Dissimilarity is calculated here for multiple species assemblages
#'   (or interaction networks) via the *CqN* generalisation of similarity
#'   indices (Chao *et al.* 2008, Jost *et al.* 2011). Increasing the value of
#'   `q` increases the 'depth' of the measure, that is, how much emphasis is
#'   placed on changes in relative abundance of the most common species. Setting
#'   `q = 0` represents the qualitative Sørensen index (Sørensen 1948), where
#'   rare and common species are treated equally. `q` > 0 is more sensitive to
#'   common species, with `q = 1` representing the Shannon-based Horn index
#'   (Horn 1966) and `q = 2` the Simpson-based Morisita-Horn index (Morisita
#'   1959, Horn 1966). For *N* > 2, indices are generalised to consider species
#'   shared across multiple assemblages (Diserud & Ødegaard 2007, eqns. 6.3-6.5
#'   in Jost *et al.* 2011). For `q` >= 2 <= *N*, common species increasingly
#'   dominate the measure, and it can then be interpreted as the ratio of two
#'   probabilities of randomly sampling `q` individuals of the same species from
#'   the *N* assemblages, where 1) the individuals came from at least one
#'   different assemblage (\eqn{^{q}G_{D}}{qGD}) and 2) they all came from the
#'   same assemblage (\eqn{^{q}G_{S}}{qGS}) (Jost *et al.* 2011). Dissimilarity
#'   is thus:
#'
#'   \deqn{1 - ^{q}G_{D} / ^{q}G_{S}}{1 - qGD / qGS}
#'
#'   Pairwise dissimilarity can be calculated for all or a subset of the
#'   assemblages (or networks) in `mat`, in which case a dissimilarity matrix is
#'   returned (one for each value of `q`). If comparing subsets, the names or
#'   indices of assemblages to compare should be supplied to `compare.sub`. Note
#'   that pairwise calculation may take a long time if *N* is large, in which
#'   case parallel processing may speed up results (e.g. `parallel = "snow"`).
#'
#'   If `shared.spp = TRUE` and `mat` is a list of interaction networks (as
#'   matrices), multiple or pairwise interaction dissimilarity will be
#'   calculated for networks of shared species only (see [netShared()]). This
#'   can be useful to help partition the different components of network
#'   dissimilarity, e.g. dissimilarity due to interaction 'rewiring' among
#'   shared species vs. that due to species turnover (Poisot *et al.* 2012).
#' @return A numeric vector of dissimilarities, or a pairwise dissimilarity
#'   matrix (or list of matrices), for the orders of `q`.
#' @references Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H., & Chazdon, R. L.
#'   (2008). A Two-Stage Probabilistic Approach to Multiple-Community Similarity
#'   Indices. *Biometrics*, **64**(4), 1178–1186. \doi{10/fcvn63}
#'
#'   Diserud, O. H., & Ødegaard, F. (2007). A multiple-site similarity measure.
#'   *Biology Letters*, **3**(1), 20–22. \doi{10/bwhfx6}
#'
#'   Horn, H. S. (1966). Measurement of “Overlap” in Comparative Ecological
#'   Studies. *The American Naturalist*, **100**(914), 419–424.
#'   \doi{10/b62ct5}
#'
#'   Jost, L., Chao, A., & Chazdon, R. L. (2011). Compositional similarity and
#'   beta diversity. In A. E. Magurran & B. J. McGill (Eds.), *Biological
#'   Diversity: Frontiers in Measurement and Assessment* (pp. 66–84). Oxford
#'   University Press.
#'
#'   Morisita, M. (1959). Measuring of interspecific association and similarity
#'   between communities. *Memoirs of the Faculty of Science, Kyushu Univ.,
#'   Series E (Biology)*, **3**, 65–80.
#'
#'   Poisot, T., Canard, E., Mouillot, D., Mouquet, N., & Gravel, D. (2012). The
#'   dissimilarity of species interaction networks. *Ecology Letters*,
#'   **15**(12), 1353–1361. \doi{10/f4dv37}
#'
#'   Sørensen, T. (1948). A method of establishing groups of equal amplitude in
#'   plant sociology based on similarity of species and its application to
#'   analyses of the vegetation on Danish commons. *Kongelige Danske
#'   Videnskabernes Selskabs Biologiske Skrifter*, **5**, 1–34.
#' @examples
#' # Sample community data from SpadeR package (three assemblages, 120 species)
#' data(SimilarityMultData, package = "SpadeR")
#' d <- SimilarityMultData$Abu
#'
#' # Multiple-assemblage dissimilarity for q = 0:2
#' (CqN <- dissCqN::dissCqN(t(d)))
#'
#' # Compare to empirical CqN values from SpadeR::SimilarityMult()
#' sim <- SpadeR::SimilarityMult(d, datatype = "abundance", nboot = 1)
#' CqN_2 <- 1 - c(
#'   "C0N" = sim$Empirical_richness["C0N(q=0,Sorensen)", "Estimate"],
#'   "C1N" = sim$Empirical_relative["C1N=U1N(q=1,Horn)", "Estimate"],
#'   "C2N" = sim$Empirical_relative["C2N(q=2,Morisita)", "Estimate"]
#' )
#' stopifnot(all.equal(CqN, CqN_2))
#'
#' # Pairwise dissimilarity matrices
#' dissCqN::dissCqN(t(d), pairwise = TRUE)
#' @export
dissCqN <- function(mat, q = 0:2, pairwise = FALSE, compare.sub = NULL,
                    shared.spp = FALSE, parallel = "no", ncpus = NULL,
                    cl = NULL) {

  m <- mat; cs <- if (pairwise) compare.sub; ss <- shared.spp; nc <- ncpus

  # Interaction networks?
  net <- isList(m)

  # Matrix of assemblages x species/interactions
  if (net) m <- intMat(m, ss, pairwise = pairwise, compare.sub = cs)
  m <- as.matrix(m)
  if (any(is.na(m))) {
    warning("Assemblages with missing data (NA) removed.")
    m <- na.omit(m)
  }

  # No. of assemblages
  N <- if (!pairwise) sum(rowSums(m) > 0) else 2

  # Orders of dissimilarity
  q <- as.integer(round(q[q >= 0 & q <= max(2, N)]))
  names(q) <- paste0("C", q, "N")

  # Species' relative abundances/frequencies
  if (any(q > 0)) {
    m <- sweep(m, 1, rowSums(m), "/")
    m[is.nan(m)] <- 0
  }

  # Function to calculate multiple assemblage dissimilarity (CqN)
  dissCqN <- function(q, m) {

    # No. of assemblages
    m <- m[rowSums(m) > 0, , drop = FALSE]
    N <- nrow(m)

    # CqN dissimilarity
    1 - if (N > 1) {
      if (q == 0) {
        S <- ncol(m)
        Sm <- mean(rowSums(m > 0))
        (N - S / Sm) / (N - 1)
      }
      else if (q == 1) {
        1 / log(N) *
          sum(apply(m, 2, function(i) {
            pi <- i[i > 0]
            sum(sapply(1:length(pi), function(j) {
              pij <- pi[j]
              pik <- pi[-j]
              pij / N * log(1 + sum(pik / pij))
            }))
          }))
      }
      else {
        qGD <- 1 / (N^q - N) * sum(colSums(m)^q - colSums(m^q))
        qGS <- 1 / N * sum(m^q)
        qGD / qGS
      }
    }
    else NA

  }

  # Calculate multiple or pairwise dissimilarity
  if (pairwise) {

    # Networks of shared species?
    net2 <- net && ss

    # Pairwise comparisons (indices, names)
    s2 <- s1 <- s <- 1:nrow(m)
    if (net2) {
      s1 <- s[c(FALSE, TRUE)]
      s2 <- s[c(TRUE, FALSE)]
    }
    n <- if (net2) names(mat) else rownames(m)
    if (is.null(n)) n <- paste0(if (net) "network" else "assemblage", s)
    n2 <- n1 <- n

    # Compare subsets?
    if (!is.null(cs)) {
      if (!isList(cs)) cs <- list(cs)
      n1 <- cs[[1]]
      n1 <- if (is.character(n1)) n[n %in% n1] else n[n1]
      n2 <- if (length(cs) == 1)  n[!n %in% n1] else cs[[2]]
      n2 <- if (is.character(n2)) n[n %in% n2] else n[n2]
      if (!net2) {
        s1 <- which(n %in% n1)
        s2 <- which(n %in% n2)
      }
    }

    # Function to calculate pairwise dissimilarity matrix
    dissMat <- function(q) {
      d <- t(pSapply(s1, function(i) {
        sapply(s2, function(j) {
          m <- m[c(i, j), , drop = FALSE]
          m <- m[, colSums(m) > 0, drop = FALSE]
          dissCqN(q, m)
        })
      }, parallel, nc, cl))
      if (nrow(d) != length(s1)) d <- t(d)
      if (net2) {
        ind <- rep(n1, each = length(n2))
        d <- tapply(diag(d), ind, eval, simplify = FALSE)
        d <- do.call(rbind, d)
      }
      dimnames(d) <- list(n1, n2)
      return(d)
    }

    # Dissimilarity matrix or matrices
    d <- sapply(q, dissMat, simplify = FALSE)
    if (isList(d) && length(d) < 2) d <- d[[1]]
    return(d)

  }
  else sapply(q, dissCqN, m)

}


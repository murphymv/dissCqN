

#' @keywords internal
#' @title Object Types
#' @description Functions to determine the 'type' of an R object using classes.
#'   Intended largely for convenience and internal use.
#' @param x An R object.
#' @return A logical value.
#' @name Object.Type
NULL
#' @describeIn Object.Type Is object a list (class `"list"` only)?
isList <- function(x) {
  class(x)[1] == "list"
}
#' @describeIn Object.Type Is object a boot object (class `"boot"`)?
isBoot <- function(x) {
  "boot" %in% class(x)
}
#' @describeIn Object.Type Is object a fitted model?
isMod <- function(x) {
  any(c("lm", "glm", "lmerMod", "glmerMod", "lmerModLmerTest", "gls", "betareg")
      %in% class(x))
}


#' @keywords internal
#' @title Recursive [mapply()]
#' @description Recursively apply a function to a list or lists.
#' @param FUN Function to apply.
#' @param ... Object(s) to which `FUN` can be applied, or lists of such objects
#'   to iterate over (defined narrowly, as of class `"list"`).
#' @param MoreArgs A list of additional arguments to `FUN`.
#' @param SIMPLIFY Logical, whether to simplify the results to a vector or
#'   array.
#' @param USE.NAMES Logical, whether to use the names of the first list object
#'   in `...` for the output.
#' @details `rMapply()` recursively applies `FUN` to the elements of the lists
#'   in `...` via [mapply()]. If only a single list is supplied, the function
#'   acts like a recursive version of [sapply()]. The particular condition that
#'   determines if the function should stop recursing is if either the first or
#'   second objects in `...` are not of class `"list"`. Thus, unlike [mapply()],
#'   it will not iterate over non-list elements in these objects, but instead
#'   returns the output of `FUN(...)`.
#'
#'   This is primarily a convenience function used internally to enable
#'   recursive application of functions to lists or nested lists. Its particular
#'   stop condition for recursing is also designed to either *a)* act as a
#'   wrapper for `FUN` if the first object in `...` is not a list, or *b)* apply
#'   a weighted averaging operation if the first object is a list and the second
#'   object is a numeric vector of weights.
#' @return The output of `FUN` in a list or nested list, or simplified to a
#'   vector or array (or list of arrays).
rMapply <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE,
                    USE.NAMES = TRUE) {
  l <- list(...)
  n <- length(l)
  i <- if (n > 0) l[[1]] else l
  j <- if (n > 1) l[[2]] else i
  if (!isList(i) || !isList(j)) {
    do.call(FUN, c(l, MoreArgs))
  } else {
    a <- list(FUN = FUN, MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY,
              USE.NAMES = USE.NAMES)
    mapply(rMapply, ..., MoreArgs = a, SIMPLIFY = SIMPLIFY,
           USE.NAMES = USE.NAMES)
  }
}


#' @keywords internal
#' @title Parallel [sapply()]
#' @description Apply a function to a vector using parallel processing.
#' @param X A vector object (numeric, character, or list).
#' @param FUN Function to apply to the elements of `X`.
#' @param parallel The type of parallel processing to use. Can be one of
#'   `"snow"` (default), `"multicore"` (not available on Windows), or `"no"`
#'   (for none). See Details.
#' @param ncpus Number of system cores to use for parallel processing. If `NULL`
#'   (default), all available cores are used.
#' @param cl Optional cluster to use if `parallel = "snow"`. If `NULL`
#'   (default), a local cluster is created using the specified number of cores.
#' @param add.obj A character vector of any additional object names to be
#'   exported to the cluster. Use if a required object or function cannot be
#'   found.
#' @param ... Additional arguments to [parSapply()],
#'   [`mcmapply()`](https://rdrr.io/r/parallel/unix/mclapply.html), or
#'   [sapply()] (note: arguments `"simplify"` and `"SIMPLIFY"` are both
#'   allowed).
#' @details This is a wrapper for [parallel::parSapply()] (`"snow"`) or
#'   [`parallel::mcmapply()`](https://rdrr.io/r/parallel/unix/mclapply.html)
#'   (`"multicore"`), enabling (potentially) faster processing of a function
#'   over a vector of objects. If `parallel = "no"`, [sapply()] is used instead.
#'
#'   Parallel processing via option `"snow"` (default) is carried out using a
#'   cluster of workers, which is automatically set up via [makeCluster()] using
#'   all available system cores or a user supplied number of cores. The function
#'   then exports the required objects and functions to this cluster using
#'   [clusterExport()], after performing a (rough) match of all objects and
#'   functions in the current global environment to those referenced in the call
#'   to `FUN` (and also any calls in `X`). Any additional required object names
#'   can be supplied using `add.obj`.
#' @return The output of `FUN` in a list, or simplified to a vector or array.
pSapply <- function(X, FUN, parallel = c("snow", "multicore", "no"),
                    ncpus = NULL, cl = NULL, add.obj = NULL, ...) {

  parallel <- match.arg(parallel); nc <- ncpus; ao <- add.obj; a <- list(...)
  if (parallel == "multicore") a$simplify <- NULL else a$SIMPLIFY <- NULL

  if (parallel != "no") {

    # No. cores to use
    if (is.null(nc)) nc <- parallel::detectCores()

    if (parallel == "snow") {

      # Create local cluster using system cores
      if (is.null(cl)) {
        cl <- parallel::makeCluster(getOption("cl.cores", nc))
      }

      # Export required objects/functions to cluster
      # (search global env. for objects in calls to X/FUN)
      P <- function(...) {
        paste(..., collapse = " ")
      }
      xc <- P(unlist(rMapply(function(i) {
        if (isMod(i) || isBoot(i)) P(getCall(i))
      }, X)))
      fa <- P(sapply(match.call(expand.dots = FALSE)$..., deparse))
      fc <- P(xc, enquote(FUN)[2], fa)
      o <- unlist(lapply(search(), ls))
      o <- o[sapply(o, function(i) grepl(i, fc, fixed = TRUE))]
      o <- c("X", o, ao)
      parallel::clusterExport(cl, o, environment())

      # Run parSapply using cluster
      out <- do.call(parallel::parSapply, c(list(cl, X, FUN), a))
      parallel::stopCluster(cl)

    } else {
      out <- parallel::mcmapply(FUN, X, mc.cores = nc, MoreArgs = a)
    }

  } else {
    out <- do.call(sapply, c(list(X, FUN), a))
  }

  return(out)

}


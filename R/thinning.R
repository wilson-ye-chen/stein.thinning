#' Perform Stein thinning
#'
#' Optimally select m points from n > m samples generated from a target
#' distribution of d dimensions.
#'
#' @param smp n-by-d matrix of n d-dimensional sample points.
#' @param scr n-by-d matrix of gradients of log target at \code{smp}.
#' @param stnd optional logical, either TRUE (default) or FALSE, indicating
#'        whether or not to standardise the colums of \code{smp}.
#' @param pre optional string, either "id" (default), "med", "sclmed", or
#'        "smpcov", specifying the preconditioner to be used. Alternatively,
#'        a numeric string can be passed as the single length-scale parameter
#'        of an isotropic kernel.
#' @return m-vector containing the row indices in \code{smp} (and \code{scr})
#'         of the selected points.
#' @export
thin <- function(smp, scr, m, stnd=TRUE, pre="id") {
    # Argument checks
    if (length(dim(smp)) != 2 || length(dim(scr)) != 2) {
        stop("smp or scr is not two-dimensional.")
    }
    n <- nrow(smp)
    d <- ncol(smp)
    if (n == 0 || d == 0) {
        stop("smp is empty.")
    }
    if (nrow(scr) != n || ncol(scr) != d) {
        stop("Dimensions of smp and scr are inconsistent.")
    }
    if (sum(is.nan(smp)) + sum(is.nan(scr)) > 0) {
        stop("smp or scr contains NaNs.")
    }
    if (sum(is.infinite(smp)) + sum(is.infinite(scr)) > 0) {
        stop("smp or scr contains Infs.")
    }

    # Standardisation
    if (stnd) {
        loc <- colMeans(smp)
        scl <- colMeans(abs(smp - matrix(rep(loc, each=n), nrow=n)))
        if (min(scl) == 0) {
            stop("Too few unique samples in smp.")
        }
        scl <- matrix(rep(scl, each=n), nrow=n)
        smp <- smp / scl
        scr <- scr * scl
    }

    # Vectorised Stein kernel function
    vfk0 <- make_imq(smp, scr, pre)

    # Pre-allocate arrays
    k0 <- matrix(NA, nrow=n, ncol=m)
    idx <- c(NA)
    length(idx) <- m

    # Populate columns of k0 as new points are selected
    k0[,1] <- vfk0(smp, smp, scr, scr)
    idx[1] <- which.min(k0[,1])
    message("THIN: 1 of ", m)
    v1 <- matrix(1, n, 1)
    for (i in 2:m) {
        smp_last <- kronecker(v1, smp[idx[i - 1],, drop=FALSE])
        scr_last <- kronecker(v1, scr[idx[i - 1],, drop=FALSE])
        k0[,i] <- vfk0(smp, smp_last, scr, scr_last)
        idx[i] <- which.min(k0[,1] + 2 * rowSums(k0[,2:i, drop=FALSE]))
        message("THIN: ", i, " of ", m)
    }
    return(idx)
}

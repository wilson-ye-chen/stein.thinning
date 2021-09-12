#' Stein kernel function
#'
#' Vectorised Stein kernel function based on the inverse multi-quadric (IMQ)
#' kernel. Each argument can be m-by-d matrices containing m d-dimensional
#' points.
#'
#' @param a first argument of the kernel, m-by-d matrix.
#' @param b second argument of the kernel, m-by-d matrix.
#' @param sa gradients of the log target density at \code{a}, m-by-d matrix.
#' @param sb gradients of the log target density at \code{b}, m-by-d matrix.
#' @param linv inverse preconditioner, d-by-d matrix.
#' @return value of the Stein kernel, m-vector.
#' @export
vfk0_imq <- function(a, b, sa, sb, linv) {
    amb <- t(a) - t(b)
    qf <- 1 + colSums(linv %*% amb * amb)
    t1 <- -3 * colSums(linv %*% linv %*% amb * amb) / (qf ^ 2.5)
    t2 <- (sum(diag(linv)) + colSums(linv %*% (t(sa) - t(sb)) * amb)) / (qf ^ 1.5)
    t3 <- colSums(t(sa) * t(sb)) / (qf ^ 0.5)
    return(t1 + t2 + t3)
}

#' Create an inverse preconditioner
#'
#' @param smp n-by-d matrix of n d-dimensional sample points.
#' @param scr n-by-d matrix of gradients of log target at \code{smp}.
#' @param pre optional string, either "id" (default), "med", "sclmed", or
#'        "smpcov", specifying the preconditioner to be used. Alternatively,
#'        a numeric string can be passed as the single length-scale parameter
#'        of an isotropic kernel.
#' @return inverse preconditioner, d-by-d matrix.
#' @export
make_precon <- function(smp, scr, pre="id") {
    sz = nrow(smp)
    dm = ncol(smp)

    med2 <- function(m) {
        # Returns the squared pairwise median
        if (sz > m) {
            sub <- smp[as.integer(seq(1, sz, length=m)),]
        } else {
            sub <- smp
        }
        d <- pdist(sub, sub)
        d <- d[lower.tri(d)]
        return(median(d) ^ 2)
    }

    # Select preconditioner
    m <- 1000
    if (pre == "id") {
        linv <- diag(dm)
    } else if (pre == "med") {
        m2 <- med2(m)
        if (m2 == 0) {
            stop("Too few unique samples in smp.")
        }
        linv <- solve(m2 * diag(dm))
    } else if (pre == "sclmed") {
        m2 <- med2(m)
        if (m2 == 0) {
            stop("Too few unique samples in smp.")
        }
        linv <- solve(m2 / log(min(m, sz)) * diag(dm))
    } else if (pre == "smpcov") {
        c <- cov(smp)
        v <- eigen(c, symmetric=T, only.values=T)$values
        if (!all(v > 0)) {
            stop("Too few unique samples in smp.")
        }
        linv <- solve(c)
    } else if (!is.na(as.double(pre))) {
        linv <- solve(as.double(pre) * diag(dm))
    } else {
        stop("incorrect preconditioner type.")
    }
    return(linv)
}

#' Create the (IMQ) Stein kernel function with a preconditioner
#'
#' @param smp n-by-d matrix of n d-dimensional sample points.
#' @param scr n-by-d matrix of gradients of log target at \code{smp}.
#' @param pre optional string, either "id" (default), "med", "sclmed", or
#'        "smpcov", specifying the preconditioner to be used. Alternatively,
#'        a numeric string can be passed as the single length-scale parameter
#'        of an isotropic kernel.
#' @return Stein kernel function.
#' @export
make_imq <- function(smp, scr, pre="id") {
    linv <- make_precon(smp, scr, pre)
    vfk0 <- function(a, b, sa, sb) {
        return(vfk0_imq(a, b, sa, sb, linv))
    }
    return(vfk0)
}

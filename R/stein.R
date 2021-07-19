#' Compute a cumulative sequence of KSD values
#'
#' @param x n-by-d matrix of n d-dimensional sample points.
#' @param s n-by-d matrix of gradients of log target at \code{x}.
#' @param vfk0 vectorised Stein kernel function.
#' @return n-vector of cumulative KSD values.
#' @export
ksd <- function(x, s, vfk0) {
    n <- nrow(x)
    ks <- c(NA)
    length(ks) <- n
    ps <- 0
    for (i in 1:n) {
        v1 <- matrix(1, i, 1)
        x_i <- kronecker(v1, x[i,, drop=F])
        s_i <- kronecker(v1, s[i,, drop=F])
        k0 <- vfk0(x_i, x[1:i,, drop=F], s_i, s[1:i,, drop=F])
        ps <- ps + 2 * sum(k0[1:(i - 1)]) + k0[i]
        ks[i] <- sqrt(ps) / i
        message("KSD: ", i, " of ", n)
    }
    return(ks)
}

#' Compute a Stein kernel matrix
#'
#' @param x n-by-d matrix of n d-dimensional sample points.
#' @param s n-by-d matrix of gradients of log target at \code{x}.
#' @param vfk0 vectorised Stein kernel function.
#' @return n-by-n Stein kernel matrix.
#' @export
kmat <- function(x, s, vfk0) {
    n <- nrow(x)
    k0 <- matrix(NA, nrow=n, ncol=n)
    for (i in 1:n) {
        for (j in 1:i) {
            k0[i, j] <- vfk0(
                x[i,, drop=F],
                x[j,, drop=F],
                s[i,, drop=F],
                s[j,, drop=F])
        }
    }
    k0 <- t(k0)
    lt <- lower.tri(k0)
    k0[lt] <- t(k0)[lt]
    return(k0)
}

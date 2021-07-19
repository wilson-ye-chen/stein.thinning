pdist <- function(a, b) {
    # Adopted from Felix Riedel's blog post:
    # https://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
    an = apply(a, 1, function(rvec) crossprod(rvec,rvec))
    bn = apply(b, 1, function(rvec) crossprod(rvec,rvec))
    m = nrow(a)
    n = nrow(b)
    tmp = matrix(rep(an, n), nrow=m)
    tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
    return(sqrt(tmp - 2 * tcrossprod(a, b)))
}

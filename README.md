# Stein Thinning for R
This R package implements an algorithm for optimally compressing
sampling algorithm outputs by minimising a kernel Stein discrepancy.
Please see the accompanying paper "Optimal Thinning of MCMC Output"
([arXiv](https://arxiv.org/pdf/2005.03952.pdf)) for details of the
algorithm.

# Installing via Github
One can install the package directly from this repository:
```r
install.packages("devtools")
devtools::install_github("wilson-ye-chen/stein.thinning")
```
The first line above is not needed if you have `devtools` installed.

# Getting Started
For example, correlated samples from a posterior distribution are
obtained using a MCMC algorithm and stored in the matrix `smpl`,
and the corresponding gradients of the log-posterior are stored in
another matrix `grad`. One can then perform Stein Thinning to
obtain a subset of 40 sample points by running the following code:
```r
idx <- thin(smpl, grad, 40)
```
The `thin` function returns a vector containing the row indices in
`smpl` (and `grad`) of the selected points. Please refer to `demo.R`
as a starting example. To run the demo:
```r
demo()
```

The default usage requires no additional user input and is based on
the identity (`id`) preconditioning matrix and standardised sample.
Alternatively, the user can choose to specify which heuristic to use
for computing the preconditioning matrix by setting the option string
to either `id`, `med`,  `sclmed`, or `smpcov`. Standardisation can be
disabled by setting `stnd=FALSE`. For example, the default setting
corresponds to:
```r
idx <- thin(smpl, grad, 40, stnd=TRUE, pre='id')
```
The details for each of the heuristics are documented in Section 2.3 of
the accompanying paper.

# RStan Example
As an illustration of how Stein Thinning can be used to post-process
output from [Stan](https://mc-stan.org/rstan/), consider the following
simple Stan script that produces correlated samples from a bivariate
Gaussian model:
```r
mc <- "
parameters {vector[2] x;}
model {x ~ multi_normal([0, 0], [[1, 0.8], [0.8, 1]]);}
"
fit <- rstan::stan(model_code=mc, iter=1000, chains=1)
```
The bivariate Gaussian model is used for illustration, but regardless of
the complexity of the model being sampled the output of Stan will always
be a `fit` object (of `stanfit` class). The sampled points and the
log-posterior gradients can be extracted from the returned `fit` object:
```r
smpl <- rstan::extract(fit, permuted=FALSE, inc_warmup=TRUE)
smpl <- smpl[,,1:2]
grad <- t(apply(smpl, 1, function(x) rstan::grad_log_prob(fit, x)))
idx <- thin(smpl, grad, 40)
```
The above example can be found in `demo.R`. To run the RStan example:
```r
demo_stan()
```

# Functions
* `thin(smp, scr, m, stnd=TRUE, pre="id")` returns indices of thinned points.
* `demo()` runs an example of post-processing MCMC output from CSV files.
* `demo_stan()` runs an example of post-processing Stan output.
* `ksd(x, s, vfk0)` returns cumulative KSD values of sample `x`.
* `kmat(x, s, vfk0)` returns a Stein kernel matrix of sample `x`.
* `make_imq(smp, scr, pre="id")` returns IMQ kernel with a predefined IPM.
* `make_precon(smp, scr, pre="id")` returns a predefined IPM.
* `vfk0_imq(a, b, sa, sb, linv)` evaluates IMQ kernel for any IPM.

Acronyms:
* IPM: inverse preconditioning matrix.
* IMQ: inverse multi-quadric.
* KSD: Kernelized Stein discrepancy.

#' Post-process MCMC output from CSV files
#' @export
demo <- function() {
    # Read MCMC output from files
    smp_file <- system.file("extdata", "smp.csv", package="stein.thinning")
    scr_file <- system.file("extdata", "scr.csv", package="stein.thinning")
    smp <- data.matrix(read.csv(smp_file, header=FALSE))
    scr <- data.matrix(read.csv(scr_file, header=FALSE))

    # Run Stein Thinning
    idx <- thin(smp, scr, 40)

    # Plot point-set over trace
    dev.new()
    plot(smp[,1], smp[,2], xlab="x_1", ylab="x_2",
        type="l", col=rgb(0, 0, 0, 0.4))
    points(smp[idx, 1], smp[idx, 2], pch=19, col="red")

    # Compute KSD
    vfk0 <- make_imq(smp, scr, pre="sclmed")
    ks_smp <- ksd(smp, scr, vfk0)
    ks_x <- ksd(smp[idx,], scr[idx,], vfk0)

    # Print out the inverse of the preconditioner
    print(make_precon(smp, scr, pre="sclmed"))

    # Visualise the Stein kernel matrix
    k0 <- kmat(smp[idx,], scr[idx,], vfk0)
    dev.new()
    image(t(k0[nrow(k0):1,]), main="Stein kernel matrix", axes=FALSE)

    # Plot KSD curves
    dev.new()
    plot(log(1:length(ks_smp)), log(ks_smp), xlab="log n", ylab="log KSD",
        type="l", col="black")
    lines(log(1:length(ks_x)), log(ks_x), col="red")
}

#' Post-process output from Stan
#' @export
demo_stan <- function() {
    # Simple bivariate Gaussian model
    mc <- "
    parameters {vector[2] x;}
    model {x ~ multi_normal([0, 0], [[1, 0.8], [0.8, 1]]);}
    "
    fit <- rstan::stan(model_code=mc, iter=1000, chains=1)

    # Extract sampled points and gradients
    smp <- rstan::extract(fit, permuted=FALSE, inc_warmup=TRUE)
    smp <- smp[,,1:2]
    scr <- t(apply(smp, 1, function(x) rstan::grad_log_prob(fit, x)))

    # Obtain a subset of 40 points
    idx <- thin(smp, scr, 40)

    # Plot point-set over trace
    dev.new()
    plot(smp[,1], smp[,2], xlab="x_1", ylab="x_2",
        type="l", col=rgb(0, 0, 0, 0.4))
    points(smp[idx, 1], smp[idx, 2], pch=19, col="red")
}

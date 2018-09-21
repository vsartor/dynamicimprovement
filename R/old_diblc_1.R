# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#' Dynamic Improvement Bayesian Lee-Carter
#'
#' This routine is kept for archiving reasons and making old methodology
#' available. See `diblc` documentation for the latest method.
#'
#' This routine in particular performs full state expansion, making every
#' parameter of the model part of the state-space, all of them being modelled
#' with fixed discount factors, and with observational variances being estimated
#' with a Gibbs sampling step.
#'
#' This routines takes too long to run and unecessarily adds dynamics to some
#' of the parameters.
#'
#' @param Y Numerical matrix with age information on rows and time information
#'          on columns. Must be log-mortality with valid and finite values.
#' @param I Number of iterations for the Gibbs sampler.
#' @param B Number of iterations to be discarded.
#'
#' @return A `diblc` object.
#'
#' @importFrom assertthat assert_that
#' @importFrom MASS mvrnorm
#' @importFrom stats var rgamma rnorm
#'
#' @export
old_diblc_1 <- function(Y, I = 3000, B = 500)
{
    #-- Input validation --#

    assert_that(is.matrix(Y))
    assert_that(is.numeric(Y))
    N <- ncol(Y)
    X <- nrow(Y)
    assert_that(all(is.finite(Y)))
    if (any(Y >= 0)) {
        warning("Expected Y to be log-mortality. Non-negative values received.")
    }

    assert_that(is.numeric(I))
    assert_that(length(I) == 1)
    assert_that(I > 10)

    assert_that(is.numeric(B))
    assert_that(length(B) == 1)
    assert_that(0 < B & B < I)

    #-- Use empirical estimates to initialize phi --#

    emp.alpha <- apply(Y, 1, mean)
    emp.beta <- (Y[ ,1] - Y[ ,N]) / sum((Y[ ,1] - Y[ ,N])^2)
    emp.kappa <- apply(Y - emp.alpha, 2, function(d) mean(d/emp.beta))
    emp.res <- apply(rbind(Y, emp.kappa), 2, function(x) x[1:X] - emp.alpha - emp.beta*x[X+1])
    emp.phi <- 1 / apply(emp.res, 1, var)

    #-- Allocate chains --#

    # Precision
    ch.phi   <- alloc(X, I)

    # Filter components
    ch.beta  <- alloc(X, N, I)
    ch.alpha <- alloc(X, N, I)
    ch.kappa <- alloc(N, I)
    ch.delta <- alloc(N, I)

    #-- Initialize chains --#

    # Run filter using empirical phi
    mod <- old_diblc_filter_1(Y, phi = emp.phi)

    # Initialize with empirical phi and filter means
    ch.phi[ ,1]     <- emp.phi
    ch.alpha[ , ,1] <- mod$mean_a
    ch.beta[ , ,1]  <- mod$mean_b
    ch.kappa[ ,1]   <- mod$mean_k
    ch.delta[ ,1]   <- mod$mean_d

    #-- Gibbs algorithm --#

    for (i in 2:I) {
        if (i %% 30 == 0)
            catf("dyn_lc: %04d/%04d", i, I)

        # Phi step
        for (x in 1:X) {
            phi.shape   <- (N + 2e-5) / 2
            phi.rate    <- (sum((Y[x, ] - mod$f_mean[x, ])^2) + 2e-5)
            ch.phi[x,i] <- rgamma(1, phi.shape, phi.rate)
        }

        # Filter step
        mod <- old_diblc_filter_1(Y, phi = ch.phi[ ,i])

        for (n in 1:N) {
            # Alpha step
            ch.alpha[ ,n,i] <- mvrnorm(1, mod$mean_a[ ,n], mod$cov_a[ , ,n])

            # Beta step
            ch.beta[ ,n,i] <- mvrnorm(1, mod$mean_b[ ,n], mod$cov_b[ , ,n])
        }

        # Kappa step
        ch.kappa[ ,i] <- rnorm(N, mod$mean_k, mod$sd_k)

        # Delta step
        ch.delta[ ,i] <- rnorm(N, mod$mean_d, mod$sd_d)
    }

    #-- Wrap-up and return --#

    # Drop burnin
    drop.msk <- -(1:B)
    ch.phi   <- ch.phi[ ,drop.msk]
    ch.alpha <- ch.alpha[ , ,drop.msk]
    ch.beta  <- ch.beta[ , ,drop.msk]
    ch.kappa <- ch.kappa[ ,drop.msk]
    ch.delta <- ch.delta[ ,drop.msk]

    # Return
    ret <- list(alpha = ch.alpha, beta = ch.beta, kappa = ch.kappa, phi = ch.phi,
                delta = ch.delta, data = Y)
    class(ret) <- "diblc"
    return(ret)
}

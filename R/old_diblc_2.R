# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#' Dynamic Improvement Bayesian Lee-Carter
#'
#' This routine is kept for archiving reasons and making old methodology
#' available. See `diblc` documentation for the latest method.
#'
#' This routine in particular considers empirical estimates of `alpha` and
#' estimates individual evolutional variances for each of the state-space
#' parameters.
#'
#' This routines gives too much freedom for the evolutional variances of the
#' improvement terms.
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
old_diblc_2 <- function(Y, I = 3000, B = 500)
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
    emp.beta  <- (Y[ ,1] - Y[ ,N]) / sum((Y[ ,1] - Y[ ,N])^2)
    emp.kappa <- apply(Y - emp.alpha, 2, function(d) mean(d/emp.beta))
    emp.scale <- sqrt(sumabs2(emp.beta)) * mode(sign(emp.beta))
    emp.beta  <- emp.beta / emp.scale
    emp.kappa <- emp.kappa * emp.scale
    emp.delta <- (emp.kappa[N] - emp.kappa[1]) / (N - 1)
    emp.res   <- apply(rbind(Y, emp.kappa), 2, function(x) x[1:X] - emp.alpha - emp.beta*x[X+1])
    emp.phi   <- 1 / apply(emp.res, 1, var)
    emp.phi_k <- 1 / sumabs2(emp.kappa[2:N] - emp.delta - emp.kappa[1:(N-1)])

    #-- Allocate chains --#

    # Static components
    ch.phi   <- alloc(X, I)
    ch.delta <- alloc(I)
    ch.phi_k <- alloc(I)
    ch.phi_b <- alloc(X, I)

    # Filter components
    ch.beta  <- alloc(X, N, I)
    ch.kappa <- alloc(N, I)

    #-- Initialize chains --#

    # Run filter using empirical estimates for phi and delta
    mod <- old_diblc_filter_2(Y, emp.phi, emp.delta, emp.phi_k, rep(10, X))

    # Initialize statics with empirical values and dynamical with filter results
    ch.phi[ ,1]     <- emp.phi
    ch.delta[1]     <- emp.delta
    ch.phi_k[1]     <- emp.phi_k
    ch.phi_b[ ,1]   <- rep(10, X)
    ch.beta[ , ,1]  <- mod$mean_b
    ch.kappa[ ,1]   <- mod$mean_k

    #-- Gibbs algorithm --#

    start <- gets()

    for (i in 2:I) {
        if (i %% 30 == 0) {
            ips = (gets() - start) / i
            etr <- (I - i) * ips
            emr <- floor(etr / 60)
            catf("[%04d/%04d] ETR = %.0f:%02d", i, I, emr, floor(etr - emr * 60))
        }

        # Phi step
        for (x in 1:X) {
            phi.shape   <- (N + 2e-5) / 2
            phi.rate    <- (sum((Y[x, ] - mod$f_mean[x, ])^2) + 2e-5)
            ch.phi[x,i] <- rgamma(1, phi.shape, phi.rate)
        }

        # Filter step
        mod <- old_diblc_filter_2(Y, ch.phi[ ,i], ch.delta[i-1], ch.phi_k[i-1], ch.phi_b[ ,i-1])

        # Beta step
        for (n in 1:N) {
            ch.beta[ ,n,i] <- mvrnorm(1, mod$mean_b[ ,n], mod$cov_b[ , ,n])
        }

        # Kappa step
        ch.kappa[ ,i] <- rnorm(N, mod$mean_k, mod$sd_k)

        # Delta step
        delta.mu <- (ch.kappa[N,i] - ch.kappa[1,i]) / (N - 1)
        delta.sigma <- 1 / sqrt(ch.phi_k[i-1] * (N-1))
        ch.delta[i] <- rnorm(1, delta.mu, delta.sigma)

        # Phi_k step
        phi_k.shape <- (N - 1 + 2e-5) / 2
        phi_k.rate  <- (sumabs2(ch.kappa[2:N,i] - ch.delta[i] - ch.kappa[1:(N-1),i]) + 2e-5) / 2
        ch.phi_k[i] <- rgamma(1, phi_k.shape, phi_k.rate)

        # Phi_b step
        for (x in 1:X) {
            phi_b.shape   <- (N - 1 + 2e-5) / 2
            phi_b.rate    <- (sumabs2(ch.beta[x,2:N,i] - ch.beta[x,1:(N-1),i]) + 2e-5) / 2
            ch.phi_b[x,i] <- rgamma(1, phi_b.shape, phi_b.rate)
        }
    }

    #-- Wrap-up and return --#

    # Drop burnin
    drop.msk <- -(1:B)
    ch.phi   <- ch.phi[ ,drop.msk]
    ch.beta  <- ch.beta[ , ,drop.msk]
    ch.kappa <- ch.kappa[ ,drop.msk]
    ch.delta <- ch.delta[drop.msk]
    ch.phi_k <- ch.phi_k[drop.msk]
    ch.phi_b <- ch.phi_b[ ,drop.msk]

    # Return
    ret <- list(alpha = emp.alpha, beta = ch.beta, kappa = ch.kappa, phi = ch.phi,
                delta = ch.delta, phi_k = ch.phi_k, phi_b = ch.phi_b, data = Y)
    class(ret) <- "diblc"
    return(ret)
}

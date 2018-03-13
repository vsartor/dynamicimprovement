# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#' Dynamic Lee-Carter Model with Fixed Precisions
#'
#' Fully-dynamic extension of the Lee-Carter model but using fixed precisions
#' for observational _and_ evolutional errors.
#'
#' @param Y Numerical matrix with age information on rows and time information
#'          on columns. Must be log-mortality with valid and finite values.
#' @param phi Observational precisions.
#' @param phi_a Evolutional precisions for alpha.
#' @param phi_b Evolutional precisions for beta.
#' @param phi_k Evolutional precision for kappa.
#' @param phi_d Evolutional precision for delta.
#'
#' @return A list with online means and variances for each parameter and
#'         forecasts. Class `dyn_lc_fixed`.
#'
#' @importFrom assertthat assert_that
#'
#' @export
dyn_lc_fixed <- function(Y, phi = NULL, phi_a = NULL, phi_b = NULL,
                         phi_k = NULL, phi_d = NULL)
{
    tic <- getms()

    #-- Input validation: Part I --#

    assert_that(is.matrix(Y) & is.numeric(Y))
    N <- ncol(Y)
    A <- nrow(Y)
    assert_that(all(is.finite(Y)))
    if (any(Y >= 0)) {
        warning("Expected Y to be log-mortality. Non-negative values received.")
    }

    #-- Optional input entry --#

    if (is.null(phi))   phi   <- rep(400, A)
    if (is.null(phi_a)) phi_a <- rep(1 / 0.005^2, A)
    if (is.null(phi_b)) phi_b <- rep(1 / 0.005^2, A)
    if (is.null(phi_k)) phi_k <- 1 / 2
    if (is.null(phi_d)) phi_d <- 1 / 0.0001

    #-- Input validation: Part II --#

    assert_that(is.numeric(phi))
    assert_that(length(phi) == A)
    assert_that(all(phi > 0))

    assert_that(is.numeric(phi_a))
    assert_that(length(phi_a) == A)
    assert_that(all(phi_a > 0))

    assert_that(is.numeric(phi_b))
    assert_that(length(phi_b) == A)
    assert_that(all(phi_b > 0))

    assert_that(is.numeric(phi_k))
    assert_that(length(phi_k) == 1)
    assert_that(phi_k > 0)

    assert_that(is.numeric(phi_d))
    assert_that(length(phi_d) == 1)
    assert_that(phi_d > 0)

    #-- Constants --#

    # State-space dimension
    P <- 2*A + 2

    # Index mapping for parameters
    idx_a <- 1:A
    idx_b <- (A+1):(2*A)
    idx_k <- 2*A + 1
    idx_d <- 2*A + 2

    # Evolution matrix
    G               <- diag(P)
    G[idx_k, idx_d] <- 1
    G_t             <- trans(G)

    # Identity matrix of size A
    I_A <- diag(A)

    # Prior mean
    m_0        <- numeric(P)
    m_0[idx_a] <- apply(Y, 1, mean) # TODO: This is theoretically wrong.
    m_0[idx_b] <- 0.02
    m_0[idx_k] <- 20
    m_0[idx_d] <- -1

    # Prior variance
    C_0        <- numeric(P)
    C_0[idx_a] <- 0.001^2
    C_0[idx_b] <- 0.015^2
    C_0[idx_k] <- 20^2
    C_0[idx_d] <- 0.05^2
    C_0        <- diag(C_0)

    # Observational covariance matrix
    V <- diag(1 / phi)

    # Evolutional covariance matrix
    W <- diag(1 / c(phi_a, phi_b, phi_k, phi_d))

    #-- Allocate arrays --#

    # Filter means and variances
    a <- alloc(P, N)
    R <- alloc(P, P, N)
    m <- alloc(P, N)
    C <- alloc(P, P, N)

    # One-step ahead forecast summaries
    f_l <- f_m <- f_u <- alloc(A, N)

    #-- Filtering --#

    #-- First step

    # State priors
    a[ ,1]   <- G %*% m_0
    R[ , ,1] <- G %*% C_0 %*% G_t + W

    # Evolutional constants (Taylor approximation values)
    k_star <- a[idx_k,1]
    b_star <- a[idx_b,1]
    Ft     <- cbind(I_A, k_star*I_A, b_star, numeric(A))
    mu_nu  <- -b_star * k_star

    # One-step ahead forecasts
    f        <- Ft %*% a[ ,1] + mu_nu
    Q        <- Ft %*% R[ , ,1] %*% trans(Ft) + V
    f_m[ ,1] <- f
    f_l[ ,1] <- f - 1.96 * sqrt(diag(Q))
    f_u[ ,1] <- f + 1.96 * sqrt(diag(Q))

    # State posteriors
    e        <- Y[ ,1] - f
    B        <- R[ , ,1] %*% trans(Ft) %*% covinv(Q)
    m[ ,1]   <- a[ ,1] + B %*% e
    C[ , ,1] <- R[ , ,1] - B %*% Q %*% trans(B)

    # Standardization
    scale            <- sqrt(sumabs2(m[idx_b, 1]))
    m[idx_b,1]       <- m[idx_b,1] / scale
    C[idx_b,idx_b,1] <- C[idx_b,idx_b,1] / scale^2
    m[idx_k,1]       <- m[idx_k,1] * scale
    C[idx_k,idx_k,1] <- C[idx_k,idx_k,1] * scale^2

    #-- Other steps

    for (t in 2:N) {
        # State priors
        a[ ,t]   <- G %*% m[ ,t-1]
        R[ , ,t] <- G %*% C[ , ,t-1] %*% G_t + W

        # Evolutional constants (Taylor approximation values)
        k_star <- a[idx_k,t]
        b_star <- a[idx_b,t]
        Ft     <- cbind(I_A, k_star*I_A, b_star, numeric(A))
        mu_nu  <- -b_star * k_star

        # One-step ahead forecasts
        f        <- Ft %*% a[ ,t] + mu_nu
        Q        <- Ft %*% R[ , ,t] %*% trans(Ft) + V
        f_m[ ,t] <- f
        f_l[ ,t] <- f - 1.96 * sqrt(diag(Q))
        f_u[ ,t] <- f + 1.96 * sqrt(diag(Q))

        # State posteriors
        e        <- Y[ ,t] - f
        B        <- R[ , ,t] %*% trans(Ft) %*% covinv(Q)
        m[ ,t]   <- a[ ,t] + B %*% e
        C[ , ,t] <- R[ , ,t] - B %*% Q %*% trans(B)

        # Standardization
        scale            <- sqrt(sumabs2(m[idx_b, 1]))
        m[idx_b,t]       <- m[idx_b,t] / scale
        C[idx_b,idx_b,t] <- C[idx_b,idx_b,t] / scale^2
        m[idx_k,t]       <- m[idx_k,t] * scale
        C[idx_k,idx_k,t] <- C[idx_k,idx_k,t] * scale^2
    }

    #-- Return results --#

    r <- list(
        mean_a = m[idx_a, ],
        mean_b = m[idx_b, ],
        mean_k = m[idx_k, ],
        mean_d = m[idx_d, ],

        sd_a = sqrt(apply(C[idx_a,idx_a, ], 3, diag)),
        sd_b = sqrt(apply(C[idx_b,idx_b, ], 3, diag)),
        sd_k = sqrt(C[idx_k,idx_k, ]),
        sd_d = sqrt(C[idx_d,idx_d, ]),

        f_lower = f_l,
        f_mean  = f_m,
        f_upper = f_u,

        data = Y,

        elapsed = getms() - tic
    )
    class(r) <- c('dyn_lc_fixed', 'dyn_lc')
    r
}

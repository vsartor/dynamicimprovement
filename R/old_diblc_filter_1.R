# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#' Dynamic Improvement Bayesian Lee-Carter Model
#'
#' See `old_diblc_1` for further documentation.
#'
#' @param Y Numerical matrix with age information on rows and time information
#'          on columns. Must be log-mortality with valid and finite values.
#' @param phi Observational precisions.
#' @param df_a Discount factor for alpha errors. Defaults to 0.99.
#' @param df_b Discount factor for beta errors. Defaults to 0.90.
#' @param df_k Discount factor for kappa error. Defaults to 0.90.
#' @param df_d Discount factor for delta error. Defaults to 0.99.
#' @return A list with online means and variances for each parameter and
#'         forecasts. Same class as `old_diblc_filter_1`.
#'
#' @importFrom assertthat assert_that
#'
#' @export
old_diblc_filter_1 <- function(Y, phi = NULL, df_a = 0.99, df_b = 0.9,
                               df_k = 0.9, df_d = 0.99)
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

    if (is.null(phi)) phi <- rep(400, A)

    #-- Input validation: Part II --#

    assert_that(is.numeric(phi))
    assert_that(length(phi) == A)
    assert_that(all(phi > 0))

    assert_that(is.numeric(df_a))
    assert_that(length(df_a) == 1)
    assert_that(df_a > 0 & df_a < 1)

    assert_that(is.numeric(df_b))
    assert_that(length(df_b) == 1)
    assert_that(df_b > 0 & df_b < 1)

    assert_that(is.numeric(df_k))
    assert_that(length(df_k) == 1)
    assert_that(df_k > 0 & df_k < 1)

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
    G_t             <- t.default(G)

    # Identity matrix of size A
    I_A <- diag(A)

    # Prior mean
    m_0        <- numeric(P)
    m_0[idx_b] <- 0.02
    m_0[idx_k] <- 30
    m_0[idx_d] <- -0.2
    m_0[idx_a] <- apply(Y, 1, mean)

    # Prior variance
    C_0        <- numeric(P)
    C_0[idx_b] <- 0.015^2
    C_0[idx_k] <- 30^2
    C_0[idx_d] <- 0.1^2
    C_0[idx_a] <- 0.001^2
    C_0        <- diag(C_0)

    # Observational covariance matrix
    V <- diag(1 / phi)

    #-- Allocate arrays --#

    # Filter means and variances
    a <- alloc(P, N)
    R <- alloc(P, P, N)
    m <- alloc(P, N)
    C <- alloc(P, P, N)
    s <- alloc(P, N)
    S <- alloc(P, P, N)

    # One-step ahead forecast summaries
    f_l <- f_m <- f_u <- alloc(A, N)
    f_Q <- alloc(A, A, N)

    #-- Filtering --#

    #-- First step

    # State priors
    a[ ,1]   <- G %*% m_0
    Pt <- symmetrize(G %*% C_0 %*% G_t)
    Wt <- matrix(0, nrow = P, ncol = P)
    Wt[idx_a,idx_a] <- (1 - df_a) * Pt[idx_a,idx_a] / df_a
    Wt[idx_b,idx_b] <- (1 - df_b) * Pt[idx_b,idx_b] / df_b
    Wt[idx_k,idx_k] <- (1 - df_k) * Pt[idx_k,idx_k] / df_k
    Wt[idx_d,idx_d] <- (1 - df_d) * Pt[idx_d,idx_d] / df_d
    R[ , ,1] <- Pt + Wt

    # Evolutional constants (Taylor approximation values)
    k_star <- a[idx_k,1]
    b_star <- a[idx_b,1]
    Ft     <- cbind(I_A, k_star*I_A, b_star, numeric(A))
    mu_nu  <- -b_star * k_star

    # One-step ahead forecasts
    f          <- Ft %*% a[ ,1] + mu_nu
    Q          <- Ft %*% R[ , ,1] %*% t.default(Ft) + V
    f_m[ ,1]   <- f
    f_l[ ,1]   <- f - 1.96 * sqrt(diag(Q))
    f_u[ ,1]   <- f + 1.96 * sqrt(diag(Q))
    f_Q[ , ,1] <- Q

    # State posteriors
    e        <- Y[ ,1] - f
    B        <- R[ , ,1] %*% t.default(Ft) %*% covinv(Q)
    m[ ,1]   <- a[ ,1] + B %*% e
    C[ , ,1] <- R[ , ,1] - B %*% Q %*% t.default(B)

    # Standardization
    scale            <- sqrt(sumabs2(m[idx_b, 1])) * mode(sign(m[idx_b,1]))
    m[idx_b,1]       <- m[idx_b,1] / scale
    C[idx_b,idx_b,1] <- C[idx_b,idx_b,1] / scale^2
    m[idx_k,1]       <- m[idx_k,1] * scale
    C[idx_k,idx_k,1] <- C[idx_k,idx_k,1] * scale^2

    #-- Other steps

    for (t in 2:N) {
        # State priors
        a[ ,t]   <- G %*% m[ ,t-1]
        Pt <- symmetrize(G %*% C[ , ,t-1] %*% G_t)
        Wt <- matrix(0, nrow = P, ncol = P)
        Wt[idx_a,idx_a] <- (1 - df_a) * Pt[idx_a,idx_a] / df_a
        Wt[idx_b,idx_b] <- (1 - df_b) * Pt[idx_b,idx_b] / df_b
        Wt[idx_k,idx_k] <- (1 - df_k) * Pt[idx_k,idx_k] / df_k
        Wt[idx_d,idx_d] <- (1 - df_d) * Pt[idx_d,idx_d] / df_d
        R[ , ,t] <- Pt + Wt

        # Evolutional constants (Taylor approximation values)
        k_star <- a[idx_k,t]
        b_star <- a[idx_b,t]
        Ft     <- cbind(I_A, k_star*I_A, b_star, numeric(A))
        mu_nu  <- -b_star * k_star

        # One-step ahead forecasts
        f          <- Ft %*% a[ ,t] + mu_nu
        Q          <- Ft %*% R[ , ,t] %*% t.default(Ft) + V
        f_m[ ,t]   <- f
        f_l[ ,t]   <- f - 1.96 * sqrt(diag(Q))
        f_u[ ,t]   <- f + 1.96 * sqrt(diag(Q))
        f_Q[ , ,t] <- Q

        # State posteriors
        e        <- Y[ ,t] - f
        B        <- R[ , ,t] %*% t.default(Ft) %*% covinv(Q)
        m[ ,t]   <- a[ ,t] + B %*% e
        C[ , ,t] <- R[ , ,t] - B %*% Q %*% t.default(B)

        # Standardization
        scale            <- sqrt(sumabs2(a[idx_b, t])) * mode(sign(a[idx_b,t]))
        a[idx_b,t]       <- a[idx_b,t] / scale
        R[idx_b,idx_b,t] <- R[idx_b,idx_b,t] / scale^2
        a[idx_k,t]       <- a[idx_k,t] * scale
        R[idx_k,idx_k,t] <- R[idx_k,idx_k,t] * scale^2
        scale            <- sqrt(sumabs2(m[idx_b, t])) * mode(sign(m[idx_b,t]))
        m[idx_b,t]       <- m[idx_b,t] / scale
        C[idx_b,idx_b,t] <- C[idx_b,idx_b,t] / scale^2
        m[idx_k,t]       <- m[idx_k,t] * scale
        C[idx_k,idx_k,t] <- C[idx_k,idx_k,t] * scale^2
    }

    #-- Smoothing --#

    s[ ,N]   <- m[ ,N]
    S[ , ,N] <- C[ , ,N]

    for (i in (N-1):1) {
        B        <- C[ , ,i] %*% G_t %*% solve(R[ , ,i+1])
        s[ ,i]   <- m[ ,i] + B %*% (s[ ,i+1] - a[ ,i+1])
        S[ , ,i] <- C[ , ,i] - B %*% (R[ , ,i+1] - S[ , ,i+1]) %*% t.default(B)

        scale            <- sqrt(sumabs2(s[idx_b, i])) * mode(sign(s[idx_b, i]))
        s[idx_b,i]       <- s[idx_b,i] / scale
        S[idx_b,idx_b,i] <- C[idx_b,idx_b,i] / scale^2
        s[idx_k,i]       <- m[idx_k,i] * scale
        S[idx_k,idx_k,i] <- C[idx_k,idx_k,i] * scale^2
    }

    #-- Static standardization --#

    #scale            <- (s[idx_k,1] - s[idx_k,N]) / (N-1)
    #s[idx_b, ]       <- s[idx_b, ] * scale
    #S[idx_b,idx_b, ] <- S[idx_b,idx_b, ] * scale^2
    #s[idx_k, ]       <- s[idx_k, ] / scale
    #S[idx_k,idx_k, ] <- S[idx_k,idx_k, ] / scale^2

    #-- Return results --#

    m <- s
    C <- S

    r <- list(
        mean_a = m[idx_a, ],
        mean_b = m[idx_b, ],
        mean_k = m[idx_k, ],
        mean_d = m[idx_d, ],

        sd_a = sqrt(apply(C[idx_a,idx_a, ], 3, diag)),
        sd_b = sqrt(apply(C[idx_b,idx_b, ], 3, diag)),
        sd_k = sqrt(C[idx_k,idx_k, ]),
        sd_d = sqrt(C[idx_d,idx_d, ]),

        cov_a = C[idx_a,idx_a, ],
        cov_b = C[idx_b,idx_b, ],

        f_lower = f_l,
        f_mean  = f_m,
        f_upper = f_u,
        f_error = f_Q,

        data = Y,

        elapsed = getms() - tic
    )
    class(r) <- c('old_diblc_filter_1')
    r
}

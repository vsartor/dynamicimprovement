# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#' Displays graphical summaries for `dyn_lc_filter` results.
#'
#' @param x The `dyn_lc` result.
#' @param y The parameter of interest.
#' @param at If applicable, the time point of interest.
#' @param ages Vector with the ages if plotting 'beta' or 'alpha'. Optional.
#' @param years Vector with the years if plotting 'kappa'. Optional.
#' @param cred The credibility for the interval. Defaults to 0.95.
#' @param age The age _index_ to check one-step ahead forecasts.
#'
#' @return A `ggplot2` object.
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_ribbon aes
#' @importFrom ggplot2 theme_bw xlab ylab
#' @importFrom latex2exp TeX
#' @importFrom stats qnorm
#'
#' @export
plot.dyn_lc <- function(x, y, at = NULL, ages = NULL, years = NULL,
                        cred = 0.95, age = NULL, ...)
{
    sig <- 1 - cred
    if (is.null(ages)) ages <- 0:(nrow(x$data)-1)
    if (is.null(years)) years <- 1:ncol(x$data)

    if (y == "beta") {
        if (is.null(age)) {
            #-- Visualize beta at some time point --#

            m <- apply(x$beta[ ,at, ], 1, median)
            u <- apply(x$beta[ ,at, ], 1, quantile, 1 - sig/2)
            l <- apply(x$beta[ ,at, ], 1, quantile, sig/2)

            ggplot(NULL, aes(x = ages, y = m)) +
                theme_bw() +
                geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
                xlab('x') + ylab(TeX(sprintf("$\\beta_{%d,x}", at))) +
                geom_ribbon(aes(ymin = l, ymax = u), fill = 'darkblue', alpha = 0.2)
        } else {
            #-- Visualize beta evolution for some age --#

            m <- apply(x$beta[age, , ], 1, median)
            u <- apply(x$beta[age, , ], 1, quantile, 1 - sig/2)
            l <- apply(x$beta[age, , ], 1, quantile, sig/2)

            ggplot(NULL, aes(x = years, y = m)) +
                theme_bw() +
                geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
                xlab('t') + ylab(TeX(sprintf("$\\beta_{%d,t}", age))) +
                geom_ribbon(aes(ymin = l, ymax = u), fill = 'darkblue', alpha = 0.2)
        }
    } else if (y == "alpha") {
        #-- Visualize alpha at some timepoint --#

        m <- apply(x$alpha[ ,at, ], 1, median)
        u <- apply(x$alpha[ ,at, ], 1, quantile, 1 - sig/2)
        l <- apply(x$alpha[ ,at, ], 1, quantile, sig/2)

        ggplot(NULL, aes(x = ages, y = m)) +
            theme_bw() +
            geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
            xlab('x') + ylab(TeX(sprintf("$\\alpha_{%d,x}", at))) +
            geom_ribbon(aes(ymin = l, ymax = u), fill = 'darkblue', alpha = 0.25)

    } else if (y == "kappa") {
        #-- Visualize kappa --#

        m <- apply(x$kappa, 1, mean)
        u <- apply(x$kappa, 1, quantile, 1 - sig/2)
        l <- apply(x$kappa, 1, quantile, sig/2)

        ggplot(NULL, aes(x = years, y = m)) +
            theme_bw() +
            geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
            xlab('t') + ylab(TeX(sprintf("$\\kappa_{t}$"))) +
            geom_ribbon(aes(ymin = l, ymax = u), fill = 'darkblue', alpha = 0.25)

    } else {
        stop(sprintf("Invalid `y` value: %s", y))
    }
}

#' @export
print.dyn_lc <- function(x, ...)
{
    catf("Dynamic Improvement Lee Carter Extension")
    catf("----------------------------------------")
    catf("  Data with %d ages and %d years.", nrow(x$data), ncol(x$data))
    catf("  Took %.2fs to compute.", x$elapsed / 1000)
}

# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#' Displays graphical summaries for `dyn_lc_filter` results.
#'
#' @param x The `dyn_lc_filter` result.
#' @param y The parameter of interest.
#' @param at If applicable, the time point of interest.
#' @param ages Vector with the ages if plotting 'beta' or 'alpha'. Optional.
#' @param years Vector with the years if plotting 'kappa'. Optional.
#' @param cred The credibility for the interval. Defaults to 0.95.
#' @param age The age _index_ to check one-step ahead forecasts.
#' @param ... Unused.
#'
#' @return A `ggplot2` object.
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_ribbon aes
#' @importFrom ggplot2 theme_bw xlab ylab
#' @importFrom latex2exp TeX
#' @importFrom stats qnorm
#'
#' @export
plot.dyn_lc_filter <- function(x, y, at = NULL, ages = NULL, years = NULL,
                              cred = 0.95, age = NULL, ...)
{
    q0 <- qnorm(.5 + .5*cred)
    if (is.null(ages)) ages <- 0:(nrow(x$data)-1)
    if (is.null(years)) years <- 1:ncol(x$data)

    if (y == "beta") {
        if (is.null(age)) {
            #-- Visualize beta at some time point --#

            m <- x$mean_b[ ,at]
            s <- x$sd_b[ ,at]

            ggplot(NULL, aes(x = ages, y = m)) +
                theme_bw() +
                geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
                xlab('x') + ylab(TeX(sprintf("$\\beta_{%d,x}", at))) +
                geom_ribbon(aes(ymin = m - q0*s, ymax = m + q0*s),
                            fill = 'darkblue', alpha = 0.2)
        } else {
            #-- Visualize beta evolution for some age --#

            m <- x$mean_b[age, ]
            s <- x$sd_b[age, ]

            ggplot(NULL, aes(x = years, y = m)) +
                theme_bw() +
                geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
                xlab('t') + ylab(TeX(sprintf("$\\beta_{%d,t}", age))) +
                geom_ribbon(aes(ymin = m - q0*s, ymax = m + q0*s),
                            fill = 'darkblue', alpha = 0.2)
        }
    } else if (y == "alpha") {
        #-- Visualize alpha at some timepoint --#

        m <- x$mean_a[ ,at]
        s <- x$sd_a[ ,at]

        ggplot(NULL, aes(x = ages, y = m)) +
            theme_bw() +
            geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
            xlab('x') + ylab(TeX(sprintf("$\\alpha_{%d,x}", at))) +
            geom_ribbon(aes(ymin = m - q0*s, ymax = m + q0*s),
                        fill = 'darkblue', alpha = 0.25)

    } else if (y == "kappa") {
        #-- Visualize kappa --#

        m <- x$mean_k
        s <- x$sd_k

        ggplot(NULL, aes(x = years, y = m)) +
            theme_bw() +
            geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
            xlab('t') + ylab(TeX(sprintf("$\\kappa_{t}$"))) +
            geom_ribbon(aes(ymin = m - q0*s, ymax = m + q0*s),
                        fill = 'darkblue', alpha = 0.25)

    } else if (y == "forecast") {
        if (is.null(age)) {
            #-- Visualize one-step ahead forecast for certain timepoint --#

            tgt <- x$data[ ,at]

            m <- x$f_mean[ ,at]
            u <- x$f_upper[ ,at]
            l <- x$f_lower[ ,at]

            ggplot(NULL, aes(x = ages)) +
                theme_bw() +
                xlab('x') + ylab(TeX(sprintf("$Y_{x,%d}$", at))) +
                geom_point(aes(y = tgt)) +
                geom_line(aes(y = m), col = 'darkblue') +
                geom_ribbon(aes(ymin = l, ymax = u),
                            fill = 'darkblue', alpha = 0.25)
        } else {
            #-- Visualize one-step ahead forecast for a certain age --#

            tgt <- x$data[age, ]

            m <- x$f_mean[age, ]
            u <- x$f_upper[age, ]
            l <- x$f_lower[age, ]

            ggplot(NULL, aes(x = years)) +
                theme_bw() +
                xlab('t') + ylab(TeX(sprintf("$Y_{%d,t}$", age))) +
                geom_point(aes(y = tgt)) +
                geom_line(aes(y = m), col = 'darkblue') +
                geom_ribbon(aes(ymin = l, ymax = u),
                            fill = 'darkblue', alpha = 0.25)
        }
    } else {
        stop(sprintf("Invalid `y` value: %s", y))
    }
}

#' @export
print.dyn_lc_filter <- function(x, ...)
{
    catf("Dynamic Improvement Lee Carter Extension")
    catf("----------------------------------------")
    catf("  Data with %d ages and %d years.", nrow(x$data), ncol(x$data))
    catf("  Took %.2fms to compute.", x$elapsed)
    catf("  NOTE: Fixed precision variation.")
}

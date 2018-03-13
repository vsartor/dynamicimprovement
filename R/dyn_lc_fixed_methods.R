# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#' Displays graphical summaries for `dyn_lc_fixed` results.
#'
#' @param x The `dyn_lc_fixed` result.
#' @param y The parameter of interest.
#' @param at If applicable, the time point of interest.
#' @param ages Vector with the ages if plotting 'beta' or 'alpha'. Optional.
#' @param years Vector with the years if plotting 'kappa'. Optional.
#' @param cred The credibility for the interval. Defaults to 0.95.
#'
#' @return A `ggplot2` object.
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_ribbon aes
#' @importFrom ggplot2 theme_bw xlab ylab
#' @importFrom latex2exp TeX
#'
#' @export
plot.dyn_lc_fixed <- function(x, y, at = NULL, ages = NULL, years = NULL,
                              cred = 0.95, ...)
{
    q0 <- qnorm(.5 + .5*cred)

    if (y == "beta") {
        if (is.null(ages)) ages <- 0:(nrow(x$data)-1)

        m <- x$mean_b[ ,at]
        s <- x$sd_b[ ,at]

        ggplot(NULL, aes(x = ages, y = m)) +
            theme_bw() +
            geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
            xlab('x') + ylab(TeX(sprintf("$\\beta_{%d,x}", at))) +
            geom_ribbon(aes(ymin = m - q0*s, ymax = m + q0*s),
                        fill = 'darkblue', alpha = 0.25)
    } else if (y == "alpha") {
        if (is.null(ages)) ages <- 0:(nrow(x$data)-1)

        m <- x$mean_a[ ,at]
        s <- x$sd_a[ ,at]

        ggplot(NULL, aes(x = ages, y = m)) +
            theme_bw() +
            geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
            xlab('x') + ylab(TeX(sprintf("$\\alpha_{%d,x}", at))) +
            geom_ribbon(aes(ymin = m - q0*s, ymax = m + q0*s),
                        fill = 'darkblue', alpha = 0.25)
    } else if (y == "kappa") {
        if (is.null(years)) years <- 1:ncol(x$data)

        m <- x$mean_k
        s <- x$sd_k

        ggplot(NULL, aes(x = years, y = m)) +
            theme_bw() +
            geom_line(col = 'darkblue') + geom_point(col = 'darkblue') +
            xlab('t') + ylab(TeX(sprintf("$\\kappa_{t}$"))) +
            geom_ribbon(aes(ymin = m - q0*s, ymax = m + q0*s),
                        fill = 'darkblue', alpha = 0.25)

    } else {
        stop(sprintf("Invalid `y` value: %s", y))
    }
}

#' @export
print.dyn_lc_fixed <- function(x, ...)
{
    catf("Dynamic Improvement Lee Carter Extension")
    catf("----------------------------------------")
    catf("  Data with %d ages and %d years.", nrow(x$data), ncol(x$data))
    catf("  Took %.2fms to compute.", x$elapsed)
    catf("  NOTE: Fixed precision variation.")
}

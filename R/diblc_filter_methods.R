# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#' @export
print.diblc_filter <- function(x, ...)
{
    catf("Dynamic Improvement Lee Carter Extension Filter")
    catf("-----------------------------------------------")
    catf("  Data with %d ages and %d years.", nrow(x$data), ncol(x$data))
    catf("  Took %.2fms to compute.", x$elapsed)
}

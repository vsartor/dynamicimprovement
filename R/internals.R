# Copyright (C) 2018 Victhor Sartório
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

alloc <- function(...) array(dim = c(...))

covinv <- function(X) chol2inv(chol(X))

trans <- t.default

sumabs2 <- function(x) sum(x^2)

catf <- function(x, ...) cat(sprintf("%s\n", sprintf(x, ...)))

getms <- function() as.numeric(format(Sys.time(), "%OS5")) * 1000

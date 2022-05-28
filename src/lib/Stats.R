## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
if (!exists("Stats") || is.environment(Stats)) {
    Stats <- new.env(parent = emptyenv())
}

local({
    .VERSION = "0.1"

    AddMatRows <- function(X, Y) {
        if (!all.equal(colnames(X), colnames(Y))) {
            stop("X and Y must have the same column names!")
        }
        X <- as.matrix(X)
        Y <- as.matrix(Y)
        a <- rownames(X)
        b <- rownames(Y)
        c <- union(a, b)
        n <- length(c)
        d <- colnames(X)
        m <- length(d)
        Z <- matrix(0, ncol = m, nrow = n, dimnames = list(c, d))
        Z[a, ] <- Z[a, ] + X[a, ]
        Z[b, ] <- Z[b, ] + Y[b, ]
        Z
    }

    summarizeMatrix <- function(X, by, FUN = "mean", ncores = 1, na.rm = TRUE) {
        if (class(X) != "matrix") { X <- as.matrix(X) }
        n <- nrow(X)
        l <- length(by)
        if (n != l) { stop("nrow of matrix doesn't match length of by!") }
        idxBy <- split(1:nrow(X), f = by)
        if (ncores > 1) { 
            library("parallel")
            matsBy <- mclapply(idxBy, function(i) apply(X[i, , drop = FALSE], 2, FUN = FUN, na.rm = na.rm), mc.cores = ncores)
        } else {
            matsBy <- lapply(idxBy, function(i) apply(X[i, , drop = FALSE], 2, FUN = FUN, na.rm = na.rm))
        }
        mat <- do.call(rbind, matsBy)
        rownames(mat) <- names(idxBy)
        mat
    }

    giniCurve <- function(v) {
        n <- seq(v)
        s <- sum(v)
        v <- sort(v)
        p <- cumsum(v) / s
    }

    giniCoef <- function(v) {
        n <- length(v)
        Curve <- giniCurve(v)
        areaUnder <- sum(Curve)
        areaHalf <- 0.5 * n
        (areaHalf-areaUnder)/areaHalf
    }

    cv <- function(x, na.rm = TRUE) {
        sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
    }

    jaccard <- function(x, y, pseudocount = 0) {
        (length(intersect(x, y)) + pseudocount) / (length(union(x, y)) + pseudocount)
    }

    jaccard <- function(x, y, cutoff = 0, pseudocount = 0) {
        i <- which(x > cutoff)
        j <- which(y > cutoff)
        (length(intersect(i, j)) + pseudocount) / (length(union(i, j)) + pseudocount)
    }

    match_coef <- function(x, y, cutoff = 0) {
        b1 <- as.integer(x > cutoff)
        b2 <- as.integer(y > cutoff)
        both_present <- sum(b1 + b2 == 2)
        both_absent <- sum(b1 + b2 == 0)
        res <- both_present + both_absent
        res / length(x)
    }

    downsample <- function(x, n) {
        if (sum(x > 0) <= n) { return(x) }
        x[sample(which(x > 0), size = sum(x > 0)-n)] <- 0
        return(x)
    }

    hamming_dist <- function(x, y, cutoff = 0, normalize = FALSE) {
        bx <- as.integer(x > cutoff)
        by <- as.integer(y > cutoff)
        d <- sum(bx + by == 1)
        if (normalize) { d <- d / length(x) }
        d
    }

    downsample_hamming_dist <- function(x, y, dsth = 100, cutoff = 0, normalize = FALSE) {
        bx <- as.integer(x > cutoff)
        by <- as.integer(y > cutoff)
        sx0 <- sum(bx)
        sy0 <- sum(by)
        if (sx0 > dsth) { bx <- downsample(bx, n = dsth) }
        if (sy0 > dsth) { by <- downsample(by, n = dsth) }
        sx <- sum(bx)
        sy <- sum(by)
        loss_x <- 1 - sx/dsth
        loss_y <- 1 - sy/dsth
        d10 <- sum(bx == 1 & by == 0) / (1-loss_x) - sum(bx == 1 & by == 1) * loss_y / ((1-loss_x) * (1-loss_y))
        d01 <- sum(bx == 0 & by == 1) / (1-loss_y) - sum(bx == 1 & by == 1) * loss_x / ((1-loss_y) * (1-loss_x))
        d <- d10 + d01
        if (normalize) { d <- d / length(x) }
        d
    }

    downsample_simple_match <- function(x, y, cutoff = 0, nperms = 300, sumfn = median) {
        res <- replicate(nperms, {
            bx <- as.integer(x > cutoff)
            by <- as.integer(y > cutoff)
            if (sum(bx) != sum(by)) {
                if (sum(bx) < sum(by)) {
                    n <- sum(bx)
                    by <- downsample(by, n)
                } else {
                    n <- sum(by)
                    bx <- downsample(bx, n)
                }
            }
            both_present <- sum(bx + by == 2)
            both_absent <- sum(bx + by == 0)
            res <- both_present + both_absent
            res / length(x)
        })
        if (nperms > 1 && !is.null(sumfn)) { sumfn(res) } else { res }
    }

    downsample_jaccard <- function(x, y, cutoff = 0, nperms = 300, sumfn = median) {
        res <- replicate(nperms, {
            bx <- as.integer(x > cutoff)
            by <- as.integer(y > cutoff)
            if (sum(bx) != sum(by)) {
                if (sum(bx) < sum(by)) {
                    n <- sum(bx)
                    by <- downsample(by, n)
                } else {
                    n <- sum(by)
                    bx <- downsample(bx, n)
                }
            }
        i <- which(bx == 1)
        j <- which(by == 1)
        (length(intersect(i, j)) + pseudocount) / (length(union(i, j)) + pseudocount)
        })
        if (nperms > 1 && !is.null(sumfn)) { sumfn(res) } else { res }
    }

    for (obj in ls()) 
        assign(obj, get(obj), envir = Stats)
})

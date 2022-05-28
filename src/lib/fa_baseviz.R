#!/usr/bin/env perl
## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
## v0.2

if (!require("optparse")) { 
    install.packages("optparse")
    requireNamespace("optparse")
}

parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, "--faFile1", action = "store", type = "character", default = NULL, help = "R1 filename (FASTA or FASTQ)")
parser <- optparse::add_option(parser, "--faFile2", action = "store", type = "character", default = NULL, help = "R2 filename (FASTA or FASTQ)")
parser <- optparse::add_option(parser, "--sortBy", action = "store", type = "integer", default = 1, help = "in which mate's order (choose from 1 or 2)")
parser <- optparse::add_option(parser, "--primerFile", action = "store", type = "character", default = NULL, help = "primer filename (FASTA)")
parser <- optparse::add_option(parser, "--outFile", action = "store", type = "character", default = "output.png", help = "output filename (PNG)")
parser <- optparse::add_option(parser, "--main", action = "store", type = "character", default = "", help = "title")
parser <- optparse::add_option(parser, "--width", action = "store", type = "numeric", default = 7, help = "width (inch)")
parser <- optparse::add_option(parser, "--height", action = "store", type = "numeric", default = 7, help = "height (inch)")
parser <- optparse::add_option(parser, "--res", action = "store", type = "integer", default = 300, help = "resolution (dpi)")


opts <- optparse::parse_args(parser, args = commandArgs(trailingOnly = TRUE))

plotfasorted <- function(faFile1, faFile2 = NULL, colorLevels = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#444444"), sortBy, primerFile, outFile, main, width, height, res) {
    baseLevels <- c("A", "C", "G", "T", "N")
    k <- NA
    if (grepl("\\.fa$|\\.fasta$", faFile1, ignore.case = TRUE)) {
        k <- 2
    } else if (grepl("\\.fq$|\\.fastq$", faFile1, ignore.case = TRUE)) {
        k <- 4
    } else { 
        stop("Unknown file exension! Should be one of '.fa', '.fasta', 'fq', 'fastq'")
    }
    Seq1 <- readLines(faFile1)
    Seq1 <- Seq1[seq(2, length(Seq1), by = k)]
    readlen <- nchar(Seq1[1])
    if (is.null(faFile2)) {
        idx <- order(Seq1)
        Seq1 <- Seq1[idx]
        X <- t(sapply(Seq1, function(x) strsplit(x, split = "")[[1]]))
    } else {
        Seq2 <- readLines(faFile2)
        Seq2 <- Seq2[seq(2, length(Seq2), by = k)]
        idx <- if (sortBy == 1) { order(Seq1) } else { order(Seq2) }
        Seq1 <- Seq1[idx]
        Seq2 <- Seq2[idx]
        X1 <- t(sapply(Seq1, function(x) strsplit(x, split = "")[[1]]))
        X2 <- t(sapply(Seq2, function(x) strsplit(x, split = "")[[1]]))
        X <- if (sortBy == 1) { 
            cbind(X1, X2[, ncol(X2):1]) 
        } else { 
            cbind(X2, X1[, ncol(X1):1]) 
            }
    }
    rownames(X) <- NULL
    l <- ncol(X)
    Y <- factor(X, levels = baseLevels)
    Z <- matrix(as.integer(Y), ncol = l)

    png(file = outFile, width = ifelse(is.null(faFile2), width, width * 2) * res, height = height * res, res = res, bg = "transparent")
    image(t(Z[nrow(Z):1, ]), breaks = c(0, 1, 2, 3, 4, 5), col = colorLevels, axes = FALSE, main = main)
    axis(side = 1, 
         at = if (is.null(faFile2)) { 
             seq(0, readlen-1, by = 5)/(readlen-1) 
             } else { 
                 c(seq(0, readlen-1, by = 5)/(l-1), seq(readlen, l, by = 5)/(l-1)) 
            }, 
        labels = if (is.null(faFile2)) { 
            seq(1, readlen, by = 5) 
            } else { 
                c(seq(1, readlen, by = 5), seq(readlen, 1, by = -5)) 
            }, 
        las = 1)
    if (!is.null(faFile2)) {
        text(x = 0.25, y = 1.02, label = ifelse(sortBy == 1, "R1", "R2"), xpd = TRUE)
        text(x = 0.75, y = 1.02, label = ifelse(sortBy == 1, "R2", "R1"), xpd = TRUE)
    }
    if (!is.null(primerFile) && file.exists(primerFile)) {
        primers <- readLines(primerFile)
        primers <- primers[seq(2, length(primers), by = 2)]
        primers <- lapply(primers, function(p) {
            p <- sapply(strsplit(p, "")[[1]], toupper)
            factor(p, levels = baseLevels)
        })
        for (i in 1:length(primers)) {
            p <- primers[[i]]
            n <- length(p)
            for (j in 1:n) { 
                rect(xleft = (j*1.4-1.5)/(l-1), 
                    xright = (j*1.4)/(l-1), 
                    ytop = 0.03*i+1.01+2/200, 
                    ybottom = 0.03*i+1.01-2/200, 
                    xpd = TRUE, col = colorLevels[as.integer(p[j])], border = NA)
                text(x = (j*1.4-0.75)/(l-1), y = 0.03*i+1.01, xpd = TRUE, label = p[j], cex = 1) 
            }
        }
    } 
    legend(x = -0.05, y = 1.03, legend = baseLevels, fill = colorLevels, border= rep(NA, length(baseLevels)), box.col = NA, xpd = TRUE, ncol = 1)
    dev.off()
}

plotfasorted(faFile1 = opts$faFile1, faFile2 = opts$faFile2, sortBy = opts$sortBy, primerFile = opts$primerFile, outFile = opts$outFile, main = opts$main, width = opts$width, height = opts$height, res = opts$res)

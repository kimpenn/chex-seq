## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
#######################################################################
if (!exists("Genome") || is.environment(Genome)) { 
    Genome <- new.env(parent = emptyenv()) 
}

local({
    .VERSION = "4.3"
    library("GenomeInfoDb")
    library("GenomicRanges")
    library("rtracklayer")

    #######################################################################
    ## Create an internal reference to human, mouse and rat genome info
    ## based on Ensembl annotation
    #######################################################################
    ## Constant parameters of human and mouse genomes (build hg38, mm10)
    ## Get chromosome names and lengths from package GenomeInfoDb
    genomeBuilds <- c(human = "hg38", mouse = "mm10", rat = "rn6")
    GInfoDbHg38SeqInfo <- Seqinfo(genome = unname(genomeBuilds["human"]))
    GInfoDbHg38SeqLevels <- seqlevels(GInfoDbHg38SeqInfo)
    GInfoDbHg38SeqLengths <- seqlengths(GInfoDbHg38SeqInfo)
    GInfoDbMm10SeqInfo <- Seqinfo(genome = unname(genomeBuilds["mouse"]))
    GInfoDbMm10SeqLevels <- seqlevels(GInfoDbMm10SeqInfo)
    GInfoDbMm10SeqLengths <- seqlengths(GInfoDbMm10SeqInfo)
    GInfoDbRn6SeqInfo <- Seqinfo(genome = unname(genomeBuilds["rat"]))
    GInfoDbRn6SeqLevels <- seqlevels(GInfoDbRn6SeqInfo)
    GInfoDbRn6SeqLengths <- seqlengths(GInfoDbRn6SeqInfo)

    ## Main chromosomes only, excluding displaced contigs
    ## filter out chromosomes unplaced or not in the main scaffolds 
    MainHg38SeqLevels <- paste("chr", c(1:22, "X", "Y", "M"), sep = "")
    MainMm10SeqLevels <- paste("chr", c(1:19, "X", "Y", "M"), sep = "")
    MainRn6SeqLevels <- paste("chr", c(1:20, "X", "Y", "M"), sep = "")
    MainHg38SeqLengths <- GInfoDbHg38SeqLengths[MainHg38SeqLevels]
    MainMm10SeqLengths <- GInfoDbMm10SeqLengths[MainMm10SeqLevels]
    MainRn6SeqLengths <- GInfoDbRn6SeqLengths[MainRn6SeqLevels]
    MainHg38SeqInfo <- GInfoDbHg38SeqInfo[MainHg38SeqLevels]
    MainMm10SeqInfo <- GInfoDbMm10SeqInfo[MainMm10SeqLevels]
    MainRn6SeqInfo <- GInfoDbRn6SeqInfo[MainRn6SeqLevels]

    ## Combine human and mouse
    MainSeqInfo <- list(human = MainHg38SeqInfo, mouse = MainMm10SeqInfo, rat = MainRn6SeqInfo)
    MainSeqLevels <- list(human = MainHg38SeqLevels, mouse = MainMm10SeqLevels, rat = MainRn6SeqLevels)
    MainSeqLengths <- list(human = MainHg38SeqLengths, mouse = MainMm10SeqLengths, rat = MainRn6SeqLengths)

    ## Autosomes only
    AutoSeqLevels <- list(
        human = MainSeqLevels$human[MainSeqLevels$human != "chrX" & MainSeqLevels$human != "chrY"  & MainSeqLevels$human != "chrM"], 
        mouse = MainSeqLevels$mouse[MainSeqLevels$mouse != "chrX" & MainSeqLevels$mouse != "chrY"  & MainSeqLevels$mouse != "chrM"], 
        rat = MainSeqLevels$rat[MainSeqLevels$rat != "chrX" & MainSeqLevels$rat != "chrY"  & MainSeqLevels$rat != "chrM"]
    )
    AutoSeqLengths <- list(human = MainSeqLengths$human[AutoSeqLevels$human], 
                           mouse = MainSeqLengths$mouse[AutoSeqLevels$mouse], 
                           rat = MainSeqLengths$rat[AutoSeqLevels$rat]
    )

    ## Autosomes + chrX
    MainNoMYSeqLevels <- list(
        human = MainSeqLevels$human[MainSeqLevels$human != "chrY" & MainSeqLevels$human != "chrM"], 
        mouse = MainSeqLevels$mouse[MainSeqLevels$mouse != "chrY" & MainSeqLevels$mouse != "chrM"], 
        rat = MainSeqLevels$rat[MainSeqLevels$rat != "chrY" & MainSeqLevels$rat != "chrM"]
    )
    MainNoMYSeqLengths <- list(                       
        human = MainSeqLengths$human[MainNoMYSeqLevels$human], 
        mouse = MainSeqLengths$mouse[MainNoMYSeqLevels$mouse], 
        rat = MainSeqLengths$rat[MainNoMYSeqLevels$rat]
    )

    ## Remove unwanted chromosomes and rewrite seqinfo
    standardizeSeqInfo <- function(GRs, seqInfo, seqLevels = NULL, genomeBuild = NULL, prune = FALSE) {
        if (!is.null(seqLevels)) { 
            if (prune) {
                GRs <- GRs[seqnames(GRs) %in% seqLevels]
                GRs <- keepSeqlevels(GRs, pruning.mode = "coarse", value = intersect(seqLevels, seqlevels(GRs)))
            }
            seqlevels(GRs) <- seqLevels 
        }
        seqinfo(GRs) <- seqInfo
        if (!is.null(genomeBuild)) {
            genome(GRs) <- genomeBuild
        }
        GRs
    }
    #######################################################################
    ## IO for additional genomic data formats
    #######################################################################
    .extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
    .extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
    .extraCols_bedMethyl <- c(startCodon = "numeric", stopCodon = "numeric", itemRgb = "character", coverage = "numeric", methylPerc = "numeric")
    import.narrowPeak <- function(...) {
        library("rtracklayer")
        rtracklayer::import(..., format = "BED", extraCols = .extraCols_narrowPeak)
    }

    import.broadPeak <- function(...) {
        library("rtracklayer")
        rtracklayer::import(..., format = "BED", extraCols = .extraCols_broadPeak)
    }

    import.bedMethyl <- function(...) {
        library("rtracklayer")
        rtracklayer::import(..., format = "BED", extraCols = .extraCols_bedMethyl)
    }

    export.narrowPeak <- function(object, con) {
        df <- as.data.frame(object)
        df <- df[, c(1:3, 6, 7, 5, 8:11)]
        df[, 2] <- df[, 2] - 1 # Remember bed start is a half open value, say nucleotide #1 is denoted as (0, 2]
        df[, 2] <- format(df[, 2], scientific = FALSE)
        df[, 3] <- format(df[, 3], scientific = FALSE)
        write.table(df, file = con, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    }

    getGRsFromBed <- function(filename, SeqInfo = NULL, seqLevels = NULL, seqLengths = NULL) {
        library("rtracklayer")
        con <- file(filename, "r") # Because it seems rtracklayer does not close connections automatically, we need to manually do it.
        GRs <- rtracklayer::import(con, format = "bed")
        close(con)
        if (!is.null(seqLevels)) GRs <- GRs[seqnames(GRs) %in% seqLevels]
        GRs <- GRanges(seqnames = seqnames(GRs), ranges = ranges(GRs), strand = strand(GRs), mcols(GRs), seqinfo = SeqInfo, seqlengths = seqLengths)
    }

    write.igv <- function(sampleIDs, genomeBuild = "hg38", inDir, inFileBase = "uniqueNoDup", ext = "bed", pos = "chr1:158325836-158334836", snapshotDir = ".", filename = "script.igv", append = FALSE, verbose = FALSE, new = TRUE) {
        sampleNames <- paste("Sample", sampleIDs, sep = "_")
        filenames <-  paste(sampleNames, "star", inFileBase, ext, sep = ".")
        filenames <- file.path(inDir, filenames)
        if (verbose) {
            message("Processing ", filenames)
        }
        fnamecmd <- paste("load", filenames, collapse = "\n")
        if (new) {
            cmd <- sprintf("new\ngenome %s\nsnapShotDirectory %s\n%s\ngoto %s", genomeBuild, snapshotDir, fnamecmd, pos)
        } else {
            cmd <- sprintf("genome %s\nsnapShotDirectory %s\n%s\ngoto %s", genomeBuild, snapshotDir, fnamecmd, pos)
        } 
        cat(cmd, file = filename, append = append)
    }
    #######################################################################
    ## Genomic Range operations
    #######################################################################
    removeDupPeaks <- function(GRs) {
        ids <- mcols(GRs)[["name"]]
        i <- grepl("[0-9][a-z]+$", ids)
        if (sum(i) == 0) {
            return(GRs)
        } else {
            dupGRs <- GRs[i]
            dupIds <- ids[i]
            dupIds <- sub("[a-z]+$", "", dupIds)
            idxList <- split(seq(length(dupGRs)), f = dupIds)
            valList <- split(mcols(dupGRs)[["signalValue"]], f = dupIds)
            keptList <- sapply(names(valList), function(id) idxList[[id]][which.max(valList[[id]])])
            keptGRs <- dupGRs[keptList]
            GRs <- c(GRs[!i], keptGRs)
            return(sort(GRs, ignore.strand = TRUE))
        }
    }

    excludeGRs <- function(GRs1, GRs2, extend = 0, ignore.strand = FALSE) {
        GRs2 <- GRs2 + extend
        DGRs <- subsetByOverlaps(GRs1, GRs2, invert = TRUE, ignore.strand = ignore.strand)
    }

    getCvgsFromGRs <- function(GRs, seqLengths, seqLevels) {
        Cvgs <- coverage(GRs)
        seqLevelsCvg <- names(Cvgs)
        commSeqLevels <- intersect(seqLevelsCvg, seqLevels)
        if (length(commSeqLevels) < 1) { 
            warning("No coverage in any chromosome")
            Cvgs <- RleList()
        } else {
            seqLevelsDiff <- setdiff(seqLevels, seqLevelsCvg)
            if (length(seqLevelsDiff) > 0) {
                CvgsDiff <- sapply(seqLevelsDiff, function(seqLevel) { l <- seqLengths[seqLevel]; Rle(0, l) }, simplify = FALSE)
                Cvgs <- c(Cvgs, CvgsDiff)
            }
            Cvgs[seqLevels]
        }
        Cvgs
    }

    nTileRange <- function(Range, n) {
        x1 <- start(Range)
        x2 <- end(Range)
        w <- x2 - x1
        s <- w/n
        y1 <- x1 + s * (0:(n-1))
        y2 <- y1 + s
        IRanges(start = y1, end = y2)
    }

    nTileRanges <- function(Ranges, n) {
        stopifnot(n > 0)
        tiles <- lapply(Ranges, nTileRange, n = n)
        tiles <- as(tiles, "IRangesList")
    }

    getTilesFromRange <- function(from = 1, to, length.out, by) {
        if (missing(length.out)) {
            length.out <- floor((to - from) / by) + 1
        } else if (missing(by)) {
            by <- floor((to - from) / length.out) + 1
        }
        starts <- seq(from = from, to = to, by = by)
        ends <- starts + by - 1
        ends[length.out] <- to
        IRanges(starts, ends)
    }

    getTilesFromSeqLength <- function(seqLength, n, width, seqInfo = NULL) {
        library("GenomicRanges")
        seqLevel <- names(seqLength)
        IR <- getTilesFromRange(from = 1, to = seqLength[seqLevel], length.out = n, by = width)
        GRanges(seqname = seqLevel, IR, seqinfo = seqInfo)
    }

    getTilesFromSeqLengths <- function(seqLengths, ...) {
        seqLevels <- names(seqLengths)
        GRTiles <- lapply(seq_along(seqLengths), function(i) getTilesFromSeqLength(seqLengths[i], ...))
        names(GRTiles) <- seqLevels
        GRTiles <- as(GRTiles, "GRangesList")
    }

    rndGRangesFromGRanges <- function(universe, width = 5000, n = 3000, seed = 0) {
        w0 <- width(universe)
        idx <- w0 < width
        universe <- universe[!idx]
        m <- length(universe)
        chrs <- seqnames(universe)
        strands <- strand(universe)
        widths <- width(universe)
        starts <- start(universe)
        ends <- end(universe)
        if (m == 0) { stop("Nothing left after filtering regions shorter than the prespecified width") } 
        k <- 0
        set.seed(seed)
        GRsList <- vector(n, mode = "list")
        ## range i being selected: pr_i = width_i / sum(width_i)
        ## a point j in range i being selected: pp_j = 1/width_i
        ## point (i, j) prob. of being selected: pr_i * pp_j = 1/sum(width_i) = const.
        while (k < n) {
            i <- sample(1:m, size = 1, prob = widths)
            s <- sample(starts[i]:ends[i], size = 1)
            e <- s + width - 1
            if (e > ends[i]) { next }
            k <- k + 1
            GRsList[[k]] <- GRanges(chrs[i], IRanges(s, e), strands[i])
            if (k %% 1000 == 0) { message(k, " regions have been sampled.") }
        }
        GRsList <- as(GRsList, "GRangesList")
        unlist(GRsList)
    }

    createRandGRs <- function(seqLevels, seqLengths, n = 1000, widthMean = 200, stranded = FALSE) {
        snames <- sample(seqLevels, n, replace = TRUE)
        nBySeqLevel <- tapply(seq_along(snames), snames, length)
        starts <- unlist(lapply(names(nBySeqLevel), function(sname) { limit <- seqLengths[sname]; sample(1:limit, nBySeqLevel[sname], replace = TRUE) }))
        widths <- rpois(lambda = widthMean, n = n)
        ends <- starts + widths
        snames <- rep(names(nBySeqLevel), unlist(nBySeqLevel))
        ends <- ifelse(ends > seqLengths[snames], seqLengths[snames], ends)
        if (stranded) {
            strands <- sample(c("+", "-"), n, replace = TRUE)
        } else {
            strands <- rep("*", n)
        }
        GRs <- GRanges(seqnames = snames, IRanges(start = starts, end = ends), strand = strands)
        GRs <- GenomicRanges::sort(GRs)
    }

    flipStrands <- function(GRs, undefined = c("-", "*", "+")) {
        undefined <- match.arg(undefined)
        strands <- strand(GRs)
        strandsFlipped <- ifelse(strands == "-", "+", ifelse(strands == "+", "-", undefined))
        strand(GRs) <- strandsFlipped
        GRs$OldStrand <- strands
        GRs
    }

    getUpstream <- function(GRs, upstream = 2500) {
        library("GenomicRanges")
        GRsUpstream <- GenomicRanges::promoters(GRs, upstream = upstream, downstream = 0)
        trim(GRsUpstream)
    }

    getDownstream <- function(GRs, downstream = 2500) {
        GRsFlipped <- flipStrands(GRs, undefined = "-")
        GRsDownstream <- getUpstream(GRsFlipped, upstream = downstream)
        strand(GRsDownstream) <- GRsDownstream$OldStrand
        GRsDownstream$OldStrand <- NULL
        trim(GRsDownstream)
    }

    getTSSFromGRs <- function(GRs, width = 250, both = TRUE) {
        flank(GRs, width = width, both = both) 
    }

    ## Because a gene has multiple TSS sites will have the left-most one listed in EnsGene, we hence have to 
    ## use EnsTranscript to derive TSS flanking regions for each gene
    getTSSFlanksFromTxGRs <- function(TxGRs, width = 250, both = TRUE, by = c("gene", "transcript")) {
        by <- match.arg(by)
        flanks <- flank(TxGRs, width = width, both = both)
        if (by == "gene") {
            flanks <- trim(flanks)
            flanks <- split(flanks, f = mcols(flanks)[["gene_id"]])
            flanks <- unique(flanks)
        } else {
            flanks <- split(flanks, f = mcols(flanks)[["tx_id"]])
        }
        flanks
    }

    ## Refer to https://stackoverflow.com/questions/29253412/finding-intergenic-regions
    ## Note, TxDb.Hsapiens.UCSC.hg38.knownGene has 7 genes longer than 1e7bp, which are not real.
    getGenomicFeaturesIntergenic <- function(GRsGene, geneLengthLimits = c(0, Inf)) {
        l <- width(GRsGene)
        idx <- l > geneLengthLimits[1] & l < geneLengthLimits[2]
        GRsGeneFiltered <- GRsGene[idx]
        genic <- reduce(GRsGeneFiltered, ignore.strand = TRUE)
        intergenic <- gaps(genic)
        intergenic[as.character(strand(intergenic)) == "*"]
    }

    extendFeatures <- function(features, locType = c("start", "end", "middle"), extType = c("point", "range"), upstream = 5000, downstream = 5000) {
        if (extType == "point") {
            featureLoc <- switch(locType,
                start = featureStarts, 
                end = featureEnds, 
                middle = featureMids
            )
            featuresUpstream <- trim(GRanges(seqnames = featureSeqnames, IRanges(start = featureLoc - upstream, end = featureLoc), strand = featureStrands, ID = featureIDs, seqinfo = featureSeqinfo))
            featuresDownstream <- trim(GRanges(seqnames = featureSeqnames, IRanges(start = featureLoc, end = featureLoc + downstream), strand = featureStrands, ID = featureIDs, seqinfo = featureSeqinfo))
            featuresExt <- list(upstream = featuresUpstream, downstream = featuresDownstream)
        } else if (extType == "range") {
            featuresUpstream <- trim(GRanges(
                seqnames = featureSeqnames, 
                IRanges(start = ifelse(featureStrands == "-", featureStarts, featureStarts - upstream),
                        end = ifelse(featureStrands == "-", featureStarts + upstream, featureStarts)), 
                strand = featureStrands, ID = featureIDs, seqinfo = featureSeqinfo
            ))
            featuresDownstream <- trim(GRanges(
                seqnames = featureSeqnames, 
                IRanges(start = ifelse(featureStrands == "-", featureEnds - downstream, featureEnds),
                        end = ifelse(featureStrands == "-", featureEnds, featureEnds + downstream)), 
                strand = featureStrands, ID = featureIDs, seqinfo = featureSeqinfo
            ))
            featuresExt <- list(upstream = featuresUpstream, geneBody = features, downstream = featuresDownstream)
        }
        featuresExt
    }

    ## Suppose there are two isoform transcripts for a given gene, '+' means exons, '-' introns
    ## Suppose we apply this function to introns
    ##      |->+++++++--------++++-------++++++++| (E1, I1, E2, I2, E3) 
    ##      |->+++++++--------+++++++------------| (E4, I3, E5, I4)
    ## none: I1, I2, I3, I4
    ## union: I1, union(I2, I4)
    ## unique: I1, I2, I4
    ## intersection: I1, intersect(I2, I4)
    ## uniqueSubID: given unique, if I1 I3 have same genomic range but different subID (e.g. exon ID), they are treated as different exons and I3 hence won't be removed. 
    ## summarizeGRsBy(GRsBy = EnsHg38ExonsByTranscript, summary = "uniqueSubID", IDs = mcols(EnsHg38Transcripts[names(EnsHg38ExonsByTranscript)])[["gene_id"]], IDcolname = NULL, subIDcolname = "exon_id")
    summarizeGRsBy <- function(GRsBy, IDs, summary = c("unique", "uniqueSubID", "intersection", "union", "none"), IDcolname = NULL, subIDcolname = NULL, ncores = 1) {
        if (ncores > 1) { library("parallel") }
        summary <- match.arg(summary)
        n1 <- length(GRsBy) 
        if (missing(IDs) && !is.null(IDcolname)) {
            if (IDcolname %in% names(mcols(GRsBy))) {
                IDs <- mcols(GRsBy)[[IDcolname]]
            } else {
                stop("IDs are missing and IDcolname is not found in the metadata!")
            }
        } else {
            n2 <- length(IDs)
            if (n1 != n2) { stop("GRsBy has different length than IDs") }
        }
        groups <- split(seq(n1), f = IDs)
        if (summary == "union") {
            message("Started mapping, summary = ", summary, "...")
            if (ncores > 1) {
                GRsByID <- mclapply(groups, function(idx) {  
                    GRs <- GRsBy[idx]
                    reduce(unlist(GRs))
                }, mc.cores = ncores)
            } else {
                GRsByID <- lapply(groups, function(idx) {  
                    GRs <- GRsBy[idx]
                    reduce(unlist(GRs))
                })
            }
        } else if (summary == "intersection") {
            message("Started mapping, summary = ", summary, "...")
            if (ncores > 1) {
                GRsByID <- mclapply(groups, function(idx) {  
                    GRs <- GRsBy[idx]
                    Reduce(intersect, GRs)
                }, mc.cores = ncores)
            } else {
                GRsByID <- lapply(groups, function(idx) {  
                    GRs <- GRsBy[idx]
                    Reduce(intersect, GRs)
                })
            }
        } else if (summary == "unique") {
            message("Started mapping, summary = ", summary, "...")
            if (ncores > 1) {
                GRsByID <- mclapply(groups, function(idx) {  
                    GRs <- GRsBy[idx]
                    unname(unique(unlist(GRs)))
                }, mc.cores = ncores)
            } else {
                GRsByID <- lapply(groups, function(idx) {  
                    GRs <- GRsBy[idx]
                    unname(unique(unlist(GRs)))
                })
            }
        } else if (summary == "uniqueSubID") {
            if (is.null(subIDcolname)) { stop("subIDcolname cannot be NULL when summary is uniqueSubID") }
            message("Started mapping, summary = ", summary, "...")
            if (ncores > 1) {
                GRsByID <- mclapply(groups, function(idx) {
                    GRs <- GRsBy[idx]
                    GRs <- unlist(GRs)
                    subIDs <- mcols(GRs)[[subIDcolname]]
                    idx <- duplicated(subIDs)
                    GRs[!idx]
                }, mc.cores = ncores)
            } else {
                GRsByID <- lapply(groups, function(idx) {
                    GRs <- GRsBy[idx]
                    GRs <- unlist(GRs)
                    subIDs <- mcols(GRs)[[subIDcolname]]
                    idx <- duplicated(subIDs)
                    GRs[!idx]
                })
            }
        } else if (summary == "none") {
            message("Started mapping, summary = ", summary, "...")
            if (ncores > 1) {
                GRsByID <- mclapply(groups, function(idx) {  
                    GRs <- GRsBy[idx]
                    unname(unlist(GRs))
                }, mc.cores = ncores)
            } else {
                GRsByID <- lapply(groups, function(idx) {  
                    GRs <- GRsBy[idx]
                    unname(unlist(GRs))
                })
            }
        }
        message("Converting it to GRangesList...")
        GRsByID <- as(GRsByID, "GRangesList")
    }

    ## Get introns for a list of transcripts
    ## ExonsByTranscript: a GRangesList object or a list of GRanges
    ## Note: mclapply does not gain us advantage because it has to convert GRangesList to list first, which is a very time-consuming step; Do not set ncores > 1!
    getIntronsByTranscript <- function(exonsByTranscript, transcripts = NULL, method = c("psetdiff", "gaps"), ncores = 1) {
        method <- match.arg(method)
        if (method == "psetdiff") {
            if (is.null(transcripts)) { stop("When psetdiff is used, transcripts should be provided!") }
            strands_by_tx <- structure(as.character(strand(transcripts)), names = mcols(transcripts)[["tx_id"]])
            transcripts <- transcripts[match(names(exonsByTranscript), mcols(transcripts)[["tx_id"]])]
            intronsByTranscript <- psetdiff(transcripts, exonsByTranscript)
            intronsByTranscript <- intronsByTranscript[lengths(intronsByTranscript) > 0]
            L <- lengths(intronsByTranscript)
            N <- names(intronsByTranscript)
            idx <- mapply(FUN = function(s, l) { ifelse(rep(s, l) == "-", l:1, 1:l)}, s = strands_by_tx[N], l = L, SIMPLIFY = FALSE)
            idx1 <- mapply(FUN = function(x, y) { x + y }, x = idx, y = c(0, cumsum(head(L, -1))), SIMPLIFY = FALSE)
            intronsByTranscript <- unlist(intronsByTranscript)[unlist(idx1)]
            mcols(intronsByTranscript)[["intron_rank"]] <- as.integer(unlist(lapply(unname(L), seq)))
            N <- names(intronsByTranscript)
            names(intronsByTranscript) <- NULL
            intronsByTranscript <- split(intronsByTranscript, f = N)
        } else if (method == "gaps") {
            if (ncores > 1) { 
                library("parallel")
                intronsByTranscript <- mclapply(exonsByTranscript, FUN = getIntrons, mc.cores = ncores)
            } else {
                intronsByTranscript <- lapply(exonsByTranscript, FUN = getIntrons)
            }
        }
        as(intronsByTranscript, "CompressedGRangesList")
    }
    
    ## union: the union of introns per transcript/isoform of a gene, requires intronsByTranscript as input
    ## intersection: the intersection of introns per transcript/isoform of a gene
    getIntronsByGene <- function(intronsByTranscript, IDs, summary = c("unique", "uniqueSubID", "intersection", "union", "none"), IDcolname = NULL, subIDcolname = NULL, standardizeSeqInfo = FALSE, seqInfo = NULL, seqLevels = NULL, ncores = 1) {
        summary <- match.arg(summary)
        intronsByGene <- summarizeGRsBy(intronsByTranscript, IDs = IDs, IDcolname = IDcolname, subIDcolname = subIDcolname, summary = summary, ncores = ncores)
        if (standardizeSeqInfo) {
            intronsByGene <- standardizeSeqInfo(intronsByGene, seqInfo = seqInfo, seqLevels = seqLevels)
        }
        intronsByGene
    }

    DFrameList2GRanges <- function(DFrame, seqLevel, sampleIDs) {
        runLengths <- DFrame[, "runLength"]
        start <- cumsum(c(0, runLengths[-length(runLengths)])) + 1
        Data <- DFrame[, -which(colnames(DFrame) == "runLength")]
        colnames(Data) <- sampleIDs
        GRs <- GRanges(seqLevel, IRanges(start, width = runLengths))
        mcols(GRs) <- Data
        GRs
    }

    #######################################################################
    ## Genomic Range statistics
    #######################################################################
    testGRsOverlap <- function(query, subject, universe, ignore.strand = TRUE) {
        x1 <- sum(width(intersect(query, subject, ignore.strand = ignore.strand)))
        x2 <- sum(width(intersect(query, universe, ignore.strand = ignore.strand)))
        y1 <- sum(width(intersect(subject, universe, ignore.strand = ignore.strand)))
        y2 <- sum(width(universe))
        log2OR <- log2(1+x1) - log2(1+x2) - log2(1+y1) + log2(1+y2)
        m <- matrix(c(x1, x2, y1, y2), ncol = 2)
        res <- chisq.test(m)
        pval <- res$p.value
        c(log2OR = log2OR, pval = pval)
    }

    #######################################################################
    ## Summarize read counts per genomic feature
    #######################################################################
    getPrimingRatesPerFeature <- function(featureRangesList, featureTypesByGene, GRsList, sampleIDs = NULL, ignore.strand = TRUE) { ## make sure we count a read as many times as the number of genes overlapping it per Junhyong's suggestion June 2019 
        if (is.null(sampleIDs)) { sampleIDs <- names(GRsList) }
        res <- sapply(featureTypesByGene, function(featureType) {
            features <- featureRangesList[[featureType]]
            featureIDs <- names(features)
            nfeatures <- length(features)
            if (inherits(features, "GRangesList")) { features <- reduce(features) }
            sapply(sampleIDs, function(sampleID) {
                message(featureType, " ", sampleID)
                GRs <- GRsList[[sampleID]]
                hits <- findOverlaps(subject = features, query = GRs, ignore.strand = ignore.strand)
                sh <- subjectHits(hits)
                qh <- queryHits(hits)
                cnts <- tapply(qh, INDEX = sh, FUN = length)
                zeros <- setdiff(as.character(1:nfeatures), names(cnts))
                l <- length(zeros)
                if (l > 0) {
                    x <- structure(rep(0, l), names = zeros)
                    cnts <- c(cnts, x)
                }
                cnts <- cnts[as.character(1:nfeatures)]
                names(cnts) <- featureIDs
                cnts
            })
        }, simplify = FALSE)
    }

    getCountsByStrandedness <- function(GRs, features, strandedness = c("sense", "antisense")) {
        strandedness <- match.arg(strandedness)
        nfeatures <- length(features)
        featureIDs <- names(features)
        if (strandedness == "antisense") { GRs <- flipStrands(GRs) }
        hits <- GenomicRanges::findOverlaps(subject = features, query = GRs, ignore.strand = FALSE)
        sh <- subjectHits(hits)
        qh <- queryHits(hits)
        cnts <- tapply(qh, INDEX = sh, FUN = length)
        zeros <- setdiff(as.character(1:nfeatures), names(cnts))
        l <- length(zeros)
        if (l > 0) {
            x <- structure(rep(0, l), names = zeros)
            cnts <- c(cnts, x)
        }
        cnts <- cnts[as.character(1:nfeatures)]
        names(cnts) <- featureIDs
        cnts
    }

    getCntsListFromBam <- function(bamfile, binsize = 5000, endType = "PE", seqLevels = NULL, ignore.strand = TRUE, bamFlags = list(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,  hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA, isFirstMateRead = NA, isSecondMateRead = NA, isSecondaryAlignment = NA, isNotPassingQualityControls = NA, isDuplicate = NA)) {
        BamFlag <- do.call(scanBamFlag, bamFlags)
        SeqInfo <- seqinfo(BamFile(bamfile))
        if (is.null(seqLevels)) { seqLevels <- seqlevels(SeqInfo) } 
        else { SeqInfo <- SeqInfo[seqLevels, ] }
        seqLengths <- seqlengths(SeqInfo)
        which <- GRanges(seqLevels, IRanges(1, seqLengths))
        BamParam <- ScanBamParam(flag = BamFlag, which = which)
        if (endType == "PE") {
            GAs <- readGAlignmentPairs(bamfile, param = BamParam)
        } else {
            GAs <- readGAlignments(bamfile, param = BamParam)
        }
        Tiles <- getTilesFromSeqLengths(seqLengths, width = binsize, seqInfo = SeqInfo)
        sapply(seqLevels, function(seqLevel) assay(summarizeOverlaps(reads = GAs, features = Tiles[[seqLevel]], ignore.strand = ignore.strand))[, 1], simplify = FALSE)
    }
    
    for (obj in c( ls() )) assign(obj, get(obj), envir = Genome)
})

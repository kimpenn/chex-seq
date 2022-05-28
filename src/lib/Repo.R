## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
if (!exists("Repo") || is.environment(Repo)) { 
    Repo <- new.env(parent = emptyenv())
}

local({
    .LOGFILE <- stderr()
    .VERSION = "0.3"
    
    catLog <- function(...) { 
        cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S %z]"), ...) 
        }
    
    getVERSEcounts <- function(exptName, sampleID, baseDir = "/lab/repo", analyzedDir = "analyzed", countsDir = "verse", fileBase = "verse", 
        features = c("exon", "exon-lev3", "anti-exon", "intron", "intron-lev1-lev2", 
                     "intron-lev3", "anti-intron", "intergenic", "lines_sines", "mito", "anti-mito")) {
        file.sample <- paste("Sample", ifelse(grepl(" ", sampleID), toupper(gsub(" ", "", sampleID)), sampleID), sep = "_")
        files.counts <- paste(file.sample, fileBase, features, "cnts", "txt", sep = ".")
        filenames <- file.path(baseDir, exptName, analyzedDir, file.sample, countsDir, files.counts)
        counts <- lapply(seq_along(features), function(i) {
            filename <- filenames[i]
            feature <- features[i]
            if (file.exists(filename)) {
                catLog("Reading", filename, "...\n", file = .LOGFILE)
                counts <- read.csv(filename, head = T, row.names = 1, sep = "\t")
            } else {
                catLog("Not found", filename, "\n", file = .LOGFILE)
                counts <- data.frame()
            }
        })
        counts <- `names<-`(counts, features)
    }
    
    getHTSeqCounts <- function(EID, SID, features = c(
            "exon", "exon-lev3", "anti-exon", 
            "intron", "intron-lev1-lev2", "intron-lev3", "anti-intron", 
            "intergenic", "lines_sines", "mito", "anti-mito"),
            base = "/lab/repo") {
        feature2filenameMap <- c(
            "exon" = "exons", "exon-lev3" = "exons-lev3", "anti-exon" = "anti-exons",
            "intron" = "introns", "intron-lev1-lev2" = "introns-lev1-lev2", 
            "intron-lev3" = "introns-lev3", "anti-intron" = "anti-introns", 
            "intergenic" = "intergenic", "lines_sines" = "lines_sines", 
            "mito" = "mito", "anti-mito" = "anti-mito")
        file.exp <- paste("E", EID, sep = ".")
        ## Fix a couple of samples whose "Sample Name" is different from the directory name, by concatenating them in upper case without any space
        file.sample <- paste("Sample", ifelse(grepl(" ", SID), toupper(gsub(" ", "", SID)), SID), sep = "_") 
        files.counts <- paste(file.sample, "htseq", feature2filenameMap[features], "cnts", "txt", sep = ".")
        filenames <- file.path(
            base, file.exp, "analyzed", file.sample, "htseq", files.counts)
        counts <- lapply(seq_along(features), function(i) {
            filename <- filenames[i]
            feature <- features[i]
            if (file.exists(filename)) {
                catLog("Reading", filename, "...\n", file = .LOGFILE)
                counts <- read.csv(filename, head = T, row.names = 1, sep = "\t")
            } else {
                catLog("Not found", filename, "\n", file = .LOGFILE)
                counts <- data.frame()
            }
        })
        counts <- `names<-`(counts, features)
    }

    getNGSCounts <- function(EID, SID, features = c(
            "exon", "exon-lev3", "anti-exon", 
            "intron", "intron-lev1-lev2", "intron-lev3", "anti-intron", 
            "intergenic", "lines_sines", "mito", "anti-mito"), 
            base = "/lab/repo") {
        feature2filenameMap <- c(
            "exon" = "exons", "exon-lev3" = "exons-lev3", "anti-exon" = "anti-exons",
            "intron" = "introns", "intron-lev1-lev2" = "introns-lev1-lev2", 
            "intron-lev3" = "introns-lev3", "anti-intron" = "anti-introns", 
            "intergenic" = "intergenic", "lines_sines" = "lines_sines", 
            "mito" = "mito", "anti-mito" = "anti-mito")
        file.exp <- paste("E", EID, sep = ".")
        file.sample <- paste("Sample", ifelse(grepl(" ", SID), toupper(gsub(" ", "", SID)), SID), sep = "_")
        counts <- lapply(seq(features), function(i) {
            feature <- features[i]
            file.counts_verse <- paste(file.sample, "verse", feature, "cnts", "txt", sep = ".")
            filename_verse <- file.path(base, file.exp, "analyzed", file.sample, "verse", file.counts_verse)
            file.counts_htseq <- paste(file.sample, "htseq", feature, "cnts", "txt", sep = ".")
            filename_htseq <- file.path(base, file.exp, "analyzed", file.sample, "htseq", file.counts_htseq)
            if (file.exists(filename_verse)) {
                catLog("Reading", filename_verse, "...\n", file = .LOGFILE)
                counts <- read.csv(filename_verse, head = T, row.names = 1, sep = "\t")
            } else if (file.exits(filename_htseq)) {
                catLog("Reading", filename_htseq, "...\n", file = .LOGFILE)
                counts <- read.csv(filename_htseq, head = T, row.names = 1, sep = "\t")
            } else {
                counts <- data.frame()
                catLog("Not found:", filename_verse, "or", filename_htseq, "...\n", file = .LOGFILE)
            }
        })
        counts <- `names<-`(counts, features)
    }

    fixStatsHeader <- function(stats) {
        colnames_lev1 <- colnames(stats)
        colnames_lev2 <- stats[1, ]
        colnames_lev1[is.na(colnames_lev1)] <- ""
        colnames_lev1 <- sub("^\\s+|\\s+$", "", colnames_lev1, perl = T)
        colnames_lev2[is.na(colnames_lev2)] <- ""
        colnames_lev2 <- sub("^\\s+|\\s+$", "", colnames_lev2, perl = T)
        colnames_fixed <- paste(unlist(rep(colnames_lev1[colnames_lev1 != ""], diff(c(which(colnames_lev1 != ""), length(colnames_lev1) + 1)))), colnames_lev2, sep = "_")
        colnames_fixed <- sub("_$", "", colnames_fixed)
        colnames(stats) <- colnames_fixed
        stats <- stats[-1, ]
    }

    ## Compared with Jamie's stats output
    ## We have 4 additional fields
    ## * Deletion-per-base Perc	
    ## * Avg Deletion Length	
    ## * Insertion-per-base Perc	
    ## * Avg Insertion Length
    getSTARstats <- function(exptName, sampleID, baseDir = "/lab/repo", analyzedDir = "analyzed", alignedDir = "star", fileBase = alignedDir, .DEBUG = FALSE) {
        sampleName <- paste("Sample", gsub(" ", "", sampleID), sep = "_")
        fname_version <- paste(sampleName, "versions", sep = ".")
        fname_stats <- paste(sampleName, fileBase, "stats", "txt", sep = ".")
        filename_version <- file.path(baseDir, exptName, analyzedDir, sampleName, alignedDir, fname_version)
        filename_stats <- file.path(baseDir, exptName, analyzedDir, sampleName, alignedDir, fname_stats)
        if (file.exists(filename_version)) {
            sect1 <- as.matrix(read.table(filename_version, sep = "\t", as.is = T, check = F))[2, -c(2, 4)]
            ver <- sub("-\\w+$", "", sect1[1])
            ## Pipeline version 2.0 has the issue that "SE?" and "readLength" are switched
            ##  pipeline                   version                   version
            ##       "2"            "STAR_2.4.0h1"                     "1.1"
            ##   species                       SE?                readLength
            ## "hg38.gencode21.stranded"     "100"                   "false"
            ## Therefore should be switched back
            if (ver == "2") {
                sect1[5:6] <- sect1[6:5]
                if (sect1[5] == "false") sect1[5] <- 0
                if (sect1[5] == "true") sect1[5] <- 1
            }
            if (ver < "2.1" && (length(sect1) < 6)) {
                # No columns "SE", "readLength"
                sect1 <- c(sect1, rep("", 6 - length(sect1))) 
            }
        } else {
            sect1 <- rep("", 6)
        }
        names(sect1) <- c( 
            "pipeline", "star.version", "samtools.version", 
            "species (genome)", "SE", "readLength")
        if (file.exists(filename_stats)) {
            x <- read.delim(filename_stats, head = F, as.is = T)
            sect2 <- c(
                as.numeric(x[grep("Number of input reads", x[, 1]), 2]), 
                as.numeric(x[grep("Average input read length", x[, 1]), 2]), 
                as.numeric(x[grep("Average mapped length", x[, 1]), 2]), 
                as.numeric(x[grep("Uniquely mapped reads number", x[, 1]), 2]), 
                as.numeric(x[grep("Number of reads mapped to multiple loci", x[, 1]), 2]) + as.numeric(x[grep("Uniquely mapped reads number", x[, 1]), 2]), 
                as.numeric(sub("%", "", x[grep("Uniquely mapped reads %", x[, 1]), 2])), 
                as.numeric(sub("%", "", x[grep("% of reads mapped to multiple loci", x[, 1]), 2])), 
                as.numeric(sub("%", "", x[grep("% of reads mapped to too many loci", x[, 1]), 2])) + as.numeric(sub("%", "", x[grep("% of reads unmapped: too short", x[, 1]), 2])) + as.numeric(sub("%", "", x[grep("% of reads unmapped: too many mismatches", x[, 1]), 2])) + as.numeric(sub("%", "", x[grep("% of reads unmapped: other", x[, 1]), 2])),
                as.numeric(sub("%", "", x[grep("Mismatch rate per base, %", x[, 1]), 2])),
                as.numeric(sub("%", "", x[grep("Deletion rate per base", x[, 1]), 2])),
                as.numeric(x[grep("Deletion average length", x[, 1]), 2]),
                as.numeric(sub("%", "", x[grep("Insertion rate per base", x[, 1]), 2])),
                as.numeric(x[grep("Insertion average length", x[, 1]), 2])
            )
        } else {
            sect2 <- rep(NA, 8 + 5)
        }
        names(sect2) <- c(
            "Tot Reads After Trim", 
            "Avg Inp Read Len", 
            "Avg Uniq Map Len",
            "Uniq-mapped Reads", 
            "Total-mapped Reads", 
            "Uniq-mapped Perc", 
            "Multi-mapped Perc", 
            "Not-mapped Perc", 
            "Mismatch-per-base Perc",
            "Deletion-per-base Perc",
            "Avg Deletion Length", 
            "Insertion-per-base Perc",
            "Avg Insertion Length"
            )
        c(sect1, sect2)
    }

    getVERSEstats_perFeature <- function(feature, exptName, sampleName, baseDir = "/lab/repo", analyzedDir = "analyzed", countsDir = "verse", fileBase = countsDir, .DEBUG = FALSE) {
        fnameFeature <- paste(sampleName, fileBase, feature, "cnts", "txt", sep = ".")
        filenameFeature <- file.path(baseDir, exptName, analyzedDir, sampleName, countsDir, fnameFeature)
        fnameSummary <- paste(sampleName, fileBase, "summary.txt", sep = ".")
        filenameSummary <- file.path(baseDir, exptName, analyzedDir, sampleName, countsDir, fnameSummary)
        if (file.exists(filenameSummary)) {
            if (.DEBUG) catLog("Reading", filenameSummary, "...\n", file = .LOGFILE)
            summary <- read.delim(filenameSummary, as.is = T, check = F, head = F)
        } else {
            catLog("Not found", filenameSummary, "...\n", file = .LOGFILE)
            summary <- data.frame()
        }
        valueHeaderMap <- c(
            "exon" = "exons Level 1,2", 
            "exon-lev3" = "exons Level 3", 
            "anti-exon" = "anti-exons", 
            "intron" = "introns Level 1,2",
            "intron-lev1-lev2" = "introns Level 1,2",
            "intron-lev3" = "introns Level 3",
            "anti-intron" = "anti-introns", 
            "mito" = "mito", 
            "anti-mito" = "anti-mito", 
            "intergenic" = "Intergenic")
        percHeaderMap <- c(
            "exon" = "Exons Level 1,2", 
            "exon-lev3" = "Exons Level 3", 
            "anti-exon" = "Anti-Exons", 
            "intron" = "Introns Level 1,2",
            "intron-lev1-lev2" = "Introns Level 1,2",
            "intron-lev3" = "Introns Level 3",
            "anti-intron" = "Anti-Introns", 
            "mito" = "Mito", 
            "anti-mito" = "Anti-Mito", 
            "intergenic" = "Intergenic")
        namesPercHeaderMap <- names(percHeaderMap)
        percHeaderMap <- paste0("Perc: ", percHeaderMap)
        names(percHeaderMap) <- namesPercHeaderMap
        stats <- NULL
        if (file.exists(filenameFeature)) {
            if (feature != "lines_sines") {
                if (.DEBUG) catLog("Reading", filenameFeature, "...\n", file = .LOGFILE)
                counts <- read.csv(filenameFeature, row.names = 1, as.is = T, check = F, sep = "\t")
                if (feature == "exon") {
                    endogenes <- grep("^spikeIn", rownames(counts), value = T, invert = T)
                    spikeins <- setdiff(rownames(counts), endogenes)
                    readsCounted <- sum(counts[endogenes, 1])
                    ## Bug: "VERSE.Perc: Exons Level 1,2" has negative values (2017-06-01)
                    ## - totalReads <- sum(counts[, 1])
                    ## +
                    totalReads <- as.numeric(summary[summary[, 1] == "TotalReadPairs", 2])
                    ## Fixed
                    spikeReadsCounted <- sum(counts[spikeins, 1])
                    pExon <- ifelse(prod(dim(summary)) == 0, NA, as.numeric(sub("%$", "", summary[grep(sprintf("^Assigned%sFraction", feature), summary[, 1], ignore.case = T), 2])))
                    if (length(spikeins) > 0) {
                        spikePercs <- as.numeric(ifelse(spikeReadsCounted > 0, sprintf("%.2f", spikeReadsCounted / totalReads * 100), 0))
                        pExon <- pExon - spikePercs
                    } else {
                        spikePercs <- NA
                    }
                    stat <- c(readsCounted, spikeReadsCounted, spikePercs, pExon)
                    names(stat) <- c(paste0(valueHeaderMap[feature], ": non-SpikeIn Reads Counted"), paste0(valueHeaderMap[feature], ": SpikeIn Reads Counted"), "Perc: Spike-In", paste0(percHeaderMap[feature]))
                    stats <- c(stats, stat)
                } else {
                    readsCounted <- ifelse(prod(dim(summary)) == 0, NA, as.numeric(summary[grep(sprintf("^Assigned%s$", feature), summary[, 1], ignore.case = T), 2]))
                    pFeature <- ifelse(prod(dim(summary)) == 0, NA, as.numeric(sub("%$", "", summary[grep(sprintf("^Assigned%sFraction", feature), summary[, 1], ignore.case = T), 2])))
                    stat <- c(readsCounted, pFeature)
                    names(stat) <- c(
                        paste0(valueHeaderMap[feature], ": Reads Counted"), 
                        percHeaderMap[feature])
                    stats <- c(stats, stat)
                }
                numGenes <- sum(counts[, 1] > 0)
                avgReadPerGene <- ifelse(numGenes > 0, readsCounted / numGenes, 0)
                maxReadsPerGene <- max(counts[, 1])
                ambiguousMapped <- ifelse(prod(dim(summary)) == 0, NA, as.numeric(summary[grep(sprintf("^Ambiguous%s$", feature), summary[, 1], ignore.case = T), 2]))
                stat <- c(numGenes, avgReadPerGene, maxReadsPerGene, ambiguousMapped)
                names(stat) <- c(
                    paste0(valueHeaderMap[feature], ifelse(feature == "intergenic", ": Num Regions", ": Num Genes")), 
                    paste0(valueHeaderMap[feature], ": Avg Read Per ", ifelse(feature == "intergenic", "Region", "Gene")), 
                    paste0(valueHeaderMap[feature], ": Max Reads Per ", ifelse(feature == "intergenic", "Region", "Gene")), 
                    paste0(valueHeaderMap[feature], ": Ambiguous Mapped"))
                stats <- c(stats, stat)
            } else {
                ## Hot fix to cope with something wrong with the file 
                ## E.201/analyzed/Sample_CM453/verse/Sample_CM453.verse.lines_sines.cnts.txt
                ## where line 1408690 and 1408691 have the same rowname
                ## The duplicate rownames in current version of VERSE 
                ## lines_sines counts are gone. 
                if (.DEBUG) catLog("Reading", filenameFeature, "...\n", file = .LOGFILE)
                counts <- read.csv(filenameFeature, as.is = T, check = F, sep = "\t")
                fnameLinesSinesSummary <- paste(sampleName, countsDir, "lines_sines", "summary.txt", sep = ".")
                filenameLinesSinesSummary <- file.path(baseDir, exptName, analyzedDir, sampleName, countsDir, fnameLinesSinesSummary)
                duplicatedRows <-duplicated(counts[, "gene"])
                counts <- counts[!duplicatedRows, ]
                linesGenes <- grep("LINE", counts[, "gene"])
                sinesGenes <- grep("SINE", counts[, "gene"])
                LINEreadsCounted <- sum(counts[linesGenes, "count"])
                SINEreadsCounted <- sum(counts[sinesGenes, "count"])
                numLINEs <- sum(counts[linesGenes, "count"] > 0)
                numSINEs <- sum(counts[sinesGenes, "count"] > 0)
                if (file.exists(filenameLinesSinesSummary)) {
                    if (.DEBUG) catLog("Reading", filenameLinesSinesSummary, "...\n", file = .LOGFILE)
                    linesSinesSummary <- read.delim(filenameLinesSinesSummary, as.is = T, check = F, head = F)
                    linesSinesNoFeature <- linesSinesSummary[linesSinesSummary[, 1] == "NoFeature", 2]
                } else {
                    catLog("Not found", filenameLinesSinesSummary, "\n", file = .LOGFILE)
                    linesSinesNoFeature <- NA
                }
                stat <- c(LINEreadsCounted, SINEreadsCounted, numLINEs, numSINEs, linesSinesNoFeature)
                names(stat) <- c("LINE Reads", "SINE Reads", "Num LINEs", "Num SINEs", "Neither LINEs nor SINEs")
                stats <- c(stats, stat)
            }
        } else {
            catLog("Not found", filenameFeature, "\n", file = .LOGFILE)
            if (feature != "lines_sines") {
                if (feature == "exon") {
                    stat <- rep(NA, 4)
                    names(stat) <- c(paste0(valueHeaderMap[feature], ": non-SpikeIn Reads Counted"), paste0(valueHeaderMap[feature], ": SpikeIn Reads Counted"), "Perc: Spike-In", paste0(percHeaderMap[feature]))
                    stats <- c(stats, stat)
                } else {
                    stat <- rep(NA, 2)
                    names(stat) <- c(
                        paste0(valueHeaderMap[feature], ": Reads Counted"),
                        percHeaderMap[feature])
                    stats <- c(stats, stat)
                }
                stat <- rep(NA, 4)
                names(stat) <- c(
                    paste0(valueHeaderMap[feature], ifelse(feature == "intergenic", ": Num Regions", ": Num Genes")), 
                    paste0(valueHeaderMap[feature], ": Avg Read Per ", ifelse(feature == "intergenic", "Region", "Gene")), 
                    paste0(valueHeaderMap[feature], ": Max Reads Per ", ifelse(feature == "intergenic", "Region", "Gene")), 
                    paste0(valueHeaderMap[feature], ": Ambiguous Mapped"))
                stats <- c(stats, stat)
            } else {
                stat <- rep(NA, 5)
                names(stat) <- c("LINE Reads", "SINE Reads", "Num LINEs", "Num SINEs", "Neither LINEs nor SINEs")
                stats <- c(stats, stat)
            }
        }
        stats
    }

    getVERSEstats <- function(exptName, sampleID, baseDir = "/lab/repo", analyzedDir = "analyzed", countsDir = "verse", fileBase = countsDir, features = c( "exon", "exon-lev3", "anti-exon", "intron-lev1-lev2", "intron-lev3", "intron", "anti-intron", "mito", "anti-mito", "intergenic", "lines_sines"), .DEBUG = FALSE) {
        sampleName <- paste("Sample", gsub(" ", "", sampleID), sep = "_")
        fname_version <- paste(sampleName, "versions", sep = ".")
        stats <- c()
        filename_version <- file.path(baseDir, exptName, analyzedDir, sampleName, countsDir, fname_version)
        if (file.exists(filename_version)) {
            if (.DEBUG) catLog("Reading", filename_version, "...\n", file = .LOGFILE)
            sect1 <- as.matrix(read.delim(filename_version, as.is = T, check = F))[1, -c(1, 2)]
        } else {
            catLog("Not found", filename_version, "\n", file = .LOGFILE)
            sect1 <- rep("", 7)
        }
        names(sect1) <- c(
            "verse.version", "transcriptome", "stranded",
            "ID_attribute", "introns", "intergenic", "lines-sines")
        sect2 <- lapply(features, FUN = getVERSEstats_perFeature, exptName = exptName, sampleName = sampleName, baseDir = baseDir, analyzedDir = analyzedDir, countsDir = countsDir, fileBase = fileBase, .DEBUG = .DEBUG)
        sect2 <- do.call(c, sect2)
        fnameSummary <- paste(sampleName, countsDir, "summary.txt", sep = ".")
        filenameSummary <- file.path(baseDir, exptName, analyzedDir, sampleName, countsDir, fnameSummary)
        if (file.exists(filenameSummary)) {
            if (.DEBUG) catLog("Reading", filenameSummary, "...\n", file = .LOGFILE)
            summary <- read.delim(filenameSummary, as.is = T, check = F, head = F)
            noFEATURE <- as.numeric(summary[summary[, 1] == "NoFeature", 2])
            missingMATES <- as.numeric(summary[summary[, 1] == "MissingMates", 2])
        } else {
            catLog("Not found", filenameSummary, "...\n", file = .LOGFILE)
            noFEATURE <- missingMATES <- NA
        }
        idcsIntron <- grep("^introns: ", names(sect2))
        idcsIntronL1L2 <- grep("^introns Level 1,2: ", names(sect2))
        if (all(is.na(sect2[idcsIntronL1L2]))) {
            sect2[idcsIntronL1L2] <- sect2[idcsIntron]
        }
        idxPercsIntron <- grep("^Perc: Introns$", names(sect2))
        idxPercsIntronL1L2 <- grep("^Perc: Introns Level 1,2", names(sect2))
        if (all(is.na(sect2[idxPercsIntronL1L2])) && length(idxPercsIntron) == length(idxPercsIntronL1L2)) {
            sect2[idxPercsIntronL1L2] <- sect2[idxPercsIntron]
        }
        if (all(is.na(sect2[idcsIntronL1L2])) && all(is.na(sect2[idxPercsIntronL1L2]))) {
            sect2 <- sect2[-c(idcsIntron, idxPercsIntron)]
        }
        sect3 <- c("No Feature" = noFEATURE, "Reads Missing Mates" = missingMATES)
        stats <- c(sect1, sect2, sect3)
        idcsPercs <- grep("^Perc: ", names(stats))
        stats <- c(stats[-idcsPercs], stats[idcsPercs][-1], stats[idcsPercs][1])
    }

    getBLASTstats <- function(exptName, sampleID, baseDir = "/lab/repo", analyzedDir = "analyzed", blastDir = "blast", fileBase = blastDir, .DEBUG = FALSE) {
        sampleName <- paste("Sample", gsub(" ", "", sampleID), sep = "_")
        fname_version <- paste(sampleName, "versions", sep = ".")
        filename_version <- file.path(baseDir, exptName, analyzedDir, sampleName, blastDir, fname_version)
        fname_stats <- paste(sampleName, fileBase, "stats.txt", sep = ".")
        filename_stats <- file.path(baseDir, exptName, analyzedDir, sampleName, blastDir, fname_stats)
        if (file.exists(filename_version)) {
            if (.DEBUG) catLog("Reading", filename_version, "...\n", file = .LOGFILE)
            sect1 <- as.matrix(read.delim(filename_version, as.is = T, check.names = F))
            if (ncol(sect1) < 9) {
                sect1 <- c(sect1[, c(3, 5, 6)], NA, NA, NA) 
            } else {
                sect1 <- sect1[, c(3, 5, 6:9)]
            }
        } else {
            if (.DEBUG) catLog("Not found", filename_version, "\n", file = .LOGFILE)
            sect1 <- rep(NA, 6)
        }
        names(sect1) <- c(
            "blastn.version", "parseBlast.py.version", "species",
            "readLength", "numReads", "kmer")
        fields2.1 <- c("numReads", "kmer") 
        fields2.2 <- c(
            "Total Hits", "Hits Not Counted", 
            "Perc Hits Not Target or SpikeIn", "Bacteria", "Fish", "Fly", 
            "Human", "Mouse", "Rat", "Yeast", "ERCC", "PhiX")
        fields2.3 <- c("TATAGTGAGT", "Perc Counted")
        if (file.exists(filename_stats)) {
            if (.DEBUG) catLog("Reading", filename_stats, "...\n", file = .LOGFILE)
            x <- readLines(filename_stats)
            n <- length(x)
            idx_numReads <- grep("^Num reads: ", x)
            idx_kMerSeq <- grep("^k-mer sequence", x)
            idx_kMerCnt <- grep("^k-mer count", x)
            sect2.1 <- c(
                ifelse(length(idx_numReads) == 1, sub("^Num reads: ", "", x[idx_numReads])), 
                ifelse(length(idx_kMerSeq) == 1, sub("k-mer sequence: ", "", x[idx_kMerSeq]), NA))
            sect2.2 <- strsplit(x[n], "\\s+")[[1]]
            names(sect2.2) <- strsplit(x[n-1], "\t")[[1]]
            sect2.2["Hits Not Target or ERCC"] <- sub("%$", "", sect2.2["Hits Not Target or ERCC"])
            sect2.2 <- sect2.2[fields2.2]
            idx_pCnt <- grep("^Number hits:", x)
            sect2.3 <- c(
                ifelse(length(idx_kMerCnt) == 1, sub("k-mer count: ", "", x[idx_kMerCnt]), NA), 
                ifelse(length(idx_pCnt) == 1, sub(".*\\( (.+)% \\)", "\\1", x[idx_pCnt]), NA))
            names(sect2.3) <- c(ifelse(length(idx_kMerSeq) == 1, sub("k-mer sequence: ", "", x[idx_kMerSeq]), "TATAGTGAGT"), 
                                "Perc Counted")
            sect2.3 <- sect2.3[fields2.3]
        } else {
            if (.DEBUG) catLog("Not found", filename_stats, "\n", file = .LOGFILE)
            sect2.1 <- rep(NA, 2)
            sect2.2 <- rep(NA, 12)
            sect2.3 <- rep(NA, 2)
        }
        names(sect2.1) <- fields2.1
        names(sect2.2) <- fields2.2
        names(sect2.3) <- fields2.3
        if (is.na(sect1["numReads"])) sect1["numReads"] <- sect2.1["numReads"]
        if (is.na(sect1["kmer"])) sect1["kmer"] <- sect2.1["kmer"]
        sect2 <- c(sect2.2, sect2.3)
        stats <- c(sect1, sect2)
    }

    getTRIMstats <- function(exptName, sampleID, baseDir = "/lab/repo", analyzedDir = "analyzed", trimDir = "trim", fileBase = trimDir, .DEBUG = FALSE) {
        sampleName <- paste("Sample", gsub(" ", "", sampleID), sep = "_")
        fname_version <- paste(sampleName, "versions", sep = ".")
        filename_version <- file.path(baseDir, exptName, analyzedDir, sampleName, trimDir, fname_version)
        fname_stats <- paste(sampleName, fileBase, "stats.txt", sep = ".")
        filename_stats <- file.path(baseDir, exptName, analyzedDir, sampleName, trimDir, fname_stats)
        fields1 <- c("trimReads.py.version", "minLen", "phredThresh", "removeN",
                     "numAT", "SE", "ContamFile")
        if (file.exists(filename_version)) {
            if (.DEBUG) catLog("Reading", filename_version, "...\n", file = .LOGFILE)
            x <- as.matrix(read.delim(filename_version, as.is = T, check.names = F))
            if (ncol(x) < 9) {
                sect1 <- c(x[1, 3], rep(NA, 6))
            } else {
                sect1 <- x[1, 3:9]
            }
        } else {
            if (.DEBUG) catLog("Not found", filename_version, "\n", file = .LOGFILE)
            sect1 <- rep(NA, 7)
        }
        names(sect1) <- fields1
        fields2 <- c("nTotalReadsBeforeTrim", 
            "nBothTrimmed", "nFirstTrimmed", "nSecondTrimmed", 
            "nBothDiscarded", "nFirstDiscarded", "nSecondDiscarded", 
            "Phred Threshold (first)", "Phred Threshold (second)", 
            "removed5N (first)", "removed5N (second)", 
            "removed3N (first)", "removed3N (second)", 
            "indexAdapter (first)", "indexAdapter (second)", 
            "univAdapterRC (first)", "univAdapterRC (second)", 
            "aRNAPrimer (first)", "aRNAPrimer (second)", 
            "aRNAPrimerRC (first)", "aRNAPrimerRC (second)", 
            "5 polyT (first)", "5 polyT (second)", 
            "3 polyA (first)", "3 polyA (second)", 
            "tivaXoligo (first)", "tivaXoligo (second)",
            "tivaXoligoRC (first)", "tivaXoligoRC (second)",
            "dT5oligo (first)", "dT5oligo (second)", 
            "dT5oligoRC (first)", "dT5oligoRC (second)", 
            "NexteraRead1ExternalAdapter (first)", "NexteraRead1ExternalAdapter (second)", 
            "NexteraRead2ExternalAdapter (second)", 
            "dT6oligo (second)")
        if (file.exists(filename_stats)) {
            if (.DEBUG) catLog("Reading", filename_stats, "...\n", file = .LOGFILE)
            x <- readLines(filename_stats)
            n <- length(x)
            sect2 <- strsplit(x[n], "\\s+")[[1]]
            names(sect2) <- strsplit(x[n-1], "\t")[[1]]
            idx_nTotalReadPairs <- grep("nTotalReadPairs", names(sect2))
            if (!is.na(sect2["nTotalReadPairs"])) names(sect2)[idx_nTotalReadPairs] <- "nTotalReadsBeforeTrim"
            sect2 <- sect2[fields2]
        } else {
            if (.DEBUG) catLog("Not found", filename_stats, "\n", file = .LOGFILE)
            sect2 <- rep(NA, 37)
        }
        names(sect2) <- fields2
        stats <- c(sect1, sect2)
    }
    
    getNGSstats <- function(exptID, exptName = NULL, sampleID, baseDir = "/lab/repo", analyzedDir = "analyzed", blastDir = "blast", trimDir = "trim", alignedDir = "star", countsDir = "verse", blastFileBase = blastDir, trimFileBase = trimDir, alignedFileBase = alignedDir, countsFileBase = countsDir, features = c("exon", "exon-lev3", "anti-exon", "intron-lev1-lev2", "intron-lev3", "intron", "anti-intron", "mito", "anti-mito", "intergenic", "lines_sines"), .DEBUG = FALSE) {
        if (missing(exptID) && is.null(exptName)) stop("You need to at least provide exptID or exptName!")
        if (!missing(exptID) && is.null(exptName)) exptName <- paste0("E.", exptID) 
        catLog("Processing", exptName, sampleID, "...\n", file = .LOGFILE)
        STARstats <- getSTARstats(exptName, sampleID, baseDir = baseDir, analyzedDir = analyzedDir, alignedDir = alignedDir, fileBase = alignedFileBase, .DEBUG = .DEBUG)
        VERSEstats <- getVERSEstats(exptName, sampleID, baseDir = baseDir, analyzedDir = analyzedDir, countsDir = countsDir, fileBase = countsFileBase, features = features, .DEBUG = .DEBUG)
        BLASTstats <- getBLASTstats(exptName, sampleID, baseDir = baseDir, analyzedDir = analyzedDir, blastDir = blastDir, fileBase = blastFileBase, .DEBUG = .DEBUG)
        TRIMstats <- getTRIMstats(exptName, sampleID, baseDir = baseDir, analyzedDir = analyzedDir, trimDir = trimDir, fileBase = trimFileBase, .DEBUG = .DEBUG)
        names(STARstats)[-1] <- paste0("STAR.", names(STARstats)[-1])
        names(VERSEstats) <- paste0("VERSE.", names(VERSEstats))
        names(BLASTstats) <- paste0("BLAST.", names(BLASTstats))
        names(TRIMstats) <- paste0("TRIM.", names(TRIMstats))
        stats <- c(STARstats, VERSEstats, BLASTstats, TRIMstats)
    }

    for (obj in ls())
        assign(obj, get(obj), envir = Repo)
})


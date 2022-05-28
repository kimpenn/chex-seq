## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.

if (!exists("CHEX") || is.environment(CHEX)) {
    CHEX <- new.env(parent = emptyenv())
}

local({
    .VERSION = "1.5"

    createBioGroupConfig <- function(bioGroups, bioGroup2species, bioGroup2cellType, bioGroup2harvest, bioGroup2drugGroup, bioGroup2hasTPA, bioGroup2TPADose, bioGroup2TPATime, bioGroup2hasMungBean, outFile = "Data/20190601/BioGroupConfigTemp.csv") {
        species <- bioGroup2species[bioGroups]
        cellTypes <- bioGroup2cellType[bioGroups]
        harvests <- bioGroup2harvest[bioGroups]
        drugGroups <- bioGroup2drugGroup[bioGroups]
        hasTPAs <- bioGroup2hasTPA[bioGroups]
        TPADoses <- bioGroup2TPADose[bioGroups]
        TPATimes <- bioGroup2TPATime[bioGroups]
        hasMungBean <- bioGroup2hasMungBean[bioGroups]
        df <- data.frame(
            BioGroup = bioGroups, 
            Species = species, 
            CellType = cellTypes, 
            Harvest = harvests, 
            DrugGroup = drugGroups,
            HasTPA = hasTPAs,
            TPADose = TPADoses,
            TPATime = TPATimes,
            HasMungBean = hasMungBean,
            stringsAsFactors = FALSE)
        write.csv(df, file = outFile, row.names = FALSE)
    }

    ## ProbeType describes whether a sample has T7 probe or not, while IsNegCtrl whether the sample is negative or not
    ## IsNegCtrl denotes the overall case/control status; it is determined by two factors:
    ## ProbeType and BioType (CellNum, HasMungBean)
    createSampleInfoFull <- function(sampleInfoFile = "Data/20190601/SampleInfo20190625CV2b.csv", sampleInfoFullFile = "Data/20190601/SampleInfoFull20190625CV2b.csv", bioGroups = c("K562", "HumanAstroCulture", "HumanNeuronCulture", "MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice"), bioGroupConfigFile = "Data/20190601/BioGroupConfig20190625CV2b.csv",  sampleIDsOutFile = NULL, sampleIDsToExclude = NULL, exptIDsToExclude = NULL) {
        if ((!is.null(sampleIDsOutFile)) && file.exists(sampleIDsOutFile)) {
            sampleIDsOut <- as.character(readLines(sampleIDsOutFile))
        } else {
            sampleIDsOut <- character(0)
        }
        
        BioGroupConfig <- read.csv(bioGroupConfigFile, as.is = TRUE, check.names = FALSE)
        rownames(BioGroupConfig) <- BioGroupConfig$BioGroup

        SampleInfo <- read.csv(sampleInfoFile, as.is = TRUE, check.names = FALSE)
        if (!is.null(sampleIDsToExclude)) { SampleInfo <- SampleInfo[!SampleInfo[["SampleID"]] %in% sampleIDsToExclude, ] }
        if (!is.null(exptIDsToExclude)) { SampleInfo <- SampleInfo[!SampleInfo[["ExptID"]] %in% exptIDsToExclude, ] }
        ## HasProbe == "N" if ProbeID == "none", "Y" otherwise
        ## ProbeType is "Positive" if HasProbe == "Y" AND HasLaser405 == "Y"
        ## ProbeType is "NegCtrl" if HasProbe == "N" OR HasLaser405 == "N"
        SampleInfo$HasProbe <- ifelse(SampleInfo$ProbeID == "none", "N", "Y")
        SampleInfo$ProbeType <- with(SampleInfo, ifelse(HasProbe == "Y" & HasLaser405 == "Y", "Positive", "NegCtrl"))
        sampleIDs <- SampleInfo$SampleID
        rownames(SampleInfo) <- sampleIDs

        if (length(sampleIDsOut) > 0) {
            if (!all(sampleIDsOut %in% sampleIDs)) {
                warning("Some outliers are not found in SampleInfo!")
            }
            sampleIDs <- setdiff(sampleIDs, sampleIDsOut)
            SampleInfo <- SampleInfo[sampleIDs, ]
        }

        sampleIDsLists <- sapply(bioGroups, function(bioGroup) {
            species <- BioGroupConfig[bioGroup, "Species"]
            cellType <- BioGroupConfig[bioGroup, "CellType"]
            harvest <- BioGroupConfig[bioGroup, "Harvest"]
            drugGroup <- BioGroupConfig[bioGroup, "DrugGroup"]
            hasProbe <- BioGroupConfig[bioGroup, "HasProbe"] ## has probe or not
            hasTPA <- BioGroupConfig[bioGroup, "HasTPA"]
            TPAdose <- BioGroupConfig[bioGroup, "TPADose"]
            TPAtime <- BioGroupConfig[bioGroup, "TPATime"]
            hasMungBean <- BioGroupConfig[bioGroup, "HasMungBean"] ## has single-strand chromatin digested or not
            ## Note, some samples have CellNum NA, but we are sure that they have multiple cells included.
            ## Suffix {Positive, NegCtrl} denotes whether a sample got positive treatment (HasProbe AND HasLaser405) or not (!HasProbe OR !HasLaser405)
            ## Here "Positive" and "NegCtrl" denote technical case/control type. To be *technically positive*, a sample must satisfy ProbeType == "Positive" AND HasLaser405 == "Y"
            ## Hence, the suffix is identical to ProbeType
            ## Note, a techincal positive sample can be biologically negative, for example, no cell or cells with mung bean digestion
            SampleInfoPositiveNone <- subset(SampleInfo, Species == species & CellType == cellType & Harvest == harvest & DrugGroup == drugGroup & HasTPA == hasTPA & TPADose == TPAdose & TPATime == TPAtime & HasMungBean == hasMungBean & ProbeType == "Positive" & CellNum == 0)
            SampleInfoNegCtrlNone <- subset(SampleInfo, Species == species & CellType == cellType & Harvest == harvest & DrugGroup == drugGroup & HasTPA == hasTPA & TPADose == TPAdose & TPATime == TPAtime & HasMungBean == hasMungBean & ProbeType == "NegCtrl" & CellNum == 0)
            SampleInfoPositiveSingle <- subset(SampleInfo, Species == species & CellType == cellType & Harvest == harvest & DrugGroup == drugGroup & HasTPA == hasTPA & TPADose == TPAdose & TPATime == TPAtime & HasMungBean == hasMungBean & ProbeType == "Positive" & CellNum == 1)
            SampleInfoNegCtrlSingle <- subset(SampleInfo, Species == species & CellType == cellType & Harvest == harvest & DrugGroup == drugGroup & HasTPA == hasTPA & TPADose == TPAdose & TPATime == TPAtime & HasMungBean == hasMungBean & ProbeType == "NegCtrl" & CellNum == 1)
            SampleInfoPositivePooled <- subset(SampleInfo, Species == species & CellType == cellType & Harvest == harvest & DrugGroup == drugGroup & HasTPA == hasTPA & TPADose == TPAdose & TPATime == TPAtime & HasMungBean == hasMungBean & ProbeType == "Positive" & (is.na(CellNum) | CellNum > 1))
            SampleInfoNegCtrlPooled <- subset(SampleInfo, Species == species & CellType == cellType & Harvest == harvest & DrugGroup == drugGroup & HasTPA == hasTPA & TPADose == TPAdose & TPATime == TPAtime & HasMungBean == hasMungBean & ProbeType == "NegCtrl" & (is.na(CellNum) | CellNum > 1))

            sampleIDsPositiveNone <- SampleInfoPositiveNone[, "SampleID"]
            sampleIDsNegCtrlNone <- SampleInfoNegCtrlNone[, "SampleID"]
            sampleIDsPositiveSingle <- SampleInfoPositiveSingle[, "SampleID"]
            sampleIDsNegCtrlSingle <- SampleInfoNegCtrlSingle[, "SampleID"]
            sampleIDsPositivePooled <- SampleInfoPositivePooled[, "SampleID"]
            sampleIDsNegCtrlPooled <- SampleInfoNegCtrlPooled[, "SampleID"]
            sampleIDsNegCtrlNone <- SampleInfoNegCtrlNone[, "SampleID"]

            ## Biological sample IDs
            sampleIDsListBiol <- list(
                PositiveNone = sampleIDsPositiveNone, NegCtrlNone = sampleIDsNegCtrlNone, 
                PositiveSingle = sampleIDsPositiveSingle, NegCtrlSingle = sampleIDsNegCtrlSingle, 
                PositivePooled = sampleIDsPositivePooled, NegCtrlPooled = sampleIDsNegCtrlPooled)

            ## Virtual sample IDs by merging no-cell samples, exclusive to NegCtrl: "Nil" 
            sampleIDsListNil <- list(
                PositiveNil = if (length(sampleIDsPositiveNone) > 1) { paste0(bioGroup, "PositiveNil") } else { character(0) },
                NegCtrlNil = if (length(sampleIDsNegCtrlNone) > 1) { paste0(bioGroup, "NegCtrlNil") } else { character(0) })
            sourceIDsListNil <- list(
                PositiveNil = if (length(sampleIDsPositiveNone) > 1) { sampleIDsPositiveNone } else { character(0) },
                NegCtrlNil = if (length(sampleIDsNegCtrlNone) > 1) { sampleIDsNegCtrlNone } else { character(0) })

            ## Virtual sample IDs by merging single cells only: "Merged"
            sampleIDsListMerged <- list(
                PositiveMerged = if (length(sampleIDsPositiveSingle) > 1) { paste0(bioGroup, "PositiveMerged") } else { character(0) },
                NegCtrlMerged = if (length(sampleIDsNegCtrlSingle) > 1) { paste0(bioGroup, "NegCtrlMerged") } else { character(0) })
            sourceIDsListMerged <- list(
                PositiveMerged = if (length(sampleIDsPositiveSingle) > 1) { sampleIDsPositiveSingle } else { character(0) },
                NegCtrlMerged = if (length(sampleIDsNegCtrlSingle) > 1) { sampleIDsNegCtrlSingle } else { character(0) })

            ## Virtual sample IDs by merging both single cells and pooled cells or multiple pooled cells; if no pooled cells but only single cells, no virtual sample ID will be designated
            sampleIDsListAll <- list(
                PositiveAll = if ((length(sampleIDsPositiveSingle) > 0 & length(sampleIDsPositivePooled) > 0) || (length(sampleIDsPositivePooled) > 1)) { paste0(bioGroup, "PositiveAll") } else { character(0) },
                NegCtrlAll = if ((length(sampleIDsNegCtrlSingle) > 0 & length(sampleIDsNegCtrlPooled) > 0) || (length(sampleIDsNegCtrlPooled) > 1)) { paste0(bioGroup, "NegCtrlAll") } else { character(0) }
            )
            sourceIDsListAll <- list(
                PositiveAll = if ((length(sampleIDsPositiveSingle) > 0 & length(sampleIDsPositivePooled) > 0) || (length(sampleIDsPositivePooled) > 1)) { c(sampleIDsPositiveSingle, sampleIDsPositivePooled) } else { character(0) },
                NegCtrlAll = if ((length(sampleIDsNegCtrlSingle) > 0 & length(sampleIDsNegCtrlPooled) > 0) || (length(sampleIDsNegCtrlPooled) > 1)) { c(sampleIDsNegCtrlSingle, sampleIDsNegCtrlPooled) } else { character(0) }
            )

            ## Full sample IDs are a union of the biological (""), no-cell ("Nil"), merged ("Merged"), and all ("All")
            sampleIDsListFull <- c(sampleIDsListBiol, sampleIDsListNil, sampleIDsListMerged, sampleIDsListAll)

            ## Virtual sample IDs are a union of the merged ("Merged") and all ("All")
            sourceIDsListVirtual <- c(sourceIDsListNil, sourceIDsListMerged, sourceIDsListAll)
            names(sourceIDsListVirtual) <- paste0(bioGroup, names(sourceIDsListVirtual))

            list(sampleIDsListBiol = sampleIDsListBiol, 
                 sampleIDsListNil = sampleIDsListNil, sourceIDsListNil = sourceIDsListNil,
                 sampleIDsListMerged = sampleIDsListMerged, sourceIDsListMerged = sourceIDsListMerged, 
                 sampleIDsListAll = sampleIDsListAll, sourceIDsListAll = sourceIDsListAll, 
                 sampleIDsListFull = sampleIDsListFull, sourceIDsListVirtual = sourceIDsListVirtual)
        }, simplify = FALSE)

        sampleIDsListBiol <- sapply(bioGroups, function(bioGroup) sampleIDsLists[[bioGroup]][["sampleIDsListBiol"]], simplify = FALSE)
        sampleIDsListNil <- sapply(bioGroups, function(bioGroup) sampleIDsLists[[bioGroup]][["sampleIDsListNil"]], simplify = FALSE)
        sampleIDsListMerged <- sapply(bioGroups, function(bioGroup) sampleIDsLists[[bioGroup]][["sampleIDsListMerged"]], simplify = FALSE)
        sampleIDsListAll <- sapply(bioGroups, function(bioGroup) sampleIDsLists[[bioGroup]][["sampleIDsListAll"]], simplify = FALSE)
        sampleIDsListFull <- sapply(bioGroups, function(bioGroup) sampleIDsLists[[bioGroup]][["sampleIDsListFull"]], simplify = FALSE)
        sourceIDsListVirtual <- sapply(bioGroups, function(bioGroup) sampleIDsLists[[bioGroup]][["sourceIDsListVirtual"]], simplify = FALSE)

        ## Append 3 additional columns to biological SampleInfo
        SampleInfo$CompType <- "Biol"
        SampleInfo$Composition <- with(SampleInfo, ifelse(is.na(CellNum) | CellNum > 1, "Pooled", ifelse(CellNum == 1, "Single", "None")))
        SampleInfo$SourceIDs <- SampleInfo$SampleID

        compositions <- c("Nil", "Merged", "All")

        ## Derive additional columns for virtual samples
        SampleInfoVirtual <- sapply(compositions, function(composition) {
            sampleIDsListName <- paste0("sampleIDsList", composition)
            sampleIDsList <- sapply(bioGroups, function(bioGroup) sampleIDsLists[[bioGroup]][[sampleIDsListName]], simplify = FALSE)
            bio.TrtGroups <- names(unlist(sampleIDsList))
            if (length(bio.TrtGroups) == 0) {
                SampleInfo <- data.frame(
                    GroupName = character(0), ExptID = character(0), 
                    SampleID = character(0), SeqType = character(0), 
                    Species = character(0),
                    EndType = character(0), Length = integer(0), Stranded = logical(0), Spikein = character(0), Contan = character(0), 
                    Barcoded = logical(0), Barcode = character(0), ProbesDegen = character(0), "2p_fwID_search" = character(0), Protocol = character(0),
                    ProbeID = character(0),
                    "2p_fwID_actual" = character(0), pc_fwID = integer(0), pcnoc_fwID = integer(0), "1p_fwID" = integer(0), "#PolyC" = integer(0),
                    CellType = character(0), 
                    HasProbe = character(0), 
                    HasLaser405 = character(0),
                    HasLaser633 = character(0), CellNum = character(0), HasLPA = character(0), 
                    IsNegCtrl = character(0),
                    ProbeType = character(0), 
                    Harvest = character(0), 
                    DrugGroup = character(0), 
                    HasTPA = character(0), 
                    TPADose = character(0), 
                    TPATime = character(0), 
                    HasMungBean = character(0), 
                    CompType = character(0), Composition = character(0), 
                    SourceIDs = character(0), 
                    stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL
                )
            } else {
                bioTrtGroups <- gsub("\\.", "",  bio.TrtGroups)
                bioGroups <- gsub("\\..*$", "", bio.TrtGroups)
                trtGroups <- gsub("^.*\\.", "", bio.TrtGroups)
                SampleInfo <- data.frame(
                    GroupName = paste0("CV2b", bioTrtGroups), ExptID = "-", 
                    SampleID = bioTrtGroups, SeqType = "CHEXSeq", 
                    Species = BioGroupConfig[bioGroups, "Species"],
                    EndType = "PE", Length = 75, Stranded = FALSE, Spikein = "none", Contan = "Default", 
                    Barcoded = TRUE, Barcode = "-", ProbesDegen = "-", "2p_fwID_search" = "-", Protocol = "C_V2b",
                    ProbeID = "-", #ifelse(grepl("Positive", trtGroups), "-", "none"),
                    "2p_fwID_actual" = "-", pc_fwID = 302, pcnoc_fwID = 303, "1p_fwID" = 303, "#PolyC" = 13,
                    CellType = BioGroupConfig[bioGroups, "CellType"], 
                    HasProbe = "-", # ifelse(grepl("Positive", trtGroups), "Y", "N"), 
                    HasLaser405 = "-", HasLaser633 = "-", CellNum = "-", HasLPA = "N", 
                    IsNegCtrl = "-", 
                    ProbeType = ifelse(grepl("Positive", trtGroups), "Positive", "NegCtrl"), 
                    Harvest = BioGroupConfig[bioGroups, "Harvest"], 
                    DrugGroup = BioGroupConfig[bioGroups, "DrugGroup"], 
                    HasTPA = BioGroupConfig[bioGroups, "HasTPA"], 
                    TPADose = BioGroupConfig[bioGroups, "TPADose"], 
                    TPATime = BioGroupConfig[bioGroups, "TPATime"], 
                    HasMungBean = BioGroupConfig[bioGroups, "HasMungBean"], 
                    CompType = "Virtual", Composition = composition, 
                    SourceIDs = mapply(function(bioGroup, bioTrtGroup) { paste(unlist(sourceIDsListVirtual[[bioGroup]][[bioTrtGroup]]), collapse = ",") }, bioGroup = bioGroups, bioTrtGroup = bioTrtGroups), 
                    stringsAsFactors = FALSE, check.names = FALSE, row.names = bioTrtGroups
                )
            }
            SampleInfo
        }, simplify = FALSE)

        SampleInfoFull <- rbind(SampleInfo, 
                                SampleInfoVirtual[["Nil"]], 
                                SampleInfoVirtual[["Merged"]], 
                                SampleInfoVirtual[["All"]])
        sampleIDsFull0 <- SampleInfoFull$SampleID
        sampleIDsFull <- unname(unlist(sampleIDsListFull))
        SampleInfoFull <- SampleInfoFull[sampleIDsFull, ]

        bio.TrtGroupsFull <- sub("[0-9]*$", "", names(unlist(sampleIDsListFull)))
        bioTrtGroupsFull <- gsub("\\.", "",  bio.TrtGroupsFull)
        bioGroupsFull <- gsub("\\..*$", "", bio.TrtGroupsFull)
        trtGroupsFull <- gsub("^.*\\.", "", bio.TrtGroupsFull)

        ## Append additional columns to full sample sheet
        SampleInfoFull$BioGroup <- bioGroupsFull
        SampleInfoFull$TrtGroup <- trtGroupsFull
        SampleInfoFull$BioTrtGroup <- bioTrtGroupsFull

        ## Fill in CellNum and IsNegCtrl for virtual samples
        SampleInfoFull$CellNum <- with(SampleInfoFull, sapply(SampleID, function(sampleID) { sourceIDs <- strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]]; cellNums <- as.numeric(SampleInfoFull[sourceIDs, "CellNum"]); sum(cellNums) }))
        SampleInfoFull$ProbeID <- with(SampleInfoFull, sapply(SampleID, function(sampleID) { ifelse(SampleInfoFull[sampleID, "CompType"] == "Biol", SampleInfoFull[sampleID, "ProbeID"], ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "ProbeID"] == "none"), "none", "-")) }))
        SampleInfoFull$HasProbe <- with(SampleInfoFull, sapply(SampleID, function(sampleID) { ifelse(SampleInfoFull[sampleID, "CompType"] == "Biol", SampleInfoFull[sampleID, "HasProbe"], ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "HasProbe"] == "Y"), "Y", ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "HasProbe"] == "N"), "N", "-"))) }))
        SampleInfoFull$HasLaser633 <- with(SampleInfoFull, sapply(SampleID, function(sampleID) { ifelse(SampleInfoFull[sampleID, "CompType"] == "Biol", SampleInfoFull[sampleID, "HasLaser633"], ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "HasLaser633"] == "Y"), "Y", ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "HasLaser633"] == "N"), "N", "-"))) }))
        SampleInfoFull$HasLaser405 <- with(SampleInfoFull, sapply(SampleID, function(sampleID) { ifelse(SampleInfoFull[sampleID, "CompType"] == "Biol", SampleInfoFull[sampleID, "HasLaser405"], ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "HasLaser405"] == "Y"), "Y", ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "HasLaser405"] == "N"), "N", "-"))) }))
        SampleInfoFull$IsNegCtrl <- with(SampleInfoFull, sapply(SampleID, function(sampleID) { ifelse(SampleInfoFull[sampleID, "CompType"] == "Biol", SampleInfoFull[sampleID, "IsNegCtrl"], ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "IsNegCtrl"] == "Y"), "Y", ifelse(all(SampleInfo[strsplit(SampleInfoFull[sampleID, "SourceIDs"], ",")[[1]], "IsNegCtrl"] == "N"), "N", "-"))) }))

        ## Append SampleName column
        compMap <- c(None = "n", Single = "s", Pooled = "p", Nil = "z", Merged = "m", All = "a")
        SampleInfoFull$SampleName <- with(SampleInfoFull, paste0(compMap[Composition], "CHEX", sprintf("%003i", seq(nrow(SampleInfoFull)))))

        SampleInfoFull <- SampleInfoFull[sampleIDsFull0, ]
        write.csv(SampleInfoFull, file = sampleInfoFullFile, row.names = FALSE)
        invisible(SampleInfoFull)
    }
    
    readQualStats <- function(readQualCategoryStatsFile) {
        QualStats <- read.csv(readQualCategoryStatsFile, sep = "\t", as.is = TRUE, check.names = FALSE, head = FALSE)
        colnames(QualStats) <- QualStats[1, ]
        QualStats <- QualStats[-1, -17]
        for (i in 4:16) QualStats[, i] <- as.numeric(QualStats[, i])
        QualStats
    }

    makeWideTableQualStatsPE <- function(QualStats) {
        QualStatsList <- list()
        sampleIDs <- subset(QualStats, Qual == "D" & OutputInPair == "either")[, "SampleID"]
        for (qual in paste(rep(c("A", "B", "C"), each = 2), 1:2, sep = "")) {
            for (outputInPair in c("read", "mate", "both", "either")) {
                Stats <- subset(QualStats, Qual == qual & OutputInPair == outputInPair)    
                if (outputInPair == "read") {
                    QualStatsList[[qual]][["ReadSingleton"]] <- data.frame(ReadSingleton = Stats[, "Singletons"])
                } else if (outputInPair == "mate") {
                    QualStatsList[[qual]][["MateSingleton"]] <- data.frame(MateSingleton = Stats[, "Singletons"])
                } else if (outputInPair == "both") {
                    QualStatsList[[qual]][["ReadMatePairs"]] <- data.frame(ReadMatePairs = Stats[, "Total"])
                } else if (outputInPair == "either") {
                    QualStatsList[[qual]][["Total"]] <- data.frame(Total = Stats[, "Total"])
                }
            }
        }
        QualStatsDf <- lapply(QualStatsList, function(X) do.call(cbind, X))
        QualStatsDf <- do.call(cbind, QualStatsDf)
        rownames(QualStatsDf) <- sampleIDs
        QualStatsDf <- data.frame(QualStatsDf, "D.R1Singleton" = with(subset(QualStats, Qual == "D" & OutputInPair == "either"), Read1 - `Properly paired`/2))
        QualStatsDf <- data.frame(QualStatsDf, "D.R2Singleton" = with(subset(QualStats, Qual == "D" & OutputInPair == "either"), Read2 - `Properly paired`/2))
        QualStatsDf <- data.frame(QualStatsDf, "D.ReadMatePairs" = with(subset(QualStats, Qual == "D" & OutputInPair == "either"), `Properly paired`))
        QualStatsDf <- data.frame(QualStatsDf, "D.Total" = subset(QualStats, Qual == "D" & OutputInPair == "either")[, "Total"])
    }

    makeWideTableQualStatsSE <- function(QualStats) {
        QualStatsList <- list()
        sampleIDs <- subset(QualStats, Qual == "D" & OutputInPair == "read")[, "SampleID"]
        for (qual in paste(rep(c("A", "B", "C"), each = 2), 1:2, sep = "")) {
            for (outputInPair in "read") {
                Stats <- subset(QualStats, Qual == qual & OutputInPair == outputInPair)    
                QualStatsList <- data.frame(ReadSingleton = Stats[, "Total"])
                QualStatsList <- data.frame(MateSingleton = integer(0))
                QualStatsList <- data.frame(ReadMatePairs = integer(0))
                QualStatsList <- data.frame(Total = Stats[, "Total"])
            }
        }
        QualStatsDf <- lapply(QualStatsList, function(X) do.call(cbind, X))
        QualStatsDf <- do.call(cbind, QualStatsDf)
        rownames(QualStatsDf) <- sampleIDs
        QualStatsDf <- data.frame(QualStatsDf, "D.Total" = subset(QualStats, Qual == "D" & OutputInPair == "read")[, "Total"])
    }

## Reorganize the long table into wide format:
## Example:
## sampleID, Qual, ReadSingletons, MateSingletons, PairedReads, Total
## sample_PE_#1, A1, 14, 24, 72, 100
## sample_PE_#1, A2, 14, 24, 72, 100
## sample_PE_#1, B1, 14, 24, 72, 100
## ...
## sample_PE_#1, D, NA, NA, NA, 100 ## For PE data, D class, there is no definition of singleton read/mate or read pairs
## sample_SE_#2, B1, 100, 0, 0, 100 ## For SE data, all reads are singletons.
## sample_SE_#2, B2, 100, 0, 0, 100 
## sample_SE_#2, C1, 100, 0, 0, 100 
## sample_SE_#2, C2, 100, 0, 0, 100 
## sample_SE_#2, D, 100, 0, 0, 100 
## Note, here all numbers are by reads (not read pairs), hence they should sum up to the "Total reads" of `samtools flagstats`
    makeWideTableQualStats <- function(QualStats, SampleInfo) {
        sampleIDsPE <- subset(SampleInfo, EndType == "PE")[, "SampleID"]
        sampleIDsSE <- subset(SampleInfo, EndType == "SE")[, "SampleID"] 
        QualStatsPE <- subset(QualStats, SampleID %in% sampleIDsPE)
        QualStatsSE <- subset(QualStats, SampleID %in% sampleIDsSE)
        QualStatsDfPE <- makeWideTableQualStatsPE(QualStatsPE)
        QualStatsDfSE <- makeWideTableQualStatsPE(QualStatsSE)
        QualStatsDf <- rbind(QualStatsDfPE, QualStatsDfSE)
        QualStatsDf$Total <- rowSums(with(QualStatsDf, cbind(A1.Total, A2.Total, B1.Total, B2.Total, C1.Total, C2.Total, D.Total)))
        QualStatsDf
    }

## Remove {A, B, C} x {1,2} x Total
    makeLongTableQualStats <- function(QualStats, SampleInfo, bioGroups = c("K562", "HumanAstroCulture", "HumanNeuronCulture", "MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "RatCardiomyoCulture"), trtGroups = c("NegCtrlSingle", "NegCtrlPooled", "PositiveSingle", "PositivePooled")) {
        library("reshape2")
        D.Total <- QualStats[, "D.Total"]
        QualStats <- QualStats[, !grepl("Total", colnames(QualStats))]
        QualStats <- cbind(QualStats, D.Total = D.Total)
        QualStatsLDf <- melt(QualStats)
        QualStats[, !grepl("Total", colnames(QualStats))]
        QualPairType <- strsplit(as.character(QualStatsLDf$variable), "\\.")
        QualStatsLDf$Qual <- sapply(QualPairType, "[", 1)
        QualStatsLDf$ReadPairType <- sapply(QualPairType, "[", 2)
        QualStatsLDf$variable <- NULL
        QualStatsLDf <- with(QualStatsLDf, data.frame(SampleID = SampleID, Qual = Qual, ReadPairType = ReadPairType, value = value, stringsAsFactors = FALSE))
        rownames(QualStatsLDf) <- NULL
        Df <- QualStatsLDf
        sampleIDs <- Df[, "SampleID"]
        SampleGroups <- SampleInfo[sampleIDs, c("BioGroup", "TrtGroup", "BioTrtGroup", "SampleName")]
        Df <- cbind(Df, SampleGroups, ExptID = SampleInfo[sampleIDs, "ExptID"])
        rownames(Df) <- NULL
        Df$BioGroup <- with(Df, factor(BioGroup, levels = bioGroups))
        Df$TrtGroup <- with(Df, factor(TrtGroup, levels = trtGroups))
        Df$BioTrtGroup <- with(Df, factor(BioTrtGroup, levels = unique(BioTrtGroup)))
        Df$ReadPairType <- with(Df, factor(ReadPairType, levels = c("ReadSingleton", "MateSingleton", "ReadMatePairs", "Total")))
        Df
    }

    for (obj in ls()) {
        assign(obj, get(obj), envir = CHEX)
    }
})

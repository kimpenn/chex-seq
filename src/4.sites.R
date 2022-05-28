source("src/lib/Genome.R")

SampleInfo <- read.csv("data/SampleInfo.csv", as.is = TRUE, check.names = FALSE)
rownames(SampleInfo) <- sampleIDs <- SampleInfo[, "SampleID"]

qualOutInPairs <- c(
    "Aread5End", "AreadAndFrag", "AbothFrag", 
    "Bread5End", "BreadAndFrag", "BbothFrag", 
    "Cmate5End", "CmateReadAndFrag", "CbothFrag", 
    "Deither5End", "DeitherReadAndFrag", "DbothFrag"
)
qualOutInPairNoAdjBases <- c(
    "Aeither.NoAdj", "Aeither.NoAdj", "Aeither.NoAdj", 
    "Beither.NoAdj", "Beither.NoAdj", "Beither.NoAdj", 
    "Ceither.NoAdj", "Ceither.NoAdj", "Ceither.NoAdj", 
    "Deither", "Deither", "Deither"
)
outTypes <- c(
    "5End", "ReadAndFrag", "Frag", 
    "5End", "ReadAndFrag", "Frag", 
    "5End", "ReadAndFrag", "Frag", 
    "5End", "ReadAndFrag", "Frag"
)
QualOutInPairConfig <- data.frame(
    QualOutInPair = qualOutInPairs, 
    QualOutInPairNoAdjBase = qualOutInPairNoAdjBases, 
    OutType = outTypes, 
    stringsAsFactors = FALSE
)

GRsQualOutInPairList <- lapply(1:nrow(QualOutInPairConfig), function(i) {
    qualOutInPair <- QualOutInPairConfig[i, "QualOutInPair"]
    qualOutInPairNoAdjBase <- QualOutInPairConfig[i, "QualOutInPairNoAdjBase"]
    outType <- QualOutInPairConfig[i, "OutType"]
    GRs <- sapply(sampleIDs, function(sampleID) {
        species <- SampleInfo[sampleID, "Species"]
        ## We treat the water samples as mouse
        if (species == "none") {
            species <- "mouse"
        }
        SeqInfo <- Genome$MainSeqInfo[[species]]
        seqLevels <- Genome$MainSeqLevels[[species]]
        seqLengths <- Genome$MainSeqLengths[[species]]
        sampleName <- paste0("Sample_", sampleID)
        filename <- sprintf("data/E.chex/analyzed/%s/star/%s.star.primaryNoDup.%s.%s.NoBlacklisted.bed", sampleName, sampleName, qualOutInPairNoAdjBase, outType)
        message(filename)
        Genome$getGRsFromBed(filename, SeqInfo = SeqInfo, seqLevels = seqLevels, seqLengths = seqLengths)
    }, simplify = FALSE)
})
names(GRsQualOutInPairList) <- qualOutInPairs

dir.create("results/PrimingSites", FALSE, TRUE)
## For each barcode/primer class, save the priming site data.
for (qualOutInPair in qualOutInPairs) {
    GRs <- GRsQualOutInPairList[[qualOutInPair]]
    filename <- sprintf("results/PrimingSites/GRs%s.RDS", qualOutInPair)
    message(filename)
    saveRDS(GRs, file = filename)
}

## Load the data
GRsQualOutInPairList <- sapply(qualOutInPairs, function(qualOutInPair) {
    message(qualOutInPair)
    filename <- sprintf("results/PrimingSites/GRs%s.RDS", qualOutInPair)
    readRDS(filename)
}, simplify = FALSE)

## Filter out (1) contaminants, (2) poor-mapping-quality reads
contamReadIDs <- sapply(sampleIDs, function(sampleID) {
    sampleName <- paste0("Sample_", sampleID)
    filename <- sprintf("data/E.chex/analyzed/%s/starAlt/%s.star.aligned.primaryNoDup.contam.readID.txt.gz", sampleName, sampleName)
    message(filename)
    readLines(filename)
}, simplify = FALSE)
dir.create("results/AlignmentArtifacts/Contam", FALSE, TRUE)
saveRDS(contamReadIDs, file = "results/AlignmentArtifacts/Contam/contamReadIDs.RDS")

maplen_mmrate_th <- "ge20_le0.1"
goodMapQualReadIDs <- sapply(sampleIDs, function(sampleID) {
    sampleName <- paste0("Sample_", sampleID)
    filename <- sprintf("data/E.chex/analyzed/%s/star/%s.star.primaryNoDup.%s.readID.txt.gz", sampleName, sampleName, maplen_mmrate_th)
    message(filename)
    readLines(filename)
}, simplify = FALSE)
dir.create("results/AlignmentArtifacts/MapQual", FALSE, TRUE)
saveRDS(goodMapQualReadIDs, file = sprintf("results/AlignmentArtifacts/MapQual/goodMapQualReadIDs_%s.RDS", maplen_mmrate_th))

GRsQualOutInPairListFiltered <- sapply(qualOutInPairs, function(qualOutInPair) {
    GRsList <- GRsQualOutInPairList[[qualOutInPair]]
    sapply(sampleIDs, function(sampleID) {
        message(maplen_mmrate_th, " ", qualOutInPair, " ", sampleID)
        GRs <- GRsList[[sampleID]]
        if (length(GRs) > 0) {
            readNames <- mcols(GRs)[, "name"]
            readIDs <- sapply(strsplit(readNames, "-"), "[", 1)
            contam <- contamReadIDs[[sampleID]]
            good <- goodMapQualReadIDs[[sampleID]]
            idx <- (!readIDs %in% contam) & (readIDs %in% good)
            return(GRs[idx])
        } else {
            return(GRs)
        }
    }, simplify = FALSE)
}, simplify = FALSE)

for (qualOutInPair in qualOutInPairs) {
    GRs <- GRsQualOutInPairListFiltered[[qualOutInPair]]
    filename <- sprintf("results/PrimingSites/GRs%sFiltered_%s.RDS", qualOutInPair, maplen_mmrate_th)
    saveRDS(GRs, file = filename)
}

qualOutInPairs_merge_map <- list(
    "ABread5End" = c("Aread5End", "Bread5End"),                   # for better specificity
    "ABreadCmate5End" = c("Aread5End", "Bread5End", "Cmate5End"), # for better sensitivity
)

for (x in names(qualOutInPairs_merge_map)) {
    GRs_each <- sapply(qualOutInPairs_merge_map[[x]], function(y) {
        filename <- sprintf("results/PrimingSites/GRs%sFiltered.RDS", y)
            message(x, " < ", y, " ", filename)
            readRDS(filename)
    }, simplify = FALSE)
    GRs_merged <- sapply(sampleIDs, function(x) {
        unname(unlist(as(lapply(GRs_each, function(GRs) GRs[[x]]), "GRangesList")))
    }, simplify = FALSE)
    filename <- sprintf("results/PrimingSites/GRs%sFiltered.RDS", x)
    message(x, " > ", filename)
    saveRDS(GRs_merged, file = filename)
}

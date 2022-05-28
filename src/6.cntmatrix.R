source("src/lib/Genome.R")

## Load genomic feature data
EnsFeatures <- readRDS("data/GenomicFeatures/EnsFeatures.RDS")
EnsFlanksByGene <- readRDS("data/GenomicFeatures/EnsFlanksByGene.RDS")

## Load sample grouping data
SampleInfo <- read.csv("data/SampleInfo.csv", as.is = TRUE, check.names = FALSE)
rownames(SampleInfo) <- sampleIDs <- SampleInfo[, "SampleID"]

Species <- c("human", "mouse", "rat")
EnsIDToSymbolMap <- sapply(Species, function(species) structure(mcols(EnsFeatures[["Gene"]][[species]])[["symbol"]], names = names(EnsFeatures[["Gene"]][[species]])), simplify = FALSE)
EnsIDToSymbolMap$none <- EnsIDToSymbolMap$mouse
## Calculate priming frequency in each gene
## 0. in GeneExt
## 1. in Gene
## 2. in TSS +/-5kb flanking region
## 3. in promoters
## 4. in 5'UTR
## 5. in exons
## 6. in introns
## 7. in 3'UTR
## 8. in Downstream
## 9. in CDS
featureTypesByGene <- c(
    "GeneExt", 
    "Gene", 
    "Flank5kByGene", 
    "PromoterByGene", 
    "FiveUTRByGene", 
    "ExonByGene", 
    "IntronByGene", 
    "ThreeUTRByGene", 
    "DownstreamByGene", 
    "CDSByGene"
)
EnsFeaturesByGene <- sapply(
    Species, function(species) 
        sapply(featureTypesByGene, function(featureTypeByGene) { 
            message(species, " ", featureTypeByGene); 
            EnsFeatures[[featureTypeByGene]][[species]] 
            }, simplify = FALSE
        ), 
    simplify = FALSE
)
EnsFeaturesByGene$none <- EnsFeaturesByGene$mouse
EnsFlanksByGene$none <- EnsFlanksByGene$mouse

## Add Species == "none" for water controls
Species <- c(Species, "none")

mapqTh <- "ge20_le0.1"
qualOutInPair <- "ABreadCmate5End"
###########################################################################
## Priming sites counting for each genomic feature type, except TSS flanking
## areas.
## Feature IDs are Ensembl. 
###########################################################################
filename <- sprintf("results/PrimingSites/GRs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
GRs <- readRDS(filename)
PFsEnsID <- sapply(Species, function(species) { 
    message(species)
    featuresByGene <- EnsFeaturesByGene[[species]]            
    SIDs <- subset(SampleInfo, Species == species)[, "SampleID"]
    Genome$getPrimingRatesPerFeature(featureRangesList = featuresByGene, featureTypesByGene = featureTypesByGene, GRsList = GRs, sampleIDs = SIDs)
}, simplify = FALSE) 
dirname <- sprintf("", .PWD)
dir.create(dirname, FALSE, TRUE)
filename <- sprintf("results/PrimingRateGene/PFs%sFilteredEnsID_%s.RDS", qualOutInPair, mapqTh)
saveRDS(PFsEnsID, file = filename)

###########################################################################
## Priming sites counting for TSS flanking areas.
## Feature IDs are Ensembl. 
###########################################################################
flankDistLabs <- paste("Flank", c("100", "200", "500", "1k", "2k", "3k", "4k", "5k"), sep = "")
FlankPFsEnsID <- sapply(Species, function(species) {
    message(species)
    FlanksByGene <- EnsFlanksByGene[[species]] 
    SIDs <- subset(SampleInfo, Species == species)[, "SampleID"]
    Genome$getPrimingRatesPerFeature(featureRangesList = FlanksByGene, featureTypesByGene = flankDistLabs, GRsList = GRs, sampleIDs = SIDs)
}, simplify = FALSE)
dir.create(dirname, FALSE, TRUE)
filename <- sprintf("results/PrimingRateGene/FlankPFs%sFilteredEnsID_%s.RDS", qualOutInPair, mapqTh)
saveRDS(FlankPFsEnsID, file = filename)

###########################################################################
## Priming sites counting for each genomic feature type, except TSS flanking
## areas.
## Feature IDs are converted to gene symbols. 
###########################################################################
ncores <- 10
filename <- sprintf("results/PrimingRateGene/PFs%sFilteredEnsID_%s.RDS", qualOutInPair, mapqTh)
PFsEnsID <- readRDS(filename)
PFs <- sapply(Species, function(species) {
    sapply(featureTypesByGene, function(featureTypeByGene) { 
        message(species, " ", featureTypeByGene)
        mat <- PFsEnsID[[species]][[featureTypeByGene]]
        mat <- Stats$summarizeMatrix(mat, by = EnsIDToSymbolMap[[species]][rownames(mat)], ncores = ncores)
        idx <- rownames(mat) == "" | is.na(rownames(mat))
        mat[!idx, ]
    }, simplify = FALSE)
}, simplify = FALSE)
filename <- sprintf("results/PrimingRateGene/PFs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
saveRDS(PFs, file = filename)

###########################################################################
## Priming sites counting for TSS flanking areas.
## Feature IDs are converted to gene symbols. 
###########################################################################
ncores <- 10
flankDistLabs <- paste("Flank", c("100", "200", "500", "1k", "2k", "3k", "4k", "5k"), sep = "")
filename <- sprintf("results/PrimingRateGene/FlankPFs%sFilteredEnsID_%s.RDS", qualOutInPair, mapqTh)
PFsEnsID <- readRDS(filename)
PFs <- sapply(Species, function(species) { 
    IDMap <- EnsIDToSymbolMap[[species]]
    sapply(flankDistLabs, function(flankDistLab) {
        message(species, " ", flankDistLab)
        X <- PFsEnsID[[species]][[flankDistLab]]
        IDs <- rownames(X)
        mat <- Stats$summarizeMatrix(X, by = IDMap[IDs], ncores = ncores)
        idx <- rownames(mat) == "" | is.na(rownames(mat))
        mat[!idx, ]
    }, simplify = FALSE)
}, simplify = FALSE)
filename <- sprintf("results/PrimingRateGene/FlankPFs%sFiltered_%s.RDS", qualOutInPair, mapqTh)
saveRDS(PFs, file = filename)

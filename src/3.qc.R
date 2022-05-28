source("src/lib/Repo.R")
source("src/lib/CHEX.R")

SampleInfo <- read.csv("data/SampleInfo.csv", as.is = TRUE, check.names = FALSE)
rownames(SampleInfo) <- sampleIDs <- SampleInfo[, "SampleID"]

## 1. NGS QC stats from the SCAP-T pipeline
## Biological samples only
SampleNGSstats <- t(mapply(function(exptID, sampleID) 
    getNGSstats(exptID = exptID, exptName = NULL, sampleID = sampleID, baseDir = "/lab/repo", 
                analyzedDir = "analyzed", blastDir = "blast", trimDir = "trim", alignedDir = "star", countsDir = "verse",
                features = c("exon", "exon-lev3", "anti-exon", "intron-lev1-lev2", "intron-lev3", "intron", "anti-intron", 
                             "mito", "anti-mito", "intergenic", "lines_sines")
        ), 
    exptID = SampleInfo[, "ExptID"], 
    sampleID = SampleInfo[, "SampleID"])
)

## Extract in-house NGS pipeline stats, plus a new field "BLAST.Target.Perc"
SampleNGSstats <- SampleNGSstats[, c("STAR.Uniq-mapped Reads", "STAR.Uniq-mapped Perc", "STAR.Multi-mapped Perc", "STAR.Perc Bases Mismatched", "STAR.Avg Inp Read Len", "STAR.Avg Uniq Map Len", "VERSE.exons Level 1,2: non-SpikeIn Reads Counted", "VERSE.exons Level 1,2: Num Genes", "VERSE.exons Level 1,2: Avg Read Per Gene",  "VERSE.exons Level 1,2: Max Reads Per Gene", "VERSE.introns Level 1,2: Reads Counted" , "VERSE.introns Level 1,2: Num Genes", "VERSE.introns Level 1,2: Avg Read Per Gene" , "VERSE.introns Level 1,2: Max Reads Per Gene", "VERSE.mito: Reads Counted", "VERSE.mito: Num Genes", "VERSE.mito: Avg Read Per Gene", "VERSE.mito: Max Reads Per Gene", "VERSE.Intergenic: Reads Counted" , "VERSE.Intergenic: Num Regions", "VERSE.Intergenic: Avg Read Per Region",  "VERSE.Intergenic: Max Reads Per Region","VERSE.LINE Reads", "VERSE.SINE Reads", "VERSE.Num LINEs", "VERSE.Num SINEs", "VERSE.Neither LINEs nor SINEs", "VERSE.Perc: Exons Level 1,2",  "VERSE.Perc: Introns Level 1,2" , "VERSE.Perc: Mito", "VERSE.Perc: Intergenic", "BLAST.Human", "BLAST.Mouse", "BLAST.Rat", "BLAST.Bacteria")]
rownames(SampleNGSstats) <- sampleIDs
SampleNGSstatsHuman <- SampleNGSstats[SampleInfo[, "Species"] == "human", ]
SampleNGSstatsHuman$BLAST.Target.Perc <- SampleNGSstatsHuman[["BLAST.Human"]] / 5000 * 100
SampleNGSstatsMouse <- SampleNGSstats[SampleInfo[, "Species"] == "mouse", ]
SampleNGSstatsMouse$BLAST.Target.Perc <- SampleNGSstatsMouse[["BLAST.Mouse"]] / 5000 * 100
SampleNGSstatsRat <- SampleNGSstats[SampleInfo[, "Species"] == "rat", ]
SampleNGSstatsRat$BLAST.Target.Perc <- SampleNGSstatsRat[["BLAST.Rat"]] / 5000 * 100
SampleNGSstatsNone <- SampleNGSstats[SampleInfo[, "Species"] == "none", ]
SampleNGSstatsNone$BLAST.Target.Perc <- NA # Water samples should have no "target" genome
SampleNGSstats <- rbind(SampleNGSstatsHuman, SampleNGSstatsMouse, SampleNGSstatsRat, SampleNGSstatsNone)[sampleIDs, ]
SampleNGSstatsDf <- data.frame(SampleID = rownames(SampleNGSstats), SampleNGSstats, stringsAsFactors = FALSE, check.names = FALSE)
write.csv(SampleNGSstatsDf, file = "results/NGSQCStats/stats.csv", row.names = FALSE)


## 2. Barcode/primer quality stats
ReadStats <- read.csv("results/ReadQualStats/overall.tsv", as.is = TRUE, check.names = FALSE, head = FALSE, sep = "\t")
colnames(ReadStats) <- ReadStats[1, ]
ReadStats <- ReadStats[-1, -15]
for (i in 2:14) { ReadStats[, i] <- as.numeric(ReadStats[, i]) }
rownames(ReadStats) <- ReadStats[, "SampleID"]
ReadStats <- ReadStats[sampleIDs, ]

## Load ReadQualCategoryStats for PE and SE (if any)
QualStatsPE <- CHEX$readQualStats("results/ReadQualStats/perclass_PE.tsv")
QualStatsSE <- CHEX$readQualStats("results/ReadQualStats/perclass_SE.tsv")

## Combine PE and SE stats
QualStats <- rbind(QualStatsPE, QualStatsSE)
write.csv(QualStats, file = "results/ReadQualStats/perclass.csv", row.names = FALSE)

## Reorganize the long table into wide format:
QualStatsDfNum <- CHEX$makeWideTableQualStats(QualStats, SampleInfo)

## Note, here all numbers are by reads (not read pairs), hence they should sum up to the "Total reads" of `samtools flagstats`
stopifnot(QualStatsDfNum[sampleIDs, "Total"] == ReadStats[sampleIDs, "Total"])
#" [1] TRUE

QualStatsDfPerc <- QualStatsDfNum / QualStatsDfNum$Total
QualStatsDfNum <- data.frame(SampleID = rownames(QualStatsDfNum), QualStatsDfNum, stringsAsFactors = FALSE)
QualStatsDfPerc <- data.frame(SampleID = rownames(QualStatsDfPerc), QualStatsDfPerc, stringsAsFactors = FALSE)
write.csv(QualStatsDfNum, file = "results/ReadQualStats/final_num.csv", row.names = FALSE)
write.csv(QualStatsDfPerc, file = "results/ReadQualStats/final_perc.csv", row.names = FALSE)
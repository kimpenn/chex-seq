## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
## v1.40
source src/lib/Repo.bash
source src/lib/NGS.bash

## Parse the CHEXTRIM produced tag and append the barcode/primer class-subclass label
## input: BAM, barcode/primer config file
## output: BAM file
function annotateReadQual {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    inFileBase="star.primaryNoDup"
    outFileBase="$inFileBase.PrimerAnnotated"
    PrimerLengthMapFile="src/lib/ChexPrimerTable.conf"
    minL2p=6
    minLpC=6
    debug="info"
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --inFileBase)
                inFileBase="$2"
                shift; shift
                ;;
            --outFileBase)
                outFileBase="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --endType)
                endType="$2"
                shift; shift
                ;;
            --PrimerLengthMapFile)
                PrimerLengthMapFile="$2"
                shift; shift
                ;;
            --primerIdx2p)
                primerIdx2p="$2"
                shift; shift
                ;;
            --primerIdxpC)
                primerIdxpC="$2"
                shift; shift
                ;;
            --minL2p)
                minL2p="$2"
                shift; shift
                ;;
            --minLpC)
                minLpC="$2"
                shift; shift
                ;;
            --debug)
                debug="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift;
                ;;
            *)
                echo "annotateReadQual [--baseDir Data] [--exptDir E.chex] --sampleID scCLTdegenNucXXX --endType PE [--analyzedDir analyzed] [--alignedDir star] [--inFileBase star.primaryNoDup] [--outFileBase star.primaryNoDup.PrimerAnnotated] [--PrimerLengthMapFile src/lib/PrimerLengthTable.conf] --primerIdx2p 505 --primerIdxpC 302 [--minL2p 6] [--minLpC 6] [--debug info] [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    if [[ $verbose == "true" ]]; then
        echo "AnnotateReadQual.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$inFileBase.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam --endType $endType --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --PrimerLengthMapFile $PrimerLengthMapFile --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam.log"
    fi
    src/lib/AnnotateReadQual.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$inFileBase.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam --endType $endType --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --PrimerLengthMapFile $PrimerLengthMapFile --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam.log
    if [[ $verbose == "true" ]]; then
        echo "samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam.bai.log 2>&1"
    fi
    samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam.bai.log 2>&1
}

## Read the barcode/primer class-subclass label and filter for reads of desired quality
## input: BAM file
## output: BAM files
function filterReadQual {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    inFileBase="star.primaryNoDup"
    outFileBase="$inFileBase"
    debug="info"
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --inFileBase)
                inFileBase="$2"
                shift; shift
                ;;
            --outFileBase)
                outFileBase="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --endType)
                endType="$2"
                shift; shift
                ;;
            --debug)
                debug="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift;
                ;;
            *)
                echo "filterReadQual [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] --sampleID scCLTdegenNucXXX --endType PE [--inFileBase star.primaryNoDup] [--outFileBase star.primaryNoDup] [--debug info] [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    for qual in A1 A2 B1 B2 C1 C2; do
        for outputInPair in "read" mate both either; do
            if [[ $verbose == "true" ]]; then
                echo "src/lib/FilterReadQual.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$inFileBase.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam --qualities $qual --endType $endType --outputInPair $outputInPair --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam.log"
            fi
            src/lib/FilterReadQual.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$inFileBase.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam --qualities $qual --endType $endType --outputInPair $outputInPair --debug $debug  2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam.log
            if [[ $verbose == "true" ]]; then
                echo "samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam.bai.log 2>&1"
            fi
            samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam.bai.log 2>&1
        done
    done
    for qual in D; do ## It does not make sense to output 'read' mate or both for D class
        for outputInPair in either; do
            if [[ $verbose == "true" ]]; then
                echo "src/lib/FilterReadQual.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$inFileBase.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam --qualities $qual --endType $endType --outputInPair $outputInPair --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam.log"
            fi
            src/lib/FilterReadQual.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$inFileBase.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam --qualities $qual --endType $endType --outputInPair $outputInPair --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam.log
            if [[ $verbose == "true" ]]; then
                echo "samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam.bai.log 2>&1"
            fi
            samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.$qual$outputInPair.bam.bai.log 2>&1
        done
    done
}

## Merged reads of different subclasses into a single BAM file
## input: BAM files
## output: BAM files
function mergeQualOutputInPairs {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup"
    qualOutputInPairs="A1either,A2either B1either,B2either C1either,C2either"
    qualOutputInPairsMerged="Aeither Beither Ceither"
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;
            --qualOutputInPairs)
                qualOutputInPairs="$2"
                shift; shift
                ;;
            --qualOutputInPairsMerged)
                qualOutputInPairsMerged="$2"
                shift; shift
                ;;
            --verbose)
                verbose=true
                shift;
                ;;
            *)
                echo "mergeQualOutputInPairs --sampleID scCLTdegenNucxxx [--baseDir baseDir] [--exptDir exptDir] [--analyzedDir analyzedDir] [--alignedDir alignedDir] [--fileBase star.primaryNoDup] [--qualOutputInPairs 'A1either,A2either B1either,B2either C1either,C2either'] [--qualOutputInPairsMerged 'Aeither Beither Ceither'] [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    qualOutputInPairs=($qualOutputInPairs)
    qualOutputInPairsMerged=($qualOutputInPairsMerged)
    n=${#qualOutputInPairs[*]}
    for (( i=0; i<=$(( $n -1 )); i++ )); do
        qualOutputInPair=${qualOutputInPairs[$i]}
        qualOutputInPairMerged=${qualOutputInPairsMerged[$i]}
        if [[ "$verbose" == "true" ]]; then
            echo "samtools merge -f $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qualOutputInPairMerged}.bam $(eval echo $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.{$qualOutputInPair}.bam) > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qualOutputInPairMerged}.bam.log 2>&1"
        fi
        samtools merge -f $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qualOutputInPairMerged}.bam $(eval echo $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.{$qualOutputInPair}.bam) > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qualOutputInPairMerged}.bam.log 2>&1

        if [[ "$verbose" == "true" ]]; then
            echo "samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qualOutputInPairMerged}.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qualOutputInPairMerged}.bam.bai.log 2>&1"
            samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qualOutputInPairMerged}.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qualOutputInPairMerged}.bam.bai.log 2>&1
        fi
    done
}

## Deduplicate the adjacent priming reads to represent a single priming event
## input: BAM file
## output: deduplicated BAM file, BED file with per-bin priming counts
function removeAdjReads {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup"
    addSoftClippedFlag=""
    endType="PE"
    ncores=1
    mindist=20
    priority="A1,A2,B1,B2,C1,C2"
    focusTag="A1:F,A2:F,B1:F,B2:F,C1:R,C2:R"
    verbose="false"
    debug="info"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;
            --endType)
                endType="$2"
                shift; shift
                ;;
            --addSoftClipped)
                addSoftClippedFlag="--addSoftClipped"
                shift;
                ;;
            --ncores)
                ncores="$2"
                shift; shift
                ;;
            --ncells)
                ncells="$2"
                shift; shift
                ;;
            --priority)
                priority="$2"
                shift; shift
                ;;
            --mindist)
                mindist="$2"
                shift; shift
                ;;
            --focusTag)
                focusTag="$2"
                shift; shift
                ;;
            --debug)
                debug="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift;
                ;;
            *)
                echo "removeAdjReads [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] [--fileBase star.primaryNoDup] [--endType PE] [--addSoftClipped] [--mindist 20] [--priority 'A1,A2,B1,B2,C1,C2'] [--focusTag 'A1:F,A2:F,B1:F,B2:F,C1:R,C2:R'] [--verbose] [--debug info]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    for qualOutputInPair in Aeither Beither Ceither; do
        if [ "$verbose" == "true" ]; then
            echo "src/lib/RemoveAdjReads.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.bam --outDupFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.AdjDup.bam --outBedFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.AdjCnts.bed --outStatsFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.tsv --endType $endType --mindist $mindist $addSoftClippedFlag --focusTag $focusTag --ncores $ncores --ncells $ncells --priority $priority --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.log"
        fi
        src/lib/RemoveAdjReads.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.bam --outDupFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.AdjDup.bam --outBedFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.AdjCnts.bed --outStatsFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.tsv --endType $endType --mindist $mindist $addSoftClippedFlag --focusTag $focusTag --ncores $ncores --ncells $ncells --priority $priority  --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.log
        if [ "$verbose" == "true" ]; then
            echo "samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.bam.bai.log 2>&1"
        fi
        samtools index $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPair.NoAdj.bam.bai.log 2>&1
    done
}

## Convert BAM to (1) priming sites, (2) primer-extension intervals, (3) both (1) and (2) preserving the strandedness
## input: BAM file
## output: BED file
function doBamToBed {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup"
    addSoftClippedFlag=""
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --endType)
                endType="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;
            --addSoftClipped)
                addSoftClippedFlag="--addSoftClipped"
                shift;
                ;;
            --verbose)
                verbose=true
                shift
                ;;
            *)
                echo "doBamToBed [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] --sampleID scCLTdegenNucXXX --endType PE [--fileBase star.primaryNoDup] [--addSoftClipped] [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    declare -A QualOutputInPairNoAdjMap=(["Aeither"]="Aeither.NoAdj" ["Beither"]="Beither.NoAdj" ["Ceither"]="Ceither.NoAdj" ["Deither"]="Deither")
    declare -A QualOutputInPairFocusTagMap=(["Aeither"]="F" ["Beither"]="F" ["Ceither"]="R" ["Deither"]="F")
    declare -A QualOutputInPairOutTypesMap=(["Aeither"]="5End ReadAndFrag Frag" ["Beither"]="5End ReadAndFrag Frag" ["Ceither"]="5End ReadAndFrag Frag" ["Deither"]="5End ReadAndFrag Frag")
    declare -A QualOutputInPairDoesOutNoFocusFlagMap=(["Aeither"]="" ["Beither"]="" ["Ceither"]="" ["Deither"]="--doesOutNoFocus")
    for qualOutputInPair in Aeither Beither Ceither Deither; do
        qualOutputInPairNoAdj=${QualOutputInPairNoAdjMap[$qualOutputInPair]}
        focusTag=${QualOutputInPairFocusTagMap[$qualOutputInPair]}
        outTypes=${QualOutputInPairOutTypesMap[$qualOutputInPair]}
        doesOutNoFocusFlag=${QualOutputInPairDoesOutNoFocusFlagMap[$qualOutputInPair]}
        for outType in $outTypes; do
            if [[ "$verbose" == "true" ]]; then
                echo "src/lib/BamToBed.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPairNoAdj.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPairNoAdj.$outType.bed $addSoftClippedFlag $doesOutNoFocusFlag --focusTag $focusTag --endType $endType --outType $outType --verbose 2>$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPairNoAdj.$outType.bed.log"
            fi
            src/lib/BamToBed.pl --inFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPairNoAdj.bam --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPairNoAdj.$outType.bed $addSoftClippedFlag $doesOutNoFocusFlag --focusTag $focusTag --endType $endType --outType $outType --verbose 2>$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qualOutputInPairNoAdj.$outType.bed.log
        done
    done
}

## Filter out unwanted genomic locations
## input: BED file
## output: BED file
function filterBed {
    baseDir="Data"
    exptDir="E.chex"
    alignedDir="star"
    blacklistFile=""
    local fileBase="star.primaryNoDup.Aeither.NoAdj.5End"
    outBase="NoBlacklisted"
    analyzedDir="analyzed"
    genomeFile="Data/Database/Genome/Human/UCSC/hg38.genome"
    extraArgs=""
    verbose="false"
    while [[ $# -gt 0 ]]; do
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;    
            --outBase)
                outBase="$2"
                shift; shift
                ;;
            --blacklistFile)
                blacklistFile="$2"
                shift; shift;
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            --extraArgs)
                shift
                extraArgs="$*"
                break
                ;;
            *)
                echo "filterBed [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] [--blacklistFile Data/Database/ENCODE/Blacklist/Kundaje/hg38.blacklist.bed] --sampleID scCLTdegenNucxxx [--fileBase star.primaryNoDup.Aeither.NoAdj.5End] [--outBase NoBlacklisted] [--verbose] [--extraArgs]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    outFileBase=$fileBase.$outBase
    ## TODO: to remove the complete query region instead of break it into multiple pieces, we should use -A
    if [[ "$verbose" == "true" ]]; then
        echo "bedtools subtract -A $extraArgs -a $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.bed -b $blacklistFile >$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bed 2>$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bed.log"
    fi
    bedtools subtract -A $extraArgs -a $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.bed -b $blacklistFile >$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bed 2>$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bed.log
}

## Filter out unwanted genomic locations by species
## input: BED file
## output: BED file
function filterBedBySpecies {
    baseDir="Data"
    exptDir="E.chex"
    alignedDir="star"
    local fileBase="star.primaryNoDup.Aeither.NoAdj.5End"
    outBase="NoBlacklisted"
    analyzedDir="analyzed"
    extraArgs=""
    verbose="false"
    verboseFlag=""
    while [[ $# -gt 0 ]]; do
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;    
            --outBase)
                outBase="$2"
                shift; shift
                ;;
            --species)
                species="$2"
                shift; shift;
                ;;
            --species)
                species="$2";
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            --extraArgs)
                shift
                extraArgs="$*"
                break
                ;;
            *)
                echo "filterBedBySpecies [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] --species human --sampleID scCLTdegenNucxxx [--fileBase star.primaryNoDup.Aeither.NoAdj.5End] [--outBase NoBlacklisted] [--verbose] [--extraArgs]"
                return 1
                ;;
        esac
    done
    declare -A BlacklistFiles=(["mouse"]="Data/ENCODE/Blacklist/Kundaje/mm10.blacklist.bed" ["human"]="Data/ENCODE/Blacklist/Kundaje/hg38.blacklist.bed" ["rat"]="Data/ENCODE/Blacklist/Kundaje/rn6.blacklist.bed" ["none"]="Data/ENCODE/Blacklist/Kundaje/mm10.blacklist.bed")
    blacklistFile=${BlacklistFiles[$species]}
    if [[ "$verbose" == "true" ]]; then
        verboseFlag="--verbose"
        echo "filterBed --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --blacklistFile $blacklistFile --sampleID $sampleID --fileBase $fileBase --outBase $outBase $verboseFlag --extraArgs $extraArgs"
    fi
    filterBed --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --blacklistFile $blacklistFile --sampleID $sampleID --fileBase $fileBase --outBase $outBase $verboseFlag --extraArgs $extraArgs
}

## Filter out unwanted genomic locations by species specific blacklists from the ENCODE Project
## input: BED file
## output: BED file
function removeBlacklisted {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup"
    outBase="NoBlacklisted"
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --species)
                species="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;
            --outBase)
                outBase="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            *)
                echo "removeBlacklisted --species human --sampleID scCLTdegenNucXXX [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] [--fileBase star.primaryNoDup] [--verbose]"
                return 1
                ;;
        esac
    done
    declare -A QualOutputInPairNoAdjMap=(["Aeither"]="Aeither.NoAdj" ["Beither"]="Beither.NoAdj" ["Ceither"]="Ceither.NoAdj" ["Deither"]="Deither")
    declare -A QualOutputInPairOutTypesMap=(["Aeither"]="5End ReadAndFrag Frag" ["Beither"]="5End ReadAndFrag Frag" ["Ceither"]="5End ReadAndFrag Frag" ["Deither"]="5End ReadAndFrag Frag")
    for qualOutputInPair in Aeither Beither Ceither Deither; do
        qualOutputInPairNoAdj=${QualOutputInPairNoAdjMap[$qualOutputInPair]}
        outTypes=${QualOutputInPairOutTypesMap[$qualOutputInPair]}
        for outType in $outTypes; do
            inFileBase=$fileBase.$qualOutputInPairNoAdj.$outType
            if [[ "$verbose" == "true" ]]; then
                echo "filterBedBySpecies --species $species --sampleID $sampleID --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $inFileBase --outBase $outBase --verbose"
            fi
            filterBedBySpecies --species $species --sampleID $sampleID --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $inFileBase --outBase $outBase --verbose
        done
    done
}

## Parse CHEXTRIM annotated read names and output the RMDUP copy number and barcode/primer trimmed length
## input: FASTQ files from the CHEXTRIM folder
## output: TSV file
function parseChexTrimAnnot {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    chextrimDir="chextrim"
    verbose="false"
    debug="info"
    PrimerLengthMapFile="src/lib/ChexPrimerTable.conf"
    minL2p=6
    minLpC=6
    inFileBase=unaligned
    inFileExt=fq.gz
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --endType)
                endType="$2"
                shift; shift
                ;;
            --chextrimDir)
                chextrimDir="$2"
                shift; shift
                ;;
            --PrimerLengthMapFile)
                PrimerLengthMapFile="$2"
                shift; shift
                ;;
            --primerIdx2p)
                primerIdx2p="$2"
                shift; shift
                ;;
            --primerIdxpC)
                primerIdxpC="$2"
                shift; shift
                ;;
            --minL2p)
                minL2p="$2"
                shift; shift
                ;;
            --minLpC)
                minLpC="$2"
                shift; shift
                ;;
            --inFileBase)
                inFileBase="$2"
                shift; shift
                ;;
            --inFileExt)
                inFileExt="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift;
                ;;
            --debug)
                debug="$2"
                shift; shift
                ;;
            *)
                echo "parseChexTrimAnnot [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--chextrimDir chextrim] --sampleID scCLTdegenNucXXX [--endType PE] [--PrimerLengthMapFile src/lib/ChexPrimerTable.conf] --primerIdx2p 505 --primerIdxpC 307 [--minL2p 6] [--minLpC 6] [--inFileBase unaligned] [--inFileExt fa.gz] [--verbose] [--debug info]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    if [[ $endType == "PE" ]]; then
        if [[ $version == "true" ]]; then
            echo "src/lib/ParseChexTrimAnnot.pl --inFile1 $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/${inFileBase}_1.$inFileExt --inFile2 $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/${inFileBase}_2.$inFileExt --endType $endType --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.annot.tsv.gz --PrimerLengthMapFile $PrimerLengthMapFile --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.annot.log"
        fi
        src/lib/ParseChexTrimAnnot.pl --inFile1 $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/${inFileBase}_1.$inFileExt --inFile2 $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/${inFileBase}_2.$inFileExt --endType $endType --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.annot.tsv.gz --PrimerLengthMapFile $PrimerLengthMapFile --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.annot.log
    elif [[ $endType == "SE" ]]; then
        if [[ $version == "true" ]]; then
            echo "src/lib/ParseChexTrimAnnot.pl --inFile1 $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/${inFileBase}_1.$inFileExt --endType $endType --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.annot.tsv.gz --PrimerLengthMapFile $PrimerLengthMapFile --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.annot.log"
        fi
        src/lib/ParseChexTrimAnnot.pl --inFile1 $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/${inFileBase}_1.$inFileExt --endType $endType --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.annot.tsv.gz --PrimerLengthMapFile $PrimerLengthMapFile --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --debug $debug 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.annot.log
    else 
        echo "--endType can only be PE or SE!"
        return 1
    fi
}

## Parse CHEXTRIM produced log files to have the distribution stats of barcode/primer locations
## input: LOG files from the CHEXTRIM folder
## output: TSV file
function parseChexTrimStats {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    chextrimDir="chextrim"
    verbose="false"
    debug="info"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --endType)
                endType="$2"
                shift; shift
                ;;
            --chextrimDir)
                chextrimDir="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift;
                ;;
            --debug)
                debug="$2"
                shift; shift
                ;;
            *)
                echo "parseChexTrimStats [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--chextrimDir chextrim] --sampleID scCLTdegenNucXXX --endType PE [--verbose] [--debug info]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    if [[ $verbose == "true" ]]; then
        echo "src/lib/ParseChexTrimStats.pl --filePrefix $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.chextrim --endType $endType --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.stats.tsv 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.stats.log"
    fi
    src/lib/ParseChexTrimStats.pl --filePrefix $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.chextrim --endType $endType --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.stats.tsv 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$sampleName.$chextrimDir.stats.log
}

## Plot bases of the random sampled reads and visualize the barcode/primer sequences
## input: primer sequences in FASTA, random sampled reads in FASTA
## output: PNG file
function plotFaPrimers {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    chextrimDir="chextrim"
    blastDir="blast"
    primerFilePrefix="chex-seq-app"
    inFile="raw.fa"
    outFileBase="baseviz"
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --chextrimDir)
                chextrimDir="$2"
                shift; shift
                ;;
            --blastDir)
                blastDir="$2"
                shift; shift
                ;;
            --primerFilePrefix)
                primerFilePrefix="$2"
                shift; shift
                ;;
            --inFile)
                inFile="$2"
                shift; shift
                ;;
            --outFileBase)
                outFileBase="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift;
                ;;
            *)
                echo "plotFaPrimers [--baseDir Data] [--exptDir E.chex] --sampleID scCLTdegenNucXXX [--analyzedDir analyzed] [--chextrimDir chextrim] [--blastDir blast] [--inFile raw.fa] [--primerFilePrefix chex-seq-app] [--outFileBase baseviz] [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    primerFile=$(echo $baseDir/$exptDir/$analyzedDir/$sampleName/$chextrimDir/$primerFilePrefix*.fa)
    if [[ $verbose == "true" ]]; then
        echo "Rscript src/lib/fa_baseviz.pl --faFile1 $baseDir/$exptDir/$analyzedDir/$sampleName/$blastDir/$inFile --primerFile $primerFile --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$blastDir/$sampleName.$outFileBase.png --main $sampleID --width 7 --height 7 --res 300 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$blastDir/$sampleName.$outFileBase.log"
    fi
    Rscript src/lib/fa_baseviz.pl --faFile1 $baseDir/$exptDir/$analyzedDir/$sampleName/$blastDir/$inFile --primerFile $primerFile --outFile $baseDir/$exptDir/$analyzedDir/$sampleName/$blastDir/$sampleName.$outFileBase.png --main $sampleID --width 7 --height 7 --res 300 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$blastDir/$sampleName.$outFileBase.log
}

## Get alignment stats per barcode/primer quality class-subclass
## intput: BAM files
## output: TXT files
function doSamFlagStats {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup"
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            *)
                echo "doSamFlagStats [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] [--fileBase star.primaryNoDup] --sampleID scCLTdegenNucXXX [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID

    if [[ $verbose == "true" ]]; then
        echo "samtools flagstat $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.stats.txt 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.stats.log"
    fi
    samtools flagstat $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.stats.txt 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.stats.log

    for qual in A1 A2 B1 B2 C1 C2; do
        for outputInPair in "read" mate both either; do
            if [[ $verbose == "true" ]]; then
                echo "samtools flagstat $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.stats.txt 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.stats.log"
            fi
            samtools flagstat $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.stats.txt 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.stats.log
        done
    done
    for qual in D; do
        for outputInPair in either; do
            if [[ $verbose == "true" ]]; then
                echo "samtools flagstat $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.stats.txt 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.stats.log"
            fi
            samtools flagstat $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.bam > $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.stats.txt 2> $baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.${fileBase}.${qual}${outputInPair}.stats.log
        done
    done
}

## Gather alignemnt stats over multiple samples
## input: BAM file
## output: TXT file
function summarizeBamFlagStats {
    baseDir="Data"
    exptName="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup"
    verbose=false
    while [[ $# -gt 0 ]]; do
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptName)
                exptName="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;    
            --outFile)
                outFile="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("$@")
                break
                ;;
            *)
                echo "summarizeBamFlagStats [--baseDir Data] [--exptName E.chex] [--analyzedDir analyzed] [--alignedDir star] [--fileBase star.primaryNoDup] --outFile stats.tsv [--verbose] --sampleIDs scCLTdegenNucxxx [scCLTdegenNucxxx...]"
                return 1
                ;;
        esac
    done
    if [[ $verbose == "true" ]]; then
        echo summarizeBamFlagStats --baseDir $baseDir --exptName $exptName --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase --outFile $outFile --sampleIDs ${sampleIDs[@]}
    fi

    echo -e "SampleID\tTotal\tSecondary\tSupplementary\tDuplicates\tMapped\tPaired\tRead1\tRead2\tProperly paired\tWith itself and mate mapped\tSingletons\tWith mate mapped to a different chr\tWith mate mapped to a different chr (mapQ>=5)" > $outFile
    for sampleID in ${sampleIDs[@]}; do
        sampleName=Sample_$sampleID
        filename=$baseDir/$exptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.stats.txt
        echo -n -e "$sampleID\t" >> $outFile
        cat $filename | cut -d "+" -f 1 | sed 's/ //' | tr '\n' '\t' >> $outFile
        echo "" >> $outFile
    done
}

## Gather alignemnt stats over multiple samples, per barcode/primer quality class-subclass
## input: BAM files
## output: TXT files
function summarizeBamFlagQualStats {
    baseDir="Data"
    exptName="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup"
    verbose="false"
    quals="A1,A2,B1,B2,C1,C2,D"
    outputInPairs="read,mate,both,either" 
    while [[ $# -gt 0 ]]; do
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptName)
                exptName="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;    
            --outFile)
                outFile="$2"
                shift; shift
                ;;
            --quals)
                quals="$2"
                shift; shift
                ;;
            --outputInPairs)
                outputInPairs="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("$@")
                break
                ;;
            *)
                echo "summarizeBamFlagQualStats [--baseDir Data] [--exptName E.chex] [--analyzedDir analyzed] [--alignedDir star] [--fileBase star.primaryNoDup] --outFile out.tsv [--quals "A1,A2,B1,B2,C1,C2,D"] [--outputInPairs "read,mate,both,either"] [--verbose] --sampleIDs scCLTdegenNucxxx [scCLTdegenNucxxx...]"
                return 1
                ;;
        esac
    done
    if [[ $verbose == "true" ]]; then
        echo summarizeBamFlagStats --baseDir $baseDir --exptName $exptName --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase --outFile $outFile --quals $quals --outputInPairs $outputInPairs --sampleIDs ${sampleIDs[@]}
    fi

    echo -e "SampleID\tQual\tOutputInPair\tTotal\tSecondary\tSupplementary\tDuplicates\tMapped\tPaired\tRead1\tRead2\tProperly paired\tWith itself and mate mapped\tSingletons\tWith mate mapped to a different chr\tWith mate mapped to a different chr (mapQ>=5)" > $outFile
    for qual in $(echo $quals | tr "," " "); do
        for outputInPair in $(echo $outputInPairs | tr "," " "); do
            for sampleID in ${sampleIDs[@]}; do
                sampleName=Sample_$sampleID
                filename=$baseDir/$exptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.$qual$outputInPair.stats.txt
                if [[ -e $filename ]]; then
                    echo -n -e "$sampleID\t$qual\t$outputInPair\t" >> $outFile
                    cat $filename | cut -d "+" -f 1 | sed 's/ //' | tr '\n' '\t' >> $outFile
                    echo "" >> $outFile
                fi
            done
        done
    done
}

## Create folder structures mimicking KimLab /lab/repo, symlink files for this pipeline to run
## input: the root folder: /lab/repo
## output: the intermidiate folder: $HOME/Datasets/lab/repo and 
##         the end user's folder: Data/E.chex (pipeline starts here)
function initNewPipe { # for both SE and PE
    mergeLanes=false
    endType="PE"
    verboseFlag=""
    baseDir="Data"
    pooledExptName="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    rawDir="raw"
    initDir="init"
    seqType="CHEXSeq"
    PrimerLengthMapFile="src/lib/ChexPrimerTable.conf"
    minL2p=6
    minLpC=6
    debug="info"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptName)
                exptName="$2"
                shift; shift
                ;;
            --pooledExptName)
                pooledExptName="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --rawDir)
                rawDir="$2"
                shift; shift
                ;;
            --initDir)
                initDir="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --seqType)
                seqType="$2"
                shift; shift
                ;;
            --mergeLanes)
                mergeLanes=true
                shift;
                ;;
            --endType)
                endType="$2"
                shift; shift
                ;;
            --PrimerLengthMapFile)
                PrimerLengthMapFile="$2"
                shift; shift
                ;;
            --primerIdx2p)
                primerIdx2p="$2"
                shift; shift
                ;;
            --primerIdxpC)
                primerIdxpC="$2"
                shift; shift
                ;;
            --minL2p)
                minL2p="$2"
                shift; shift
                ;;
            --minLpC)
                minLpC="$2"
                shift; shift
                ;;
            --debug)
                debug="$2"
                shift; shift
                ;;
            --verbose)
                verboseFlag='--verbose'
                shift;
                ;;
            *)
                echo "initNewPipe --exptName E.xxx --sampleID scCLTdegenNucXXX [--baseDir Data] [--pooledExptName E.chex] [--analyzedDir analyzed] [--alignedDir star] [--seqType CHEXSeq] [--rawDir raw] [--initDir init] [--endType PE] [--mergeLanes] [--PrimerLengthMapFile src/lib/PrimerLengthTable.conf] --primerIdx2p 505 --primerIdxpC 302 [--minL2p 6] [--minLpC 6] [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    if [[ $verbose == "true" ]]; then
        echo "initRepo --inBaseDir /lab/repo --outBaseDir Data/Datasets/lab/repo --analyzedDir $analyzedDir $verboseFlag --seqType $seqType --exptName $exptName --sampleIDs $sampleID"
    fi
    initRepo --inBaseDir /lab/repo --outBaseDir Data/Datasets/lab/repo --analyzedDir $analyzedDir $verboseFlag --seqType $seqType --exptName $exptName --sampleIDs $sampleID
    if [[ "$mergeLanes" == "true" ]]; then
        if [[ $verbose == "true" ]]; then
            echo "mergeRawLanes --baseDir Data/Datasets/lab/repo --exptName $exptName --endType $endType $verboseFlag --sampleIDs $sampleID"
        fi
        mergeRawLanes --baseDir Data/Datasets/lab/repo --exptName $exptName --endType $endType $verboseFlag --sampleIDs $sampleID
    fi
    if [[ $verbose == "true" ]]; then
        echo "poolRawDir --baseDir Data/Datasets/lab/repo --pooledExptName $pooledExptName --rawDir $rawDir $verboseFlag --exptName $exptName --sampleIDs $sampleID"
    fi
    poolRawDir --baseDir Data/Datasets/lab/repo --pooledExptName $pooledExptName --rawDir $rawDir $verboseFlag --exptName $exptName --sampleIDs $sampleID
    if [[ $verbose == "true" ]]; then
        echo "poolAnalyzedDir --baseDir Data/Datasets/lab/repo --pooledExptName $pooledExptName --analyzedDir $analyzedDir $verboseFlag --exptName $exptName --seqType $seqType --sampleIDs $sampleID"
    fi
    poolAnalyzedDir --baseDir Data/Datasets/lab/repo --pooledExptName $pooledExptName --analyzedDir $analyzedDir $verboseFlag --exptName $exptName --seqType $seqType --sampleIDs $sampleID
    if [[ $verbose == "true" ]]; then
        echo "symlinkRawToAnalyzedInit --baseDir Data/Datasets/lab/repo --exptName $pooledExptName --rawDir $rawDir --analyzedDir $analyzedDir --initDir $initDir $verboseFlag --sampleIDs $sampleID"
    fi   
    symlinkRawToAnalyzedInit --baseDir Data/Datasets/lab/repo --exptName $pooledExptName --rawDir $rawDir --analyzedDir $analyzedDir --initDir $initDir $verboseFlag --sampleIDs $sampleID
    if [[ $verbose == "true" ]]; then
        echo "annotateReadQual --baseDir $baseDir --exptDir $pooledExptName --sampleID $sampleID --endType $endType --PrimerLengthMapFile $PrimerLengthMapFile --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --inFileBase star.posSorted --outFileBase star.PrimerAnnotated --debug $debug $verboseFlag"
    fi    
    annotateReadQual --baseDir $baseDir --exptDir $pooledExptName --sampleID $sampleID --endType $endType --PrimerLengthMapFile $PrimerLengthMapFile --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --inFileBase star.posSorted --outFileBase star.PrimerAnnotated --debug $debug $verboseFlag
    if [[ $verbose == "true" ]]; then
        echo "markDup --baseDir $baseDir --exptName $pooledExptName --analyzedDir $analyzedDir --alignedDir $alignedDir --inFileBase star.PrimerAnnotated --outFileBase star.DupMarked $verboseFlag --index --sampleIDs $sampleID"
    fi
    markDup --baseDir $baseDir --exptName $pooledExptName --analyzedDir $analyzedDir --alignedDir $alignedDir --inFileBase star.PrimerAnnotated --outFileBase star.DupMarked $verboseFlag --index --sampleIDs $sampleID
    if [[ $verbose == "true" ]]; then
        echo "rmDup --baseDir $baseDir --exptName $pooledExptName --analyzedDir $analyzedDir --alignedDir $alignedDir --inFileBase star.DupMarked --outFileBase star.NoDup $verboseFlag --index --sampleIDs $sampleID"
    fi
    rmDup --baseDir $baseDir --exptName $pooledExptName --analyzedDir $analyzedDir --alignedDir $alignedDir --inFileBase star.DupMarked --outFileBase star.NoDup $verboseFlag --index --sampleIDs $sampleID
    if [[ $verbose == "true" ]]; then
        echo "selectReads --baseDir $baseDir --exptName $pooledExptName --analyzedDir $analyzedDir --alignedDir $alignedDir --inFileBase star.NoDup --outFileBase star.primaryNoDup --flags '-F 0x100' --verbose --index --sampleIDs $sampleID"
    fi
    selectReads --baseDir $baseDir --exptName $pooledExptName --analyzedDir $analyzedDir --alignedDir $alignedDir --inFileBase star.NoDup --outFileBase star.primaryNoDup --flags '-F 0x100' --verbose --index --sampleIDs $sampleID
    ## in star, strip the read annotation
    if [[ $verbose == "true" ]]; then
        echo "samtools view -F 0x100 -h $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.posSorted.bam | perl -ane 'print && next if /^@/; \$F[0] = (split(q(-), \$F[0]))[0]; print join(qq(\t), @F), qq(\n)' | samtools view -Sb - > $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.bam 2> $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.bam.log"
    fi
    samtools view -F 0x100 -h $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.posSorted.bam | perl -ane 'print && next if /^@/; $F[0] = (split(q(-), $F[0]))[0]; print join(qq(\t), @F), qq(\n)' | samtools view -Sb - > $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.bam 2> $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.bam.log
    if [[ $verbose == "true" ]]; then  
        echo "samtools view $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.bam | perl -ane 'BEGIN {%R} \$F[0]=(split(\"-\", \$F[0]))[0]; print \$F[0], \"\n\" if ++\$R{\$F[0]}==1' | gzip -c > $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.readID.txt.gz 2>$baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.readID.log"
    fi
    samtools view $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.bam | perl -ane 'BEGIN {%R} $F[0]=(split("-", $F[0]))[0]; print $F[0], "\n" if ++$R{$F[0]}==1' | gzip -c > $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.readID.txt.gz 2>$baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primary.readID.log
    if [[ $verbose == "true" ]]; then
        echo "samtools view $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primaryNoDup.bam | perl -ane 'BEGIN {%R} \$F[0]=(split(\"-\", \$F[0]))[0]; print \$F[0], \"\n\" if ++\$R{\$F[0]}==1' | gzip -c > $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primaryNoDup.readID.txt.gz 2>$baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primaryNoDup.readID.log"
    fi
    samtools view $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primaryNoDup.bam | perl -ane 'BEGIN {%R} $F[0]=(split("-", $F[0]))[0]; print $F[0], "\n" if ++$R{$F[0]}==1' | gzip -c > $baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primaryNoDup.readID.txt.gz 2>$baseDir/$pooledExptName/$analyzedDir/$sampleName/$alignedDir/$sampleName.star.primaryNoDup.readID.log
}

## Align reads in the TRIM folder to alternative species and detect contaminant reads
## input: FASTQ files
## output: BAM files, TSV files with mapping scores, best hits and most likely source of species
function alignAltSpecies {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    trimDir="trim"
    alignedDir="star"
    altAlignedDir="starAlt"
    inFilePrefix="unaligned"
    inFileExt="fq"
    endType="PE"
    readLength=75
    genomeBaseRepo="/lab/repo/resources/lib/star/2.4.0h1"
    genomeBaseCustom="Data/Database/Indices/STAR"
    verbose="false"
    ncores=5
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --trimDir)
                trimDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --altAlignedDir)
                altAlignedDir="$2"
                shift; shift
                ;;
            --genomeBaseRepo)
                genomeBaseRepo=$2
                shift; shift
                ;;
            --genomeBaseCustom)
                genomeBaseCustom=$2
                shift; shift
                ;;
            --inFilePrefix)
                inFilePrefix=$2
                shift; shift
                ;;
            --inFileExt)
                inFileExt=$2
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --species)
                species=$2
                shift; shift
                ;;
            --endType)
                endType=$2
                shift; shift
                ;;
            --readLength)
                readLength=$2
                shift; shift
                ;;
            --ncores)
                ncores=$2
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            *) 
                echo "alignAltSpecies [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--trimDir trim] [--alignedDir star] [--altAlignedDir starAlt] [--genomeBaseRepo /lab/repo/resources/lib/star/2.4.0h1] [--genomeBaseCustom Data/Database/Indices/STAR] [--inFilePrefix unaligned] [--inFileExt fq] [--endType PE] [--readLength 75] [--ncores 5] [--verbose] --sampleID scCLTdegenNucXXX"
                return 1
                ;;
        esac
    done
    declare -A SpeciesAltMap=(["human"]="mouse bacteria mycoplasma" ["mouse"]="human bacteria mycoplasma" ["rat"]="humam bacteria mycoplasma" ["none"]="human bacteria mycoplasma")
    declare -A GenomePathMap=(["human"]="$genomeBaseRepo/hg38.gencode25.$readLength" ["mouse"]="$genomeBaseRepo/mm10.$readLength" ["rat"]="$genomeBaseRepo/rn6.$readLength" ["bacteria"]="$genomeBaseCustom/bacteria/top20" ["mycoplasma"]="$genomeBaseCustom/mycoplasma/mycoplasma")
    sampleName=Sample_$sampleID
    AltSpecies=${SpeciesAltMap[$species]}
    inDir=$baseDir/$exptDir/$analyzedDir/$sampleName/$trimDir
    outDir=$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir
    altOut=$baseDir/$exptDir/$analyzedDir/$sampleName/$altAlignedDir
    for altSpecies in $AltSpecies; do
        genomePath=${GenomePathMap[$altSpecies]}
        if [[ $genomePath == "" ]]; then
            echo "The species specified does not have a predefined STAR index! Exiting..."
            return 2
        fi
        altOutDir=$altOut/$altSpecies
        runlogPath="$altOutDir/Aligned.out.bam.log"
        logPath="$altOutDir/Log.final.out"
        statFile="$sampleName.star.stats.txt"
        statPath=$altOutDir/$statFile
        mkdir -p $altOutDir

        inFile1=$inDir/${inFilePrefix}_1.$inFileExt
        if [[ $endType == "PE" ]]; then
            inFile2=$inDir/${inFilePrefix}_2.$inFileExt
        fi
        if [[ ! -e $inFile1 ]]; then
            gzip -d ${inFile1}.gz -c > $inFile1
        fi
        if [[ $endType == "PE" ]]; then
            if [[ ! -e $inFile2 ]]; then
                gzip -d ${inFile2}.gz -c > $inFile2
            fi
        fi

        if [[ $verbose == "true" ]]; then
            echo "STAR --outFilterScoreMin 0 --outFilterScoreMinOverLread 0 --outFilterMismatchNmax 100 --outFilterMismatchNoverLmax 0.1 --outReadsUnmapped None --outSAMattributes All --outSAMunmapped 'Within KeepPairs' --genomeLoad LoadAndRemove --outSAMtype BAM Unsorted --outFilterMatchNmin 10 --outFilterMatchNminOverLread 0.4 --readFilesIn $inFile1 $inFile2 --alignIntronMax 1 --alignSJDBoverhangMin 999 --alignMatesGapMax 800 --genomeDir $genomePath --runThreadN $ncores --outFileNamePrefix $altOutDir/ > $runlogPath 2>&1"
        fi
        STAR --outFilterScoreMin 0 --outFilterScoreMinOverLread 0 --outFilterMismatchNmax 100 --outFilterMismatchNoverLmax 0.1 --outReadsUnmapped None --outSAMattributes All --outSAMunmapped 'Within KeepPairs' --genomeLoad LoadAndRemove --outSAMtype BAM Unsorted --outFilterMatchNmin 10 --outFilterMatchNminOverLread 0.4 --readFilesIn $inFile1 $inFile2 --alignIntronMax 1 --alignSJDBoverhangMin 999 --alignMatesGapMax 800 --genomeDir $genomePath --runThreadN $ncores --outFileNamePrefix $altOutDir/ > $runlogPath 2>&1

        if [[ $verbose == "true" ]]; then
            echo "mv $logPath $statPath"
        fi
        mv $logPath $statPath
    done

    ## Clean the read IDs in alternative alignments
    for altSpecies in $AltSpecies; do
        altOutDir=$altOut/$altSpecies
        if [[ $verbose == "true" ]]; then
            echo "samtools view -F 0x04 -h $altOutDir/Aligned.out.bam | perl -lane 'print && next if /^@/; \$F[0] = (split(\"-\", \$F[0]))[0]; print join(\"\t\", @F)' | samtools view -Sb - > $altOutDir/$sampleName.star.aligned.bam 2> $altOutDir/$sampleName.star.aligned.bam.log"
        fi
        samtools view -F 0x04 -h $altOutDir/Aligned.out.bam | perl -lane 'print && next if /^@/; $F[0] = (split("-", $F[0]))[0]; print join("\t", @F)' | samtools view -Sb - > $altOutDir/$sampleName.star.aligned.bam 2> $altOutDir/$sampleName.star.aligned.bam.log
        if [[ $verbose == "true" ]]; then
            echo "samtools sort -@ $ncores -o $altOutDir/$sampleName.star.aligned.posSorted.bam $altOutDir/$sampleName.star.aligned.bam 2> $altOutDir/$sampleName.star.aligned.posSorted.bam.log"
        fi
        samtools sort -@ $ncores -o $altOutDir/$sampleName.star.aligned.posSorted.bam $altOutDir/$sampleName.star.aligned.bam 2> $altOutDir/$sampleName.star.aligned.posSorted.bam.log
        if [[ $verbose == "true" ]]; then
            echo "samtools view -F 0x100 -h $altOutDir/$sampleName.star.aligned.posSorted.bam | samtools view -Sb > $altOutDir/$sampleName.star.aligned.primary.bam 2> $altOutDir/$sampleName.star.aligned.primary.bam.log"
        fi
        samtools view -F 0x100 -h $altOutDir/$sampleName.star.aligned.posSorted.bam | samtools view -Sb > $altOutDir/$sampleName.star.aligned.primary.bam 2> $altOutDir/$sampleName.star.aligned.primary.bam.log
    done

    ## Compare target and alternative alignments
    labels=$species
    infiles=$outDir/$sampleName.star.primary.bam
    for altSpecies in $AltSpecies; do
        altOutDir=$altOut/$altSpecies
        labels="$labels,$altSpecies"
        infiles="$infiles,$altOutDir/$sampleName.star.aligned.primary.bam"
    done
    if [[ $verbose == "true" ]]; then
        echo "src/lib/sam_best_hits.pl --inFiles $infiles --labels $labels -o $altOut/$sampleName.star.aligned.primary.best.tsv.gz 2> $altOut/$sampleName.star.aligned.primary.best.tsv.log"
    fi
    src/lib/sam_best_hits.pl --ncores $ncores --inFiles $infiles --labels $labels -o $altOut/$sampleName.star.aligned.primary.best.tsv.gz 2> $altOut/$sampleName.star.aligned.primary.best.tsv.log

    ## For each species (target or alternative), separate corresponding reads that best align to it
    for s in $species $AltSpecies; do
        if [[ $verbose == "true" ]]; then
            echo "zcat $altOut/$sampleName.star.aligned.primary.best.tsv.gz | perl -sne 'print && next if \$. == 1; chomp; @F=split(\"\t\", \$_); print \$_, \"\n\" if \$F[6] eq \$s' -- -s=$s | gzip -c > $altOut/$sampleName.star.aligned.primary.$s.tsv.gz 2> $altOut/$sampleName.star.aligned.primary.$s.tsv.log"
        fi
        zcat $altOut/$sampleName.star.aligned.primary.best.tsv.gz | perl -sne 'print && next if $. == 1; chomp; @F=split("\t", $_); print $_, "\n" if $F[6] eq $s' -- -s=$s | gzip -c > $altOut/$sampleName.star.aligned.primary.$s.tsv.gz 2> $altOut/$sampleName.star.aligned.primary.$s.tsv.log
    done
    
    ## Separate ties and define them as ambiguous
    if [[ $verbose == "true" ]]; then
        echo "zcat $altOut/$sampleName.star.aligned.primary.best.tsv.gz | perl -ne 'print && next if \$. == 1; chomp; @F=split(\"\t\", \$_); print \$_, \"\n\" if \$F[6] =~ /,/' | gzip -c > $altOut/$sampleName.star.aligned.primary.ambiguous.tsv.gz 2> $altOut/$sampleName.star.aligned.primary.ambiguous.tsv.log"
    fi
    zcat $altOut/$sampleName.star.aligned.primary.best.tsv.gz | perl -ne 'print && next if $. == 1; chomp; @F=split("\t", $_); print $_, "\n" if $F[6] =~ /,/' | gzip -c > $altOut/$sampleName.star.aligned.primary.ambiguous.tsv.gz 2> $altOut/$sampleName.star.aligned.primary.ambiguous.tsv.log

    ## Count read number for species best and ambiguous mapping
    for s in $species $AltSpecies ambiguous; do
        if [[ $verbose == "true" ]]; then
            echo "zcat $altOut/$sampleName.star.aligned.primary.$s.tsv.gz | sed 1d | wc -l > $altOut/$sampleName.star.aligned.primary.$s.nReadIDs.txt 2> $altOut/$sampleName.star.aligned.primary.$s.nReadIDs.txt.log"
        fi
        zcat $altOut/$sampleName.star.aligned.primary.$s.tsv.gz | sed 1d | wc -l > $altOut/$sampleName.star.aligned.primary.$s.nReadIDs.txt 2> $altOut/$sampleName.star.aligned.primary.$s.nReadIDs.txt.log
    done

    ## Separate reads that align better to alternative species and define them as contaminants
    if [[ $verbose == "true" ]]; then
        echo "zcat $altOut/$sampleName.star.aligned.primary.best.tsv.gz | perl -ne 'print && next if \$. == 1; chomp; @F=split(\"\t\", \$_); print \$_, \"\n\" unless \$F[6] =~ /human/' | gzip -c > $altOut/$sampleName.star.aligned.primary.contam.tsv.gz 2> $altOut/$sampleName.star.aligned.primary.contam.log"
    fi
    zcat $altOut/$sampleName.star.aligned.primary.best.tsv.gz | perl -ne 'print && next if $. == 1; chomp; @F=split("\t", $_); print $_, "\n" unless $F[6] =~ /human/' | gzip -c > $altOut/$sampleName.star.aligned.primary.contam.tsv.gz 2> $altOut/$sampleName.star.aligned.primary.contam.log
    if [[ $verbose == "true" ]]; then
        echo "zcat $altOut/$sampleName.star.aligned.primary.contam.tsv.gz | sed 1d | cut -f1 | gzip -c > $altOut/$sampleName.star.aligned.primary.contam.readID.txt.gz 2> $altOut/$sampleName.star.aligned.primary.contam.readID.log"
    fi
    zcat $altOut/$sampleName.star.aligned.primary.contam.tsv.gz | sed 1d | cut -f1 | gzip -c > $altOut/$sampleName.star.aligned.primary.contam.readID.txt.gz 2> $altOut/$sampleName.star.aligned.primary.contam.readID.log

    ## Compare target and alternative alignments (primaryNoDup)
    if [[ $verbose == "true" ]]; then
        echo "zcat $altOut/$sampleName.star.aligned.primary.best.tsv.gz | perl -MIO::Zlib -sane 'BEGIN {%R; @readIDs; \$FH=IO::Zlib->new(\"\$outDir/\$sampleName.star.primaryNoDup.readID.txt.gz\", \"r\"); while(<\$FH>) { chomp; push @readIDs, \$_; }} if (\$. == 1) {print} else {\$R{\$F[0]} = \$_} END { for \$r (@readIDs) { print \$R{\$r} } }'  -- -outDir=$outDir -sampleName=$sampleName | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.best.tsv.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.best.log"
    fi
    zcat $altOut/$sampleName.star.aligned.primary.best.tsv.gz | perl -MIO::Zlib -sane 'BEGIN {%R; @readIDs; $FH=IO::Zlib->new("$outDir/$sampleName.star.primaryNoDup.readID.txt.gz", "r"); while(<$FH>) { chomp; push @readIDs, $_; }} if ($. == 1) {print} else {$R{$F[0]} = $_} END { for $r (@readIDs) { print $R{$r} } }' -- -outDir=$outDir -sampleName=$sampleName | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.best.tsv.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.best.log

    ## For each species (target or alternative), separate corresponding reads that best align to it (primaryNoDup)
    for s in $species $AltSpecies; do
        if [[ $verbose == "true" ]]; then
            echo "zcat $altOut/$sampleName.star.aligned.primaryNoDup.best.tsv.gz | perl -sne 'print && next if \$. == 1; chomp; @F=split(\"\t\", \$_); print \$_, \"\n\" if \$F[6] eq \$s' -- -s=$s | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.$s.tsv.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.$s.tsv.log"
        fi
        zcat $altOut/$sampleName.star.aligned.primaryNoDup.best.tsv.gz | perl -sne 'print && next if $. == 1; chomp; @F=split("\t", $_); print $_, "\n" if $F[6] eq $s' -- -s=$s | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.$s.tsv.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.$s.tsv.log
    done

    ## Separate ties and define them as ambiguous (primaryNoDup)
    if [[ $verbose == "true" ]]; then
        echo "zcat $altOut/$sampleName.star.aligned.primaryNoDup.best.tsv.gz | perl -ne 'print && next if \$. == 1; chomp; @F=split(\"\t\", \$_); print \$_, \"\n\" if \$F[6] =~ /,/' | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.ambiguous.tsv.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.ambiguous.tsv.log"
    fi
    zcat $altOut/$sampleName.star.aligned.primaryNoDup.best.tsv.gz | perl -ne 'print && next if $. == 1; chomp; @F=split("\t", $_); print $_, "\n" if $F[6] =~ /,/' | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.ambiguous.tsv.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.ambiguous.tsv.log

    ## Separate reads that align better to alternative species and define them as contaminants (primaryNoDup)
    if [[ $verbose == "true" ]]; then
        echo "zcat $altOut/$sampleName.star.aligned.primaryNoDup.best.tsv.gz | perl -ne 'print && next if \$. == 1; chomp; @F=split(\"\t\", \$_); print \$_, \"\n\" unless \$F[6] =~ /human/' | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.contam.tsv.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.contam.log"
    fi
    zcat $altOut/$sampleName.star.aligned.primaryNoDup.best.tsv.gz | perl -ne 'print && next if $. == 1; chomp; @F=split("\t", $_); print $_, "\n" unless $F[6] =~ /human/' | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.contam.tsv.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.contam.log
    if [[ $verbose == "true" ]]; then
        echo "zcat $altOut/$sampleName.star.aligned.primaryNoDup.contam.tsv.gz | sed 1d | cut -f1 | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.contam.readID.txt.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.contam.readID.log"
    fi
    zcat $altOut/$sampleName.star.aligned.primaryNoDup.contam.tsv.gz | sed 1d | cut -f1 | gzip -c > $altOut/$sampleName.star.aligned.primaryNoDup.contam.readID.txt.gz 2> $altOut/$sampleName.star.aligned.primaryNoDup.contam.readID.log

    ## Count read number for species best and ambiguous mapping (primaryNoDup)
    for s in $species $AltSpecies ambiguous; do
        if [[ $verbose == "true" ]]; then
            echo "zcat $altOut/$sampleName.star.aligned.primaryNoDup.$s.tsv.gz | sed 1d | wc -l > $altOut/$sampleName.star.aligned.primaryNoDup.$s.nReadIDs.txt 2> $altOut/$sampleName.star.aligned.primaryNoDup.$s.nReadIDs.txt.log"
        fi
        zcat $altOut/$sampleName.star.aligned.primaryNoDup.$s.tsv.gz | sed 1d | wc -l > $altOut/$sampleName.star.aligned.primaryNoDup.$s.nReadIDs.txt 2> $altOut/$sampleName.star.aligned.primaryNoDup.$s.nReadIDs.txt.log
    done
}

## Get alignment stats regarding mapped length and mismatch bases
## intput: BAM file
## output: TSV file
function getMappedLengthMismatches {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup"
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            *)
                echo "getMappedLengthMismatches [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] [--fileBase star.primaryNoDup] --sampleID scCLTdegenNucXXX [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    inFile=$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.bam
    outFile=$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.mappedLength_mismatches.tsv.gz
    logFile=$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.mappedLength_mismatches.log

    if [[ $verbose == "true" ]]; then
        echo "sam_mlen_mm.pl -i $inFile 2> $logFile | perl -ane '\$F[0] = (split(qq(-), \$F[0]))[0]; print join(qq(\t), @F), qq(\n)' | gzip -c > $outFile"
    fi
    sam_mlen_mm.pl -i $inFile 2> $logFile | perl -ane '$F[0] = (split(qq(-), $F[0]))[0]; print join(qq(\t), @F), qq(\n)' | gzip -c > $outFile
}

## Extact read IDs with good mapping quality
## intput: TSV file
## output: TXT file
function filterMappedLengthMismatches {
    baseDir="Data"
    exptDir="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    fileBase="star.primaryNoDup.mappedLength_mismatches"
    mlen=30
    mm=0.1
    verbose="false"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;
            --mlen)
                mlen="$2"
                shift; shift
                ;;
            --mm)
                mm="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            *)
                echo "filterMappedLengthMismatches [--baseDir Data] [--exptDir E.chex] [--analyzedDir analyzed] [--alignedDir star] [--fileBase star.primaryNoDup.mappedLength_mismatches] [--mlen 30] [--mm 0.1] --sampleID scCLTdegenNucXXX [--verbose]"
                return 1
                ;;
        esac
    done
    sampleName=Sample_$sampleID
    inFile=$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.tsv.gz
    outFile=$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.ge${mlen}_le${mm}.readID.txt.gz
    logFile=$baseDir/$exptDir/$analyzedDir/$sampleName/$alignedDir/$sampleName.$fileBase.ge${mlen}_le${mm}.readID.log

    if [[ $verbose == "true" ]]; then
        echo "gzip -cd $inFile | sed 1d | perl -sane 'print \$F[0], qq(\n) if \$F[2] >= \$mlen && \$F[3]/\$F[1] <= \$mm' -- -mlen=$mlen -mm=$mm 2> $logFile | gzip -c > $outFile"
    fi
    gzip -cd $inFile | sed 1d | perl -sane 'print $F[0], qq(\n) if $F[2] >= $mlen && $F[3]/$F[1] <= $mm' -- -mlen=$mlen -mm=$mm 2> $logFile | gzip -c > $outFile
}

## Run the CHEX-seq pipeline
## input: unprocessed files in the folder Data/E.chex
## output: processed files in the folder Data/E.chex
function doNewPipe {
    baseDir="Data"
    exptDir="E.chex"
    rawDir="raw"
    initDir="init"
    analyzedDir="analyzed"
    chextrimDir="chextrim"
    trimDir="trim"
    alignedDir="star"
    altAlignedDir="starAlt"
    fileBase="star.primaryNoDup"
    PrimerLengthMapFile="src/lib/ChexPrimerTable.conf"
    primerIdx2p="505"
    primerIdxpC="302"
    species="human"
    endType="PE"
    readLength=75
    minL2p=6
    minLpC=6
    mindist=20
    ncells=1
    priority="A1,A2,B1,B2,C1,C2"
    focusTag="A1:F,A2:F,B1:F,B2:F,C1:R,C2:R"
    mergeLanesFlag=""
    addSoftClippedFlag=""
    genomeBaseRepo="/lab/repo/resources/lib/star/2.4.0h1"
    genomeBaseCustom="Data/Database/Indices/STAR"
    ncores=5
    verbose="false"
    verboseFlag=""
    debug="info"
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --exptDir)
                exptDir="$2"
                shift; shift
                ;;
            --rawDir)
                rawDir="$2"
                shift; shift
                ;;
            --initDir)
                initDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --chextrimDir)
                chextrimDir="$2"
                shift; shift
                ;;
            --trimDir)
                trimDir="$2"
                shift; shift
                ;;
            --alignedDir)
                alignedDir="$2"
                shift; shift
                ;;
            --altAlignedDir)
                altAlignedDir="$2"
                shift; shift
                ;;
            --genomeBaseRepo)
                genomeBaseRepo=$2
                shift; shift
                ;;
            --genomeBaseCustom)
                genomeBaseCustom=$2
                shift; shift
                ;;
            --fileBase)
                fileBase="$2"
                shift; shift
                ;;
            --exptID)
                exptID="$2"
                shift; shift
                ;;
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --endType)
                endType="$2"
                shift; shift
                ;;
            --readLength)
                readLength=$2
                shift; shift
                ;;
            --species)
                species="$2"
                shift; shift
                ;;
            --PrimerLengthMapFile)
                PrimerLengthMapFile="$2"
                shift; shift
                ;;
            --primerIdx2p)
                primerIdx2p="$2"
                shift; shift
                ;;
            --primerIdxpC)
                primerIdxpC="$2"
                shift; shift
                ;;
            --minL2p)
                minL2p="$2"
                shift; shift
                ;;
            --minLpC)
                minLpC="$2"
                shift; shift
                ;;
            --mergeLanes)
                mergeLanesFlag="--mergeLanes"
                shift;
                ;;
            --mindist)
                mindist="$2"
                shift; shift
                ;;
            --ncells)
                ncells="$2"
                shift; shift
                ;;
            --addSoftClipped)
                addSoftClippedFlag="--addSoftClipped"
                shift;
                ;;
            --priority)
                priority="$2"
                shift; shift
                ;;
            --focusTag)
                focusTag="$2"
                shift; shift
                ;;
            --ncores)
                ncores=$2
                shift; shift
                ;;
            --verbose)
                verbose="true"
                verboseFlag="--verbose"
                shift;
                ;;
            --debug)
                debug="$2"
                shift; shift
                ;;
            *)
                echo "doNewPipe --exptID 605 --sampleID scCLTdegenNucXXX [--baseDir Data] [--exptDir E.chex] [--rawDir raw] [--initDir init] [--analyzedDir analyzed] [--alignedDir star] [--altAlignedDir starAlt] [--genomeBaseRepo /lab/repo/resources/lib/star/2.4.0h1] [--genomeBaseCustom Data/Database/Indices/STAR] [--fileBase star.primaryNoDup] [--endType PE] [--readLength 75] [--species human] [--PrimerLengthMapFile src/lib/ChexPrimerTable.conf] [--ncores 5] [--primerIdx2p 505] [--primerIdxpC 302] [--minL2p 6] [--minLpC 6] [--addSoftClipped] [--mindist 20] [--ncells 1] [--debug info] [--verbose]"
                return 1
                ;;
        esac
    done

    ## Initialize repo
    echo "initNewPipe --baseDir $baseDir --exptName E.$exptID --pooledExptName $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --rawDir $rawDir --initDir $initDir --seqType CHEXSeq $mergeLanesFlag --sampleID $sampleID --endType $endType --PrimerLengthMapFile $PrimerLengthMapFile --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --debug $debug $verboseFlag"
    initNewPipe --baseDir $baseDir --exptName E.$exptID --pooledExptName $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --rawDir $rawDir --initDir $initDir --seqType CHEXSeq $mergeLanesFlag --sampleID $sampleID --endType $endType --PrimerLengthMapFile $PrimerLengthMapFile --primerIdx2p $primerIdx2p --primerIdxpC $primerIdxpC --minL2p $minL2p --minLpC $minLpC --debug $debug $verboseFlag

    ## Annotate primer quality in CHEXTRIM
    echo "parseChexTrimAnnot --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --chextrimDir $chextrimDir --sampleID $sampleID --endType $endType --PrimerLengthMapFile $PrimerLengthMapFile --debug $debug $verboseFlag"
    parseChexTrimAnnot --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --chextrimDir $chextrimDir --sampleID $sampleID --endType $endType --PrimerLengthMapFile $PrimerLengthMapFile --debug $debug $verboseFlag
    echo "parseChexTrimStats --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --chextrimDir $chextrimDir --sampleID $sampleID --endType $endType --debug $debug $verboseFlag"
    parseChexTrimStats --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --chextrimDir $chextrimDir --sampleID $sampleID --endType $endType --debug $debug $verboseFlag
    
    ## Separate reads into A, B, C, D and sub-categories
    echo "filterReadQual --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --sampleID $sampleID --endType $endType --inFileBase $fileBase --outFileBase $fileBase --debug $debug $verboseFlag"
    filterReadQual --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --sampleID $sampleID --endType $endType --inFileBase $fileBase --outFileBase $fileBase --debug $debug $verboseFlag
    
    ## Read stats for {qual x outputInPair}
    echo "doSamFlagStats --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase --sampleID $sampleID $verboseFlag"
    doSamFlagStats --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase --sampleID $sampleID $verboseFlag
    
    ## Merge subclasses as well as major classes A, B and C (A and B for the stringent criterion, or A, B and C for less stringent)
    echo "mergeQualOutputInPairs --sampleID $sampleID --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase --qualOutputInPairs 'A1either,A2either B1either,B2either C1either,C2either' --qualOutputInPairsMerged 'Aeither Beither Ceither' $verboseFlag"
    mergeQualOutputInPairs --sampleID $sampleID --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase --qualOutputInPairs 'A1either,A2either B1either,B2either C1either,C2either' --qualOutputInPairsMerged 'Aeither Beither Ceither' $verboseFlag
    
    ## Remove duplicate priming from class A, B and C 
    echo "removeAdjReads --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --sampleID $sampleID --fileBase $fileBase --endType $endType $addSoftClippedFlag --mindist $mindist --ncells $ncells --priority $priority --focusTag $focusTag $verboseFlag --debug $debug"
    removeAdjReads --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --sampleID $sampleID --fileBase $fileBase --endType $endType $addSoftClippedFlag --mindist $mindist --ncells $ncells --priority $priority --focusTag $focusTag $verboseFlag --debug $debug
    
    ## BAM to BED, keep tracking of the read orientation so that the strandedness is preserved in BED 
    echo "doBamToBed --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --sampleID $sampleID --endType $endType --fileBase $fileBase $addSoftClippedFlag $verboseFlag"
    doBamToBed --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --sampleID $sampleID --endType $endType --fileBase $fileBase $addSoftClippedFlag $verboseFlag
    
    ## Remove ENCODE blacklisted regions from BED
    echo "removeBlacklisted --species $species --sampleID $sampleID --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase $verboseFlag"
    removeBlacklisted --species $species --sampleID $sampleID --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase $verboseFlag
    
    ## Align reads to alternative genomes and detect contaminants
    echo "alignAltSpecies --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --trimDir $trimDir --alignedDir $alignedDir --altAlignedDir $altAlignedDir --genomeBaseRepo $genomeBaseRepo --genomeBaseCustom $genomeBaseCustom --sampleID $sampleID --endType $endType --readLength $readLength --species $species $verboseFlag"
    alignAltSpecies --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --trimDir $trimDir --alignedDir $alignedDir --altAlignedDir $altAlignedDir --genomeBaseRepo $genomeBaseRepo --genomeBaseCustom $genomeBaseCustom --sampleID $sampleID --endType $endType --readLength $readLength --species $species $verboseFlag

    ## Get alignment quality in terms of mapped read length and number of mismatches
    echo "getMappedLengthMismatches --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase --sampleID $sampleID $verboseFlag"
    getMappedLengthMismatches --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase $fileBase --sampleID $sampleID $verboseFlag

    ##  Extract read IDs with good mapping quality
    echo "filterMappedLengthMismatches --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase star.primaryNoDup.mappedLength_mismatches --mlen 20 --mm 0.1 --sampleID $sampleID $verboseFlag"
    filterMappedLengthMismatches --baseDir $baseDir --exptDir $exptDir --analyzedDir $analyzedDir --alignedDir $alignedDir --fileBase star.primaryNoDup.mappedLength_mismatches --mlen 20 --mm 0.1 --sampleID $sampleID $verboseFlag
}
###########################################################################
export -f annotateReadQual
export -f filterReadQual
export -f mergeQualOutputInPairs
export -f removeAdjReads 
export -f doBamToBed
export -f filterBed
export -f filterBedBySpecies
export -f removeBlacklisted
export -f doSamFlagStats
export -f parseChexTrimAnnot
export -f parseChexTrimStats
export -f summarizeBamFlagStats
export -f summarizeBamFlagQualStats
export -f initNewPipe
export -f alignAltSpecies
export -f getMappedLengthMismatches
export -f filterMappedLengthMismatches
export -f doNewPipe
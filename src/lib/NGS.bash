## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
## v1.0

## User Picard to detect and mark PCR duplicates
## input: BAM file
## output: BAM file
function markDup {
    baseDir="Data"
    exptName="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    inFileBase="star.posSorted"
    outFileBase="star.DupMarked"
    index=false
    verbose=false
    parallelGCThreads=5
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
            --inFileBase)
                inFileBase="$2"
                shift; shift
                ;;
            --outFileBase)
                outFileBase="$2"
                shift; shift
                ;;
            --index)
                index=true;
                shift
                ;;
            --verbose)
                verbose=true
                shift;
                ;;
            --parallelGCThreads)
                parallelGCThreads="$2"
                shift; shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            *)
                echo "markDup [--baseDir Data] [--exptName E.chex] [--analyzedDir analyzed] [--alignedDir star] [--inFileBasestar.posSorted] [--outFileBase star.uniquePosSorted] [--verbose] [--index] [--parallelGCThreads 5] [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done

    if [[ -z $sampleIDs ]]; then 
        sampleIDs=($(echo $baseDir/$exptName/raw/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi

    inDir=$baseDir/$exptName/$analyzedDir
    for sampleID in ${sampleIDs[@]}; do
        sampleName=Sample_$sampleID
        bamFileIn=$inDir/$sampleName/$alignedDir/$sampleName.$inFileBase.bam
        bamFileOut=$inDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam
        metricFile=$inDir/$sampleName/$alignedDir/$sampleName.$outFileBase.stats.txt
        logFile=$inDir/$sampleName/$alignedDir/$sampleName.$outFileBase.log
        logIndex=$inDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam.bai.log
        if [[ "$verbose" == "true" ]]; then
            ## Note, -XX:ParallelGCThreads=<Num> cannot be placed after -jar.
            echo "java -XX:ParallelGCThreads=$parallelGCThreads -jar $HOME/Applications/Picard/2.17.0/picard.jar MarkDuplicates INPUT=$bamFileIn OUTPUT=$bamFileOut METRICS_FILE=$metricFile REMOVE_DUPLICATES=false USE_JDK_DEFLATER=true USE_JDK_INFLATER=true 2> $logFile"
        fi
        java -XX:ParallelGCThreads=$parallelGCThreads -jar $HOME/Applications/Picard/2.17.0/picard.jar MarkDuplicates INPUT=$bamFileIn OUTPUT=$bamFileOut METRICS_FILE=$metricFile REMOVE_DUPLICATES=false USE_JDK_DEFLATER=true USE_JDK_INFLATER=true 2> $logFile
        if [[ "$index" == "true" ]]; then
            if [[ "$verbose" == "true" ]]; then
                echo "samtools index $bamFileOut 2> $logIndex"
            fi
            samtools index $bamFileOut 2> $logIndex
        fi
    done
}

## Sort the BAM by queryname or position
## input: BAM file
## output: BAM file (, BAI file if --index)
function sortUniqBam {
    baseDir="Data"
    exptName="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    inFileBase="star.posSorted"
    outFileBase="star.uniquePosSorted"
    mode="FromPosSorted"
    ncores=5
    verbose=false
    index=false
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
            --inFileBase)
                inFileBase="$2"
                shift; shift
                ;;
            --outFileBase)
                outFileBase="$2"
                shift; shift
                ;;
            --mode)
                mode="$2"
                shift; shift
                ;;
            --ncores)
                ncores="$2"
                shift; shift
                ;;
            --index)
                index=true
                shift
                ;;
            --verbose)
                verbose=true
                shift;
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            *)
                echo "sortUniqBam [--baseDir Data] [--exptName E.chex] [--analyzedDir analyzed] [--alignedDir star] [--inFileBasestar.posSorted] [--outFileBase star.uniquePosSorted] [--mode FromPosSorted] [--ncores 5] [--index] [--verbose] [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done

    if [[ -z $sampleIDs ]]; then 
        sampleIDs=($(echo $baseDir/$exptName/raw/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi

    inDir=$baseDir/$exptName/$analyzedDir
    for sampleID in ${sampleIDs[@]}; do
        sampleName=Sample_$sampleID
        bamFileIn=$inDir/$sampleName/$alignedDir/${sampleName}.${inFileBase}.bam
        bamFileOut=$inDir/$sampleName/$alignedDir/${sampleName}.${outFileBase}.bam
        logIndex=$inDir/$sampleName/$alignedDir/${sampleName}.${outFileBase}.bam.bai.log
        logFile=$inDir/$sampleName/$alignedDir/${sampleName}.${outFileBase}.log
        if [[ "$mode" == "FromPosSorted" ]]; then
            if [[ "$verbose" == "true" ]]; then 
                echo "samtools view -h -F 0x4 $bamFileIn | perl -ane 'print if (\$F[0] =~ /^@/ or \$_ =~ /\"NH:i:1\t/\")' | samtools view -Sb - > $bamFileOut 2>$logFile"
            fi
            samtools view -h -F 0x4 $bamFileIn | perl -ane 'print if ($F[0] =~ /^@/ or $_ =~ /NH:i:1\t/)' | samtools view -Sb - > $bamFileOut 2>$logFile
        elif [[ "$mode" == "FromUnique" ]]; then
            if [[ "$verbose" == "true" ]]; then 
                echo "samtools sort -@ $ncores -O bam -T sortTmp -o $bamFileOut $bamFileIn 2>$logFile"
            fi
            samtools sort -@ $ncores -O bam -T sortTmp -o $bamFileOut $bamFileIn 2>$logFile
        fi
        if [[ "$index" == "true" ]]; then
            if [[ "$verbose" == "true" ]]; then
                echo "samtools index $bamFileOut 2> $logIndex"
            fi
            samtools index $bamFileOut 2> $logIndex
        fi
    done
}

## Detect and mark PCR duplicates in uniquely mapped reads
## input: BAM file
## output: BAM file (, BAI file if --index)
function markDupUniq {
    analyzedDir="analyzed"
    verboseFlag=""
    while [[ $# -gt 0 ]]; do 
        case $1 in 
            --sampleID)
                sampleID="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift;
                ;;
            --verbose)
                verboseFlag="--verbose"
                shift;
                ;;
            *)
                echo "markDupUniq --sampleID scCLTdegenNucxxx [--analyzedDir analyzed] [--verbose]"
                return 1
                ;;
        esac
    done
    if [[ "$verbose" == "true" ]]; then
        echo "markDup --baseDir Data --exptName E.chex --analyzedDir $analyzedDir --alignedDir star --inFileBase star.posSorted --outFileBase star.DupMarked $verboseFlag --index --sampleIDs $sampleID"
    fi
    markDup --baseDir Data --exptName E.chex --analyzedDir $analyzedDir --alignedDir star --inFileBase star.posSorted --outFileBase star.DupMarked $verboseFlag --index --sampleIDs $sampleID
    if [[ "$verbose" == "true" ]]; then
        echo "rmDup --baseDir Data --exptName E.chex --analyzedDir $analyzedDir --alignedDir star --inFileBase star.DupMarked --outFileBase star.NoDup $verboseFlag --index --sampleIDs $sampleID"
    fi
    rmDup --baseDir Data --exptName E.chex --analyzedDir $analyzedDir --alignedDir star --inFileBase star.DupMarked --outFileBase star.NoDup $verboseFlag --index --sampleIDs $sampleID
    if [[ "$verbose" == "true" ]]; then
        echo "sortUniqBam --baseDir Data --exptName E.chex --analyzedDir $analyzedDir --alignedDir star --inFileBase star.NoDup --outFileBase star.uniqueNoDup --mode FromPosSorted $verboseFlag --index --sampleIDs $sampleID"
    fi
    sortUniqBam --baseDir Data --exptName E.chex --analyzedDir $analyzedDir --alignedDir star --inFileBase star.NoDup --outFileBase star.uniqueNoDup --mode FromPosSorted $verboseFlag --index --sampleIDs $sampleID
}

## Filter BAM for reads of particular flag(s)
## input: BAM file
## output: BAM file (, BAI file if --index)
function selectReads {
    baseDir="Data"
    exptName="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    inFileBase=""
    outFileBase=""
    verbose=false
    flags="-F 0x4" # by default we just keep the mapped reads
    index=false
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
            --inFileBase)
                inFileBase="$2"
                shift; shift
                ;;
            --outFileBase)
                outFileBase="$2"
                shift; shift
                ;;
            --flags)
                flags="$2"
                shift; shift
                ;;
            --index)
                index=true
                shift
                ;;
            --verbose)
                verbose=true
                shift;
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            *)
                echo "selectReads [--baseDir Data] [--flags "-F 0x4"] [--exptName E.chex] [--analyzedDir analyzed] [--alignedDir star] [--inFileBasestar.posSorted] [--outFileBase star.uniquePosSorted] [--verbose] [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done

    if [[ -z $sampleIDs ]]; then 
        sampleIDs=($(echo $baseDir/$exptName/raw/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi

    if [[ "$inFileBase" == "" ]]; then 
        echo "inFileBase not set!"
        return 2
    fi

    if [[ "$inFileBase" == "$outFileBase" ]]; then
        echo "inFileBase identical to outFileBase!"
        return 3
    fi

    inDir=$baseDir/$exptName/$analyzedDir
    for sampleID in ${sampleIDs[@]}; do 
        sampleName=Sample_$sampleID
        bamFileIn=$inDir/$sampleName/$alignedDir/${sampleName}.${inFileBase}.bam
        bamFileOut=$inDir/$sampleName/$alignedDir/${sampleName}.${outFileBase}.bam
        logFile=$inDir/$sampleName/$alignedDir/${sampleName}.${outFileBase}.bam.log
        logIndex=$inDir/$sampleName/$alignedDir/${sampleName}.${outFileBase}.bam.bai.log
        if [[ "$verbose" == "true" ]]; then 
            echo "samtools view -h $flags $bamFileIn | samtools view -Sb - > $bamFileOut 2> $logFile"
        fi
        samtools view -h $flags $bamFileIn | samtools view -Sb - > $bamFileOut 2> $logFile
        if [[ "$index" == "true" ]]; then 
            if [[ "$verbose" == "true" ]]; then 
                echo "samtools index $bamFileOut 2>$logIndex"
            fi
            samtools index $bamFileOut 2> $logIndex
        fi
    done
}

## Filter out PCR duplicates from BAM
## input: BAM file
## output: BAM file
function rmDup {
    selectReads --flags '-F 0x400' $@
}

## Insert a line "@RG     ID:id   LB:library      PL:platform     SM:sample       PU:machine" to BAM
## input: BAM file
## output: BAM file
function addReadGroup {
    baseDir="Data"
    exptName="E.chex"
    analyzedDir="analyzed"
    alignedDir="star"
    rawDir="raw"
    inFileBase="star.posSorted"
    outFileBase="star.RGAdded"
    SO="coordinate"
    RGID="id"
    RGLB="library"
    RGPL="platform"
    RGPU="machine"
    RGSM="sample"
    verbose=false
    parallelGCThreads=5
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
            --rawDir)
                rawDir="$2"
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
            --verbose)
                verbose=true
                shift;
                ;;
            --parallelGCThreads)
                parallelGCThreads="$2"
                shift; shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            *)
                echo "addReadGroup [--baseDir Data] [--exptName E.chex] [--analyzedDir analyzed] [--alignedDir star] [--inFileBasestar.posSorted] [--outFileBase star.RGAdded] [--verbose] [--index] [--parallelGCThreads 5] [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done

    if [[ -z $sampleIDs ]]; then 
        sampleIDs=($(echo $baseDir/$exptName/$analyzedDir/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi

    inDir=$baseDir/$exptName/$analyzedDir
    for sampleID in ${sampleIDs[@]}; do
        sampleName=Sample_$sampleID
        bamFileIn=$inDir/$sampleName/$alignedDir/$sampleName.$inFileBase.bam
        bamFileOut=$inDir/$sampleName/$alignedDir/$sampleName.$outFileBase.bam
        logFile=$inDir/$sampleName/$alignedDir/$sampleName.$outFileBase.log
        if [[ "$verbose" == "true" ]]; then 
            ## Note, -XX:ParallelGCThreads=<Num> cannot be placed after -jar.
            echo "java -XX:ParallelGCThreads=$parallelGCThreads -jar $HOME/Applications/Picard/2.17.0/picard.jar AddOrReplaceReadGroups I=$bamFileIn O=$bamFileOut SO=$SO RGID=$RGID RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM 2> $logFile"
        fi
        java -XX:ParallelGCThreads=$parallelGCThreads -jar $HOME/Applications/Picard/2.17.0/picard.jar AddOrReplaceReadGroups I=$bamFileIn O=$bamFileOut SO=$SO RGID=$RGID RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM 2> $logFile
    done
}
###########################################################################
export -f markDup
export -f sortUniqBam
export -f markDupUniq
export -f selectReads
export -f rmDup

## Youtao Lu <luyoutao@sas.upenn.edu>
## Copyright (c) 2017-2022, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2022, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
###########################################################################
## v1.1

## The CHEX-seq data structure is based on the SCAP-T NGS pipeline by Fisher et al.
## https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/document.cgi?study_id=phs000833.v1.p1&phd=4779
## For more details refer to https://github.com/safisher/ngs/wiki.
## 
## The basic data structure consists of the folders:
## 
## E.XXX
## |--src
##    |--[FLOWCELL_ID]
## |--raw
##    |--Reports
##    |--Sample_XXX
##       |--XXX_1_S1_L001_R1_001.fastq.gz
##       |--XXX_1_S1_L001_R2_001.fastq.gz
##       |--XXX_1_S1_L002_R1_001.fastq.gz
##       |--XXX_1_S1_L002_R2_001.fastq.gz
##       |--...
##    |--Sample_YYY
##       |--XXX_1_S1_L001_R1_001.fastq.gz
##       |--...
##    |--...
##    |--Stats
## |--info
##    |--config.csv
##    |--SampleSheet.csv
##    |--...
## |--analyzed
##    |--Sample_XXX
##       |--blast
##       |--fastqc
##       |--trim
##       |--fastqc.trim
##       |--star
##       |--verse
##       |--log
## 
## 
## CHEX-seq extends the above basic folders by adding these specific ones:
## |--analyzed
##    |--Sample_XXX
##       |--rmdup
##       |--fastqc.rmdup
##       |--chextrim
##       |--fastqc.chextrim
##       |--starAlt

## Create a folder mirroring the root folder /lab/repo in the end user's space
## input: the KimLab /lab/repo/E.xxx folder
## output: $HOME/Datasets/lab/repo/E.xxx
function symlinkDir {
    # Note, by default, bash variables are non-local, which means outside variables are accessible within stack calls--definitely not going to work for recursion.
    # We need to set them `local`.
    local inDir=""
    local outDir="./"
    hidden=false
    verbose=false
    while [[ $# -gt 0 ]]; do
        case $1 in
            --inDir)
                inDir="$2"
                shift; shift
                ;;
            --outDir)
                outDir="$2"
                shift; shift
                ;;
            --hidden)
                hidden=true
                shift;
                ;;
            --verbose)
                verbose=true
                shift;
                ;;
            *)
                echo "symlinkDir --inDir src --outDir dest [--hidden] [--verbose]"
                return 1;
                ;;
        esac
    done
    if [[ "$inDir" == "" ]] || ! [[ -e "$inDir"  ]]; then
        echo "$inDir does not exist!"
        return 2;
    fi
    mkdir -p $outDir

    if [[ "$hidden" == "true" ]]; then
        Files=$(/bin/ls -aF $inDir)
    else 
        Files=$(/bin/ls -F $inDir)
    fi
    for File in $Files; do
        if [[ "$File" == './' ]] || [[ "$File" == '../' ]]; then 
            continue
        fi
        echo $File | grep -q -E '/$'
        if [[ $? -eq 0 ]]; then # this is a directory
            File=$(echo $File | sed 's/\/$//')
            local subInDir=$inDir/$File  
            local subOutDir=$outDir/$File
            if [[ "$hidden" == "true" ]]; then
                if [[ "$verbose" == "true" ]]; then
                    echo "symlinkDir --inDir $subInDir --outDir $subOutDir --hidden --verbose"
                    symlinkDir --inDir $subInDir --outDir $subOutDir --hidden --verbose
                else 
                    symlinkDir --inDir $subInDir --outDir $subOutDir --hidden
                fi
            else
                if [[ "$verbose" == "true" ]]; then
                    echo "symlinkDir --inDir $subInDir --outDir $subOutDir --verbose"
                    symlinkDir --inDir $subInDir --outDir $subOutDir --verbose
                else 
                    symlinkDir --inDir $subInDir --outDir $subOutDir
                fi
            fi
        else 
            File=$(basename $File)
            if [[ "$verbose" == "true" ]]; then
                echo "ln -sf $inDir/$File $outDir/$File"
            fi
            ln -sf $inDir/$File $outDir/$File
        fi
    done
}

## Symbolic link /lab/repo/E.xxx/{src,info}/* to ~/Datasets/lab/repo/E.xxx/{src,info}/*
## Symbolic link /lab/repo/E.xxx/raw/Sample_xxx/* to ~/Datasets/lab/repo/E.xxx/raw/Sample_xxx/*
## Symbolic link /lab/repo/E.xxx/analyzed/Sample_xxx/{blast,chextrim,fastqc,fastqc.chextrim,fastqc.rmdup,fastqc.trim,log,rmdup,star,trim,verse}/* to ~/Datasets/lab/repo/E.xxx/analyzed/Sample_xxx/{blast,chextrim,fastqc,fastqc.chextrim,fastqc.rmdup,fastqc.trim,log,rmdup,star,trim,verse}/*
function initRepo {
    inBaseDir="/lab/repo"
    outBaseDir="Data/Datasets/lab/repo"
    analyzedDir="analyzed"
    rawDir="raw"
    seqType="CHEXSeq"
    verbose="false"
    declare -A seqTypeSubfolders=([WGS]="blast bowtie fastqc fastqc.trim log trim" [RNASeq]="blast fastqc fastqc.trim log star trim verse" [CHEXSeq]="blast chextrim fastqc fastqc.chextrim fastqc.rmdup fastqc.trim log rmdup star trim verse")
    ## star-strict removed from seqTypeSubfolders
    ## default subfolders is set for CHeX-seq module
    while [[ $# -gt 0 ]]; do
        case $1 in
            --inBaseDir) 
                inBaseDir="$2"
                shift; shift
                ;;
            --outBaseDir)
                outBaseDir="$2"
                shift; shift
                ;;
            --analyzedDir)
                analyzedDir="$2"
                shift; shift
                ;;
            --rawDir)
                rawDir="$2"
                shift; shift
                ;;
            --exptName)
                exptName="$2"
                shift; shift
                ;;
            --seqType)
                seqType="$2"
                shift; shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            *)
                echo "initRepo [--inBaseDir /lab/repo] [--outBaseDir ~/Datasets/lab/repo] --exptName E.xxx [--analyzedDir analyzed] [--rawDir raw] [--seqType CHEXSeq] [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done
    if [[ -z $sampleIDs ]]; then 
        sampleIDs=($(echo /lab/repo/$exptName/$analyzedDir/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi
    if [[ "$verbose" == "true" ]]; then 
        echo "initRepo --inBaseDir $inBaseDir --outBaseDir $outBaseDir --exptName $exptName --analyzedDir $analyzedDir --seqType $seqType --sampleIDs ${sampleIDs[@]}"
    fi

    for folder in info src; do 
        if [[ "$verbose" == "true" ]]; then
            echo "mkdir -p $outBaseDir/$exptName/$folder"
        fi 
        mkdir -p $outBaseDir/$exptName/$folder
        files=$(ls --color=no $inBaseDir/$exptName/$folder)
        for file in $files; do
            if [[ "$verbose" == "true" ]]; then
                echo "ln -fs $inBaseDir/$exptName/$folder/$file $outBaseDir/$exptName/$folder/$file"
            fi 
            ln -fs $inBaseDir/$exptName/$folder/$file $outBaseDir/$exptName/$folder/$file
        done
    done
    for sampleID in ${sampleIDs[@]}; do 
        sampleName=Sample_$sampleID
        if [[ "$verbose" == "true" ]]; then
            echo "mkdir -p $outBaseDir/$exptName/$rawDir/$sampleName"
        fi
        files=$(ls --color=no $inBaseDir/$exptName/$rawDir/$sampleName)
        for file in $files; do
            mkdir -p $outBaseDir/$exptName/$rawDir/$sampleName
            if [[ "$verbose" == "true" ]]; then
                echo "ln -fs $inBaseDir/$exptName/$rawDir/$sampleName/$file $outBaseDir/$exptName/$rawDir/$sampleName/$file"
            fi
            ln -fs $inBaseDir/$exptName/$rawDir/$sampleName/$file $outBaseDir/$exptName/$rawDir/$sampleName/$file
        done
        for folder in ${seqTypeSubfolders[$seqType]}; do 
            if [[ "$verbose" == "true" ]]; then
                echo "mkdir -p $outBaseDir/$exptName/$analyzedDir/$sampleName/$folder"
            fi
            mkdir -p $outBaseDir/$exptName/$analyzedDir/$sampleName/$folder
            files=$(ls --color=no $inBaseDir/$exptName/$analyzedDir/$sampleName/$folder)
            for file in $files; do
                if [[ "$verbose" == "true" ]]; then
                    echo "ln -fs $inBaseDir/$exptName/$analyzedDir/$sampleName/$folder/$file $outBaseDir/$exptName/$analyzedDir/$sampleName/$folder/$file"
                fi
                ln -fs $inBaseDir/$exptName/$analyzedDir/$sampleName/$folder/$file $outBaseDir/$exptName/$analyzedDir/$sampleName/$folder/$file
            done
        done
    done
}

## Merge sequencing lanes of a sample, specificially, concatenate files 
##  ~/Datasets/lab/repo/E.xxx/raw/Sample_xxx/{xxx_S3_L001_R1_001.fastq.gz,xxx_S3_L002_R1_001.fastq.gz,xxx_S3_L003_R1_001.fastq.gz,xxx_S3_L004_R1_001.fastq.gz}
## into ~/Datasets/lab/repo/E.xxx/raw/Sample_xxx/unaligned_1.fq.gz
## If paired-end sequencing, same processing for R2.
function mergeRawLanes {
    baseDir="Data/Datasets/lab/repo"
    rawDir="raw"
    endType="PE"
    verbose="false"
    while [ $# -gt 0 ]; do
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
            --endType)
                endType="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            *)
                echo "mergeRawLanes [--baseDir Data/Datasets/lab/repo] --exptName E.xxx [--rawDir raw] [--endType PE] [--verbose] [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done

    if [[ "$endType" == "PE" ]]; then
        isPE=true
    elif [[ "$endType" == "SE" ]]; then
        isPE=false
    else
        echo "endType should be either PE or SE!"
        return 2
    fi

    if [[ -z $sampleIDs ]]; then 
        sampleIDs=($(echo $baseDir/$exptName/$analyzedDir/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi
    if [[ "$verbose" == "true" ]]; then 
        echo mergeRawLanes --baseDir $baseDir --rawDir $rawDir --endType $endType --verbose $verbose --exptName $exptName --sampleIDs ${sampleIDs[@]} 
    fi
    for sampleID in ${sampleIDs[@]}; do 
        sampleName=Sample_$sampleID
        Dir=$baseDir/$exptName/$rawDir/$sampleName
        R1out=$Dir/unaligned_1.fq.gz
        R1log=$Dir/unaligned_1.log
        if [[ "$verbose" == "true" ]]; then
            echo "cat $Dir/${sampleID}_*_L*_R1_001.fastq.gz 1>$R1out 2>$R1log"
        fi
        cat $Dir/${sampleID}_*_L*_R1_001.fastq.gz 1>$R1out 2>$R1log
        if [[ "$isPE" == "true" ]]; then
            R2out=$Dir/unaligned_2.fq.gz
            R2log=$Dir/unaligned_2.log
            if [[ "$verbose" == "true" ]]; then
                echo "cat $Dir/${sampleID}_*_L*_R2_001.fastq.gz 1>$R2out 2>$R2log"
                cat $Dir/${sampleID}_*_L*_R2_001.fastq.gz 1>$R2out 2>$R2log
            fi
        fi
    done
}

## Symlink raw folder to init inside analyzed folder, here all exptName folders have been pooled together already
## ../../../raw/Sample_xxx/{unaligned_1.fq.gz, unaligned_2.fq.gz} to Data/Datasets/lab/repo/E.123/analyzed/Sample_xxx/init/{unaligned_1.fq.gz, unaligned_2.fq.gz}
function symlinkRawToAnalyzedInit {
    baseDir="Data/Datasets/lab/repo"
    exptName="E.chex"
    rawDir="raw"
    analyzedDir="analyzed"
    initDir="init"
    inRawFqHasSampleID="false"
    inRawFqBase="unaligned_1.fq.gz"
    outInitFqBase="unaligned_1.fq.gz"
    verbose="false"
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
            --initDir)
                initDir="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            --inRawFqHasSampleID)
                inRawFqHasSampleID="true"
                shift
                ;;
            --inRawFqBase)
                inRawFqBase="$2"
                shift; shift
                ;;
            --outInitFqFile)
                outInitFqFile="$2"
                shift; shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            *)
                echo "symlinkRawToAnalyzedInit [--baseDir Data/Datasets/lab/repo] [--exptName E.chex] [--rawDir raw] [--analyzedDir analyzed] [--initDir init] [--inRawFqHasSampleID] [--inRawFqBase _R1_001.fastq.gz] [--outInitFqFile unaligned_1.fq.gz] [--verbose] [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done

    if [[ -z $sampleIDs ]]; then 
        sampleIDs=($(echo $baseDir/$exptName/$analyzedDir/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi

    if [[ "$verbose" == "true" ]]; then 
        echo symlinkRawToAnalyzedInit --baseDir $baseDir --exptName $exptName --rawDir $rawDir --analyzedDir $analyzedDir --initDir $initDir --verbose $verbose --sampleIDs ${sampleIDs[@]}
    fi

    for sampleID in ${sampleIDs[@]}; do
        sampleName=Sample_$sampleID
        mkdir -p $baseDir/$exptName/$analyzedDir/$sampleName/$initDir
        if [[ "$verbose" == "true" ]]; then 
            echo "mkdir -p $baseDir/$exptName/$analyzedDir/$sampleName/$initDir"
            if [[ "$inRawFqHasSampleID" == "true" ]]; then
                echo "ln -fs ../../../$rawDir/$sampleName/$sampleID$inRawFqBase $baseDir/$exptName/$analyzedDir/$sampleName/$initDir/$outInitFqFile"
            else 
                echo "ln -fs ../../../$rawDir/$sampleName/$inRawFqBase $baseDir/$exptName/$analyzedDir/$sampleName/$initDir/$outInitFqFile"
            fi
        fi
        if [[ "$inRawFqHasSampleID" == "true" ]]; then
            ln -fs ../../../$rawDir/$sampleName/$sampleID$inRawFqBase $baseDir/$exptName/$analyzedDir/$sampleName/$initDir/$outInitFqFile
        else 
            ln -fs ../../../$rawDir/$sampleName/$inRawFqBase $baseDir/$exptName/$analyzedDir/$sampleName/$initDir/$outInitFqFile
        fi
    done
}

## Pool individual exptName (e.g. E.123) raw folders into a common E.project (e.g. E.chex) directory by symlink
## inside current project's working directory:
## ../../E.123/raw/Sample_xxx to Data/Datasets/lab/repo/E.chex/raw/Sample_xxx
function poolRawDir {
    baseDir=Data/Datasets/lab/repo
    pooledExptName="E.chex"
    rawDir="raw"
    verbose="false"
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir="$2"
                shift; shift
                ;;
            --pooledExptName)
                pooledExptName="$2"
                shift; shift
                ;;
            --rawDir)
                rawDir="$2"
                shift; shift
                ;;
            --exptName)
                exptName="$2"
                shift; shift
                ;;
            --verbose)
                verbose="true"
                shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            *)
                echo "poolRawDir [--baseDir Data/Datasets/lab/repo] [--pooledExptName E.proj] [--rawDir raw] [--verbose] --exptName E.xxx [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done

    if [[ -z $sampleIDs ]]; then
        sampleIDs=($(echo $baseDir/$exptName/$analyzedDir/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi

    if [[ "$verbose" == "true" ]]; then
        echo "mkdir -p $baseDir/$pooledExptName/$rawDir"
    fi
    mkdir -p $baseDir/$pooledExptName/$rawDir

    for sampleID in ${sampleIDs[@]}; do
        sampleName=Sample_$sampleID
        if [[ "$verbose" == "true" ]]; then
            echo "ln -fs ../../$exptName/$rawDir/$sampleName $baseDir/$pooledExptName/$rawDir/$sampleName"
        fi
        ln -fs ../../$exptName/$rawDir/$sampleName $baseDir/$pooledExptName/$rawDir/$sampleName
    done
}

## Pool individual ExptName analyzed folders into a common E.project directory by symlink. e.g.
## ../../E.123/analyzed/Sample_xxx to Data/Datasets/lab/repo/E.chex/analyzed/Sample_xxx
function poolAnalyzedDir {
    baseDir="Data/Datasets/lab/repo"
    pooledExptName="E.chex"
    analyzedDir="analyzed"
    verbose="false"
    seqType="CHEXSeq"
    declare -A seqTypeSubfolders=([WGS]="blast bowtie fastqc fastqc.trim log trim" [RNASeq]="blast fastqc fastqc.trim log star trim verse" [CHEXSeq]="blast chextrim fastqc fastqc.chextrim fastqc.rmdup fastqc.trim log rmdup star trim verse")
    ## star-strict removed from seqTypeSubfolder
    ## default subfolders is set for CHeX-seq module
    while [[ $# -gt 0 ]]; do
        case $1 in
            --baseDir)
                baseDir="$2"
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
            --verbose)
                verbose="true"
                shift;
                ;;
            --exptName)
                exptName="$2"
                shift; shift
                ;;
            --seqType)
                seqType="$2"
                shift; shift
                ;;
            --sampleIDs)
                shift
                sampleIDs=("${@}")
                break
                ;;
            *)
                echo "poolAnalyzedDir [--baseDir Data/Datasets/lab/repo] [--pooledExptName E.chex] [--verbose] --exptName E.xxx [--analyzedDir analyzed] [--sampleIDs sampleID1 sampleID2 ...]"
                return 1
                ;;
        esac
    done

    if [[ -z $sampleIDs ]]; then 
        sampleIDs=($(echo $baseDir/$exptName/$analyzedDir/Sample_* | sed -e 's/ /\n/g' | sed -e 's/.*\///g; s/Sample_//g;' | xargs))
    fi

    for sampleID in ${sampleIDs[@]}; do
        sampleName=Sample_$sampleID
        if [[ "$verbose" == "true" ]]; then 
            echo "mkdir -p $baseDir/$pooledExptName/$analyzedDir/$sampleName"
        fi
        mkdir -p $baseDir/$pooledExptName/$analyzedDir/$sampleName
        subdirs=$(ls --color=no $baseDir/$exptName/$analyzedDir/$sampleName)
        for subdir in $subdirs; do
            if [[ "$verbose" == "true" ]]; then 
                echo "ln -fs ../../../$exptName/$analyzedDir/$sampleName/$subdir $baseDir/$pooledExptName/$analyzedDir/$sampleName/$subdir"
            fi
            ln -fs ../../../$exptName/$analyzedDir/$sampleName/$subdir $baseDir/$pooledExptName/$analyzedDir/$sampleName/$subdir
        done
    done
}
###########################################################################
export -f initRepo
export -f mergeRawLanes
export -f symlinkRawToAnalyzedInit
export -f poolRawDir
export -f poolAnalyzedDir

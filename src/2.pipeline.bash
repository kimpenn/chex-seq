source src/lib/CHEX.bash

## 1. Run the SCAP-T general pipeline per 12 samples in paralell
parallel --progress -j 1 --pipe -N1 launcher_parallel.sh - 12 < data/E.chex/info/config.csv

## 2. Run CHEX-seq pipeline per 12 samples in paralell
sed 1d data/SampleInfo.csv | parallel --progress --colsep ',' -j 12 "doNewPipe --baseDir Data --exptDir E.chex --exptID {2} --sampleID {3} --endType {6} --primerIdx2p {14} --primerIdxpC {18} --PrimerLengthMapFile Data/20190801/ChexPrimerTable.conf --species {5} --ncells {26} --mergeLanes --ncores 5 --debug info --verbose"

## 3.1 Gather all samples' stats regardless of barcode/primer quality
summarizeBamFlagStats --fileBase star.primaryNoDup --alignedDir star --outFile results/ReadQualStats/overall.tsv --sampleIDs $(sed '1d' data/SampleInfo.csv | cut -d "," -f3 | tr "\n" " ")

## 3.2 Gather all samples' stats for each barcode/primer quality class-subclass
## Note, we don't need to read singletons, as they can be inferred from the overall stats. 
## for paired-end
summarizeBamFlagQualStats --fileBase star.primaryNoDup --alignedDir star --outFile results/ReadQualStats/perclass_PE.tsv --quals "A1,A2,B1,B2,C1,C2,D" --outputInPairs "read,mate,both,either" --verbose --sampleIDs $(sed '1d' data/SampleInfo.csv | grep ",PE," | cut -d "," -f3 | tr "\n" " ") 
## for single-end (if any)
summarizeBamFlagQualStats --fileBase star.primaryNoDup --alignedDir star --outFile results/ReadQualStats/perclass_SE.tsv --quals "A1,A2,B1,B2,C1,C2,D" --outputInPairs "read" --verbose --sampleIDs $(sed '1d' data/SampleInfo.csv | grep ",SE," | cut -d "," -f3 | tr "\n" " ")
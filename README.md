# CHEX-seq: CHromatin EXposed Sequencing for Single-Stranded Open-Chromatin Analysis

## Overview
Photoactivatable primer based Transcriptome In Vivo Analysis (TIVA) has shown successful application in profiling the transcriptome in live single cells [1]. As one of its successors, CHEX-seq employs 3' blocked, light-activatable probes which, upon annealing to the single-stranded DNA, is laser activated in selected cells (fixed or primary culture) at a specified foci, then serves as a primer for *in situ* copying of the single-stranded DNA into amplicons that can be finally sequenced. In this way, CHEX-seq enables a genome-wide map of the single-stranded open-chromatin landscape preserving the local tissue context.

This is the computational pipeline for converting the reads from next-generation sequencing to single-stranded chromatin sites (in BED format) as well as the gene-wise count matrix in different subgenic regions. 

This pipeline takes in the [PennSCAP-T NGS pipeline](https://github.com/kimpenn/ngs-pipeline) preprocessed data, evaluates the reads based on the barcode/primer quality as well as the mapping quality, then tracks the start site of a primer-extension event (i.e the primer annealing site) and removes multiple start sites in a small genomic range (20bp by default) (as biochemically unlikely). Next, it summarizes the primer annealing sites (as single spots) as well as the primer extending regions (as genomic intervals) in BED format. Finally, it searches and removes various technical artifacts, including the ENCODE low-mappability blacklists, the cross-species contaminants, and the low-mapping-quality reads. 

For complete manual please refer to the [Wiki Page](doc/wiki.md).

## Dependencies
* [PennSCAP-T pipeline](https://github.com/safisher/ngs) (>=2.3)
* GNU parallel (>=20190122)
* Perl (>=5.010), Log::Log4perl (>1.01)
* R (>=3.5.0), GenomicRanges (>=1.34), GenomicAlignments (>=1.18), rtracklayer (>=1.42), GenomeInfoDb (>=1.29), ensembldb (>=2.6), EnsDb.Hsapiens.v86 (>=2.99), EnsDb.Mmusculus.v79 (>=2.99)
* Picard (>=2.17)
* bedtools (>=2.26)

## Installation
1. Install the dependencies as listed above;
2. This package is designed to work as a `lib` hence doesn't need to be installed to your system or local `bin`. You just need to download and place the source code folder (`src`) to the workplace where you want to run the pipeline. Make sure the folder `data` is the same level as `src`. If you want to use some module as standalone, you may simply source the corresponding file in bash or R environment. 

## Usage
1. Prepare the Sequencing Configuration (`data/info/config.csv`), the sample sheet (`data/SampleInfo.csv`);
2. Generate the sample grouping info in R session
```
source("src/1.groupinfo.R")
``` 
2. Run the pipeline by
```
$ bash src/2.pipeline.bash
``` 
3. Get QC statistics in R session
```
source("src/3.qc.R")
```
4. Extract the priming sites of desired quality in R session
```
source("src/4.sites.R")
```
5. Generate genomic features and store for count quantification in R session
```
source("src/5.features.R")
```
6. Generate priming count matrices per genomic feature type in R session
```
source("src/6.cntmatrix.R")
```

## Copyright
```
Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
All Rights Reserved.
```

## License
[MPL-2](https://www.mozilla.org/en-US/MPL/2.0/)

## References
1. Ditte Lovatt, Brittani K Ruble, Jaehee Lee, Hannah Dueck, Tae Kyung Kim, Stephen Fisher, Chantal Francis, Jennifer M Spaethling, John A Wolf, M Sean Grady, Alexandra V Ulyanova, Sean B Yeldell, Julianne C Griepenburg, Peter T Buckley, Junhyong Kim, Jai-Yoon Sul, Ivan J Dmochowski and James Eberwine. "Transcriptome in Vivo Analysis (TIVA) of Spatially Defined Single Cells in Live Tissue." ([Nature Methods 11, no. 2 (February 2014): 190-96] (https://doi.org/10.1038/nmeth.2804))
2. Youtao Lu, Jaehee Lee, Jifen Li, Srini Rao Allu, Jinhui Wang, Hyun-Bum Kim, Kevin L. Bullaghey, Stephen A. Fisher, C. Erik Nordgren, Jean G. Rosario, Stewart A. Anderson, Alexandra V. Ulyanova, Steven Brem, H. Isaac Chen, John A. Wolf, M. Sean Grady, Sergei A. Vinogradov, Junhyong Kim and James Eberwine. "CHEX-seq Detects Single-Cell Single-Strand DNA With Catalytical Capacity." (under review)

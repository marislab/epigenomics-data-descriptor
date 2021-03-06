# epigenomics-data-descriptor
Scripts to call super enhancers, filter and generate ChIP-seq heatmaps for [neuroblastoma epigenomics data descriptor paper](https://www.biorxiv.org/content/10.1101/829754v2), and to calculate IDR.

[![DOI](https://zenodo.org/badge/208844258.svg)](https://zenodo.org/badge/latestdoi/208844258)

## Table of Contents
- [1. LILY (Super Enhancer calling)](#1.%20LILY%20(Super%20Enhancer%20calling))
- [2. ChipSeq Heatmaps](#2.%20ChipSeq%20Heatmaps)
- [3. IDR Analysis](#3.%20IDR%20Analysis)
- [Additional info](#Additional%20info)

<a name="1. LILY (Super Enhancer calling)"></a>
## 1. LILY (Super Enhancer calling)
LILY is a pipeline for detection of super-enhancers using H3K27ac ChIP-seq data, which includes explicit correction for copy number variation inherent to cancer samples. The pipeline is based on the ROSE algorithm originally developed by the the Young lab. 

### Before running
Follow steps 1-3 provided in the LILY's github documentation (https://github.com/BoevaLab/LILY). Clone the repository to run LILY scripts.

### Prerequisites
You will need following files to run LILY:
1. narrowPeak, regions.bed and .wig files for a particular sample should be present in the data folder
2. hg19_refseq.ucsc (transcriptome information which can be found at https://github.com/linlabbcm/rose2/tree/master/rose2/annotation)
3. hg19.chrom.sizes (file with chromosome lengths)

### How to run 
Script calls super enhances iteratively for all lines. Make sure to change data directory and result directory paths before running the script.
```R
Rscript lily.R
```

<a name="2. ChipSeq Heatmaps"></a>
## 2. ChipSeq Heatmaps
ChipSeq heatmaps were generated for MYCN (annotating top 5K peaks present in at least 5 mycn amplified cell lines) and all histone marks for COGN415 and NB69 line (annotating promoters, enhancers and super enhancers called from step 1)

### Prerequisites
[deepTools 3.2.0](https://deeptools.readthedocs.io/en/develop/content/installation.html)

### How to run
* To generate MYCN heatmaps (with top 5K peaks)
```R
Rscript makeHeatmaps.R
```

To generate heatmaps of histone-marks for COGN415 and NB69:
* To plot promoters
```R
Rscript promoterPlot.R
```
* To plot enhancers
```R
Rscript enhancerPlot.R
```
* To plot Super Enhancers
```R
Rscript SEPlot.R
```

<a name="3. IDR Analysis"></a>
## 3. IDR Analysis

The number of peaks passing the IDR threshold is calculated between true replicates (Nt), pooled pseudoreplicates (Np), self-pseudoreplicates1 (N1) and self-pseudoreplicates2 (N2).

* Nt: iqr_selfreps_noinput_merged_Final.sh
* Np: shuffle_bam_pooled_pseudoreps_noinput_Final.sh
* N1: shuffle_bam_pseudo_selfreps_noinput_REP1_Final.sh
* N2: shuffle_bam_pseudo_selfreps_noinput_REP2_Final.sh

Example run
```bash
qsub iqr_selfreps_noinput.sh Cell_REP1 Cell_REP2 Cell_Merged_Self_IDR
```

<a name="Additional info"></a>
## Additional info
> Authors: Khushbu Patel (@kpatel427), Apexa Modi (@apexamodi), Jo Lynne Rokita (@jharenza)

> Organization: The Children's Hospital of Philadelphia (CHOP)


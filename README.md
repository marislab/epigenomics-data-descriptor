# epigenomics-data-descriptor
Scripts to call super enhancers, filter and generate ChIP-seq heatmaps for epigenomics data descriptor paper

## Table of Contents
- [1. LILY (Super Enhancer calling)](#1.%20LILY%20(Super%20Enhancer%20calling))
- [2. Filtering Super Enhancers (SEs)](#2.%20Filtering%20Super%20Enhancers%20(SEs))
- [3. ChipSeq Heatmaps](3.%20ChipSeq%20Heatmaps)
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

<a name="2. Filtering Super Enhancers (SEs)"></a>
## 2. Filtering Super Enhancers (SEs)
This step filters for super enhances present in two or more MYCN amplified/non-amplified lines.
Annotate all SEs called by LILY before running the filtering script.

### Prerequisites
Annotated SEs called by LILY (Annotations were performed using [Homer](http://homer.ucsd.edu/homer/ngs/annotation.html))

### How to run
Change paths to directories to reflect paths to your files.

```R
Rscript SEfilter.R
```

<a name="3. ChipSeq Heatmaps"></a>
## 3. ChipSeq Heatmaps
ChipSeq heatmaps were generated for MYCN (annotating top 5K peaks) and all histone marks for COGN415 line (annotating filtered SEs from step 2)

### Prerequisites
[deepTools 3.2.0](https://deeptools.readthedocs.io/en/develop/content/installation.html)

### How to run
To generate MYCN heatmaps (with top 5K peaks)
```R
Rscript makeHeatmaps.R
```

To generate COGN415-histone-marks (with filtered SEs)
```bash
computeMatrix reference-point -S bigwigs/COGN415-H3K27Ac.bw bigwigs/COGN415-H3K27me3.bw bigwigs/COGN415-H3K4me1.bw bigwigs/COGN415-H3K4me3.bw -R SE_filtered/SE_filtered.bed -a 4000 -b 4000 --sortUsing max --skipZeros -o COGN415.mat.gz

plotHeatmap -m COGN415.mat.gz --colorList 'white,#ff7400' 'white,#004000' 'white,#00007F' 'white,#6F326F' --whatToShow 'heatmap and colorbar' --sortRegions descend --sortUsing max --zMin 0 --zMax 4 --regionsLabel Super_Enhancers -out COGN415.png

plotProfile -m COGN415.mat.gz -out COGN415.profile.png --perGroup  --colors orange green blue purple --regionsLabel Super_Enhancers --plotHeight 10 --plotWidth 15

```

<a name="Additional info"></a>
## Additional info
> Author: patelk26@email.chop.edu

> Organization: The Children's Hospital of Philadelphia (CHOP)


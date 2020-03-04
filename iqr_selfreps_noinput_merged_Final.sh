#!/bin/bash

#$ -cwd
#$ -l mem_free=25G
#$ -l h_vmem=25G
#$ -S /bin/bash
#$ -sync n
#$ -pe smp 1

#Tutorial/code source can be found on: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html 

# USAGE: sh pseudorep_idr.sh <chip BAM rep1> <chip BAM rep2> <NAME for IDR output>

# This script will take the BAM files and perform the following steps: 
    ## Merge BAMs for ChiP files,
    ## Shuffle reads and split into two new BAM files (pseudo-replicates), 
    ## Merge BAMs for Input files,
    ## Shuffle reads and split into two new BAM files (pseudo-replicates), 
    ## Call peaks on pseudo-replicates with MACS2 , 
    ## Sort peaks called on pseudo-replicates,
    ## IDR analysis using pseudo-replicate peak calls

treatFile1=`basename $1`
treatFile2=`basename $2`
treatFile3=`basename $3`
EXPT=$4

NAME1=`basename $treatFile1 _full.bam`
NAME2=`basename $treatFile2 _full.bam`
NAME3=`basename $treatFile3 _full.bam`

# Make Directories
#mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/idr_chipseq/macs
#mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/pooled_pseudoreps
#mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/idr_chipseq/tmp_${EXPT}

# Set paths
baseDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/aligned
macsDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Self_Rep_IDR_Standardized/macs
outputDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Self_Rep_IDR_Standardized/idr_results
#tmpDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/peakfiles_hg19rev/sorted_peakfiles/idr_chipseq/tmp_${EXPT}

#Peak calling on real replicates

#Peak calling on pseudoreplicates
echo "Calling peaks for replicate1"
macs2 callpeak -t ${baseDir}/${NAME1}.bam -f BAMPE -g hs -n $macsDir/${NAME1} --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME1}_peaks_macs2.log

echo "Calling peaks for replicate2"
macs2 callpeak -t ${baseDir}/${NAME2}.bam -f BAMPE -g hs -n $macsDir/${NAME2} --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME2}_peaks_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${NAME1}_peaks.narrowPeak > $macsDir/${NAME1}_peaks_sorted.narrowPeak
sort -k8,8nr $macsDir/${NAME2}_peaks.narrowPeak > $macsDir/${NAME2}_peaks_sorted.narrowPeak

echo "Calling peaks on merged files."
macs2 callpeak -t ${baseDir}/${NAME3}.merge.sort.bam -f BAMPE -g hs -n $macsDir/${NAME3}_merged --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME3}_peaks_merged_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${NAME3}_merged_peaks.narrowPeak > $macsDir/${NAME3}_merged_peaks_sorted.narrowPeak

#Independent replicate IDR
echo "Running IDR on replicates..."
idr --samples $macsDir/${NAME1}_peaks_sorted.narrowPeak $macsDir/${NAME2}_peaks_sorted.narrowPeak --peak-list $macsDir/${NAME3}_merged_peaks_sorted.narrowPeak --input-file-type narrowPeak --output-file ${EXPT} --rank p.value --plot


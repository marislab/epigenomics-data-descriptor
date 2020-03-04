#!/bin/bash

#$ -cwd
#$ -l mem_free=75G
#$ -l h_vmem=75G
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
EXPT=$3

NAME1=`basename $treatFile1 _full.bam`
NAME2=`basename $treatFile2 _full.bam`

# Make Directories
mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/idr_chipseq/macs
mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/idr_chipseq/pooled_pseudoreps
mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/idr_chipseq/tmp_${EXPT}

# Set paths
baseDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/aligned
macsDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Pseduo_Rep_IDR/idr_chipseq/macs
outputDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Pseduo_Rep_IDR/idr_chipseq/pooled_pseudoreps
tmpDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Pseduo_Rep_IDR/idr_chipseq/tmp_${EXPT}

#Merge treatment BAMS
echo "Merging BAM files for pseudoreplicates..."
samtools merge -u ${tmpDir}/${NAME1}_${NAME2}_merged.bam $baseDir/${treatFile1}.bam $baseDir/${treatFile2}.bam
samtools view -H ${tmpDir}/${NAME1}_${NAME2}_merged.bam > ${tmpDir}/${EXPT}_header.sam

#Split merged treatments
echo "Split merged treatments..."
nlines=$(samtools view ${tmpDir}/${NAME1}_${NAME2}_merged.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number

echo "Shuffle and split into two same files"
samtools view ${tmpDir}/${NAME1}_${NAME2}_merged.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${EXPT}_" 
#This will shuffle the lines in the file and split it into two SAM files

echo "Shuffle completed"
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}_00 | samtools view -bS - > ${outputDir}/${EXPT}_00.bam
cat ${tmpDir}/${EXPT}_header.sam ${tmpDir}/${EXPT}_01 | samtools view -bS - > ${outputDir}/${EXPT}_01.bam

#Merge bams for merged peak calling
echo "Merging pseudorep BAM files for peakcalling..."
samtools merge -u ${outputDir}/${NAME1}_${NAME2}_merged_pseudo.bam ${outputDir}/${EXPT}_00.bam ${outputDir}/${EXPT}_01.bam


#Peak calling on pseudoreplicates
echo "Calling peaks for pseudoreplicate1"
macs2 callpeak -t ${outputDir}/${EXPT}_00.bam -f BAMPE -g hs -n $macsDir/${NAME1}_pr --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME1}_pr_macs2.log

echo "Calling peaks for pseudoreplicate2"
macs2 callpeak -t ${outputDir}/${EXPT}_01.bam -f BAMPE -g hs -n $macsDir/${NAME2}_pr --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME2}_pr_macs2.log

echo "Calling peaks for pseudoreplicate merged"
macs2 callpeak -t ${outputDir}/${NAME1}_${NAME2}_merged_pseudo.bam -f BAMPE -g hs -n $macsDir/${NAME1}_${NAME2}_pr_merged --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME1}_${NAME2}_pr_merged_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${NAME1}_pr_peaks.narrowPeak > $macsDir/${NAME1}_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${NAME2}_pr_peaks.narrowPeak > $macsDir/${NAME2}_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${NAME1}_${NAME2}_pr_merged_peaks.narrowPeak > $macsDir/${NAME1}_${NAME2}_pr_merged_sorted.narrowPeak

#Independent replicate IDR
echo "Running IDR on pseudoreplicates..."
idr --samples $macsDir/${NAME1}_pr_sorted.narrowPeak $macsDir/${NAME2}_pr_sorted.narrowPeak --peak-list $macsDir/${NAME1}_${NAME2}_pr_merged_sorted.narrowPeak --input-file-type narrowPeak --output-file ${EXPT}_pseudorep-idr-merged --rank p.value --plot

# Remove the tmp directory
rm -r $tmpDir

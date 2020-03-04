#!/bin/bash

#$ -cwd
#$ -l mem_free=60G
#$ -l h_vmem=60G
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
#treatFile2=`basename $2`
EXPT=$2

NAME2=`basename $treatFile1 _full.bam`
#NAME2=`basename $treatFile2 _full.bam`

# Make Directories
mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Pseduo_Self_Rep2_IDR/macs
mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Pseduo_Self_Rep2_IDR/pooled_pseudoreps
mkdir -p /mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Pseduo_Self_Rep_IDR/tmp_${EXPT}

# Set paths
baseDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/aligned
#baseDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/peakfiles_hg19rev/sorted_peakfiles
macsDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Pseduo_Self_Rep2_IDR/macs
outputDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/revision/IDR/Pseduo_Self_Rep2_IDR/pooled_pseudoreps
#tmpDir=/mnt/isilon/maris_lab/target_nbl_ngs/AtacSeqAnalysis/peakfiles_hg19rev/sorted_peakfiles/Pseduo_Self_Rep_IDR/tmp_${EXPT}

#Merge treatment BAMS
echo "Merging BAM files for pseudoreplicates..."
samtools merge -u ${tmpDir}/${NAME2}_${NAME2}_merged.bam $baseDir/${treatFile1}.bam $baseDir/${treatFile2}.bam
samtools view -H ${baseDir}/${NAME2}.bam > ${tmpDir}/${EXPT}_REP2_header.sam
samtools view -H ${baseDir}/${NAME2}.bam > ${tmpDir}/${EXPT}_REP2_header.sam

#Split merged treatments
echo "Split each replicate..."
nlines=$(samtools view ${baseDir}/${NAME2}.bam | wc -l ) # Number of reads in the BAM file
nlines=$(( (nlines + 1) / 2 )) # half that number
nlines2=$(samtools view ${baseDir}/${NAME2}.bam | wc -l )
nlines2=$(( (nlines2 + 1) / 2 )) 

echo "Shuffle and split into two same files (per replicate)"
samtools view ${baseDir}/${NAME2}.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${EXPT}_REP2_"
samtools view ${baseDir}/${NAME2}.bam | shuf - | split -d -l ${nlines} - "${tmpDir}/${EXPT}_REP2_" 
#This will shuffle the lines in the file and split it into two SAM files
samtools view ${baseDir}/${NAME2}.bam | shuf - | split -d -l ${nlines2} - "${tmpDir}/${EXPT}_REP2_"
samtools view ${baseDir}/${NAME2}.bam | shuf - | split -d -l ${nlines2} - "${tmpDir}/${EXPT}_REP2_"


echo "Shuffle completed"
cat ${tmpDir}/${EXPT}_REP2_header.sam ${tmpDir}/${EXPT}_REP2_00 | samtools view -bS - > ${outputDir}/${EXPT}_REP2_00.bam
cat ${tmpDir}/${EXPT}_REP2_header.sam ${tmpDir}/${EXPT}_REP2_01 | samtools view -bS - > ${outputDir}/${EXPT}_REP2_01.bam
cat ${tmpDir}/${EXPT}_REP2_header.sam ${tmpDir}/${EXPT}_REP2_00 | samtools view -bS - > ${outputDir}/${EXPT}_REP2_00.bam
cat ${tmpDir}/${EXPT}_REP2_header.sam ${tmpDir}/${EXPT}_REP2_01 | samtools view -bS - > ${outputDir}/${EXPT}_REP2_01.bam

echo "Merged self pseudo rep bams"
samtools merge -u ${outputDir}/${NAME2}_merged_pseudo.bam ${outputDir}/${NAME2}_Self_REP2_00.bam ${outputDir}/${NAME2}_Self_REP2_01.bam

#Peak calling on pseudoreplicates
echo "Calling peaks for pseudoreplicate1"
macs2 callpeak -t ${outputDir}/${NAME2}_Self_REP2_00.bam -f BAMPE -g hs -n $macsDir/${NAME2}_00_pr --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME2}_00_pr_macs2.log

#macs2 callpeak -t ${outputDir}/${EXPT}_REP2_00.bam -f BAM -g hs -n $macsDir/${NAME2}_00_pr -B -p 1e-3  2> $macsDir/${NAME2}_00_pr_macs2.log

echo "Calling peaks for pseudoreplicate2"
macs2 callpeak -t ${outputDir}/${NAME2}_Self_REP2_01.bam -f BAMPE -g hs -n $macsDir/${NAME2}_01_pr --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME2}_01_pr_macs2.log

#macs2 callpeak -t ${outputDir}/${EXPT}_REP2_01.bam -f BAM -g hs -n $macsDir/${NAME2}_01_pr -B -p 1e-3  2> $macsDir/${NAME2}_01_pr_macs2.log

echo "Calling peaks for merged pseuedorep"
macs2 callpeak -t ${outputDir}/${NAME2}_merged_pseudo.bam -f BAMPE -g hs -n $macsDir/${NAME2}_merged_self_pr --nomodel -B -p 1e-3 --verbose 3 --extsize 200 --shift -100 --SPMR > $macsDir/${NAME2}_pr_merged_macs2.log

#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/${NAME2}_00_pr_peaks.narrowPeak > $macsDir/${NAME2}_00_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${NAME2}_01_pr_peaks.narrowPeak > $macsDir/${NAME2}_01_pr_sorted.narrowPeak
sort -k8,8nr $macsDir/${NAME2}_01_pr_peaks.narrowPeak > $macsDir/${NAME2}_merged_self_pr.narrowPeak

#sort -k8,8nr $macsDir/${NAME2}_00_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME2}_00_pr_sorted.narrowPeak
#sort -k8,8nr $macsDir/${NAME2}_01_pr_peaks.narrowPeak | head -n 100000 > $macsDir/${NAME2}_01_pr_sorted.narrowPeak

#Independent replicate IDR
echo "Running IDR on pseudoreplicates..."
idr --samples $macsDir/${NAME2}_00_pr_sorted.narrowPeak $macsDir/${NAME2}_01_pr_sorted.narrowPeak --peak-list $macsDir/${NAME2}_merged_self_pr.narrowPeak --input-file-type narrowPeak --output-file ${EXPT}_REP2_self_pseudo-idr --rank p.value --plot
#idr --samples $macsDir/${NAME2}_00_pr_sorted.narrowPeak $macsDir/${NAME2}_01_pr_sorted.narrowPeak --input-file-type narrowPeak --output-file ${EXPT}_REP2_self_pseudo-idr --rank p.value --plot

# Remove the tmp directory
rm -r $tmpDir

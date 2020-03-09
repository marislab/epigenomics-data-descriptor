# script to generate heatmaps for MYCN and CMYC using new MYCN top 5K peaks
# setwd("/Volumes/target_nbl_ngs/KP/chipseq_heatmaps/2020-02-26")
# with new MYCN_top5K_peaks, n = 1335

mycn_bw <- '/mnt/isilon/maris_lab/target_nbl_ngs/ChipSeqAnalysis/mycn-myc-full-dataset/peakfiles/bigwig/'
mycn_alt_bw <- '/mnt/isilon/maris_lab/target_nbl_ngs/ChipSeqAnalysis/mycn-chip-testdir/peakfiles/bigwig/'
cl <- c('COGN415-NMYC-20171205','Kelly-MYCN-20150914', 'LAN5-NMYC-20171205', 'NB1643-NMYC-20171205', 'NB69-NMYC-20171205', 'NGP-MYCN-20150914')

top_5K_mycn_peaks <- '/mnt/isilon/maris_lab/target_nbl_ngs/KP/chipseq_heatmaps/2020-02-26/files/mycnAmp_top5K_peaks.bed'
MYCN_heatmap_dir <- '/mnt/isilon/maris_lab/target_nbl_ngs/KP/chipseq_heatmaps/2020-02-26/MYCN/'

#.. For MYCN-----

for(x in cl){
  if(x == 'Kelly-MYCN-20150914' | x == 'NGP-MYCN-20150914') {

    print(paste0('Processing...', x))
    # compute matrix
    system(paste0('computeMatrix reference-point -S ',mycn_alt_bw,x,'.macs2.SPMR.IP.sorted.bw -R ',top_5K_mycn_peaks,' -a 4000 -b 4000 -o ', MYCN_heatmap_dir, x, '.mat.gz'))

    # make heatmaps
    system(paste0('plotHeatmap -m ', MYCN_heatmap_dir, x, '.mat.gz --colorList "white,red" --heatmapHeight 16 -out ', MYCN_heatmap_dir,x,'.png --sortUsing max --sortRegions descend --zMin 0 --zMax 4')     )

    # make profiles
    system(paste0('plotProfile -m ', MYCN_heatmap_dir, x,'.mat.gz -out ', MYCN_heatmap_dir,x,'.profile.png --colors red --regionsLabel MYCN --plotHeight 10 --plotWidth 15'))

    } else {
    print(paste0('Processing...', x))
    # compute matrix
    system(paste0('computeMatrix reference-point -S ',mycn_bw,x,'.macs2.SPMR.IP.sorted.bw -R ',top_5K_mycn_peaks,' -a 4000 -b 4000 -o ', MYCN_heatmap_dir, x, '.mat.gz'))

    # make heatmaps
    system(paste0('plotHeatmap -m ', MYCN_heatmap_dir, x, '.mat.gz --colorList "white,red" --heatmapHeight 16 -out ', MYCN_heatmap_dir,x,'.png --sortUsing max --sortRegions descend --zMin 0 --zMax 4')     )

    # make profiles
    system(paste0('plotProfile -m ', MYCN_heatmap_dir, x,'.mat.gz -out ', MYCN_heatmap_dir,x,'.profile.png --colors red --regionsLabel MYCN --plotHeight 10 --plotWidth 15'))
  }
}


#.. For CMYC-----
cmyc_bw <- '/mnt/isilon/maris_lab/target_nbl_ngs/ChipSeqAnalysis/mycn-myc-full-dataset/peakfiles/bigwig/'
cmyc_cl <- c('KELLY-CMYC-20171205','SKNAS-CMYC-20171205','NB69-CMYC-20171205')
CMYC_heatmap_dir <- '/mnt/isilon/maris_lab/target_nbl_ngs/KP/chipseq_heatmaps/2020-02-26/CMYC/'

for(x in cmyc_cl){
  print(paste0('Processing...', x))
  # compute matrix
  system(paste0('computeMatrix reference-point -S ',cmyc_bw,x,'.macs2.SPMR.IP.sorted.bw -R ',top_5K_mycn_peaks,' -a 4000 -b 4000 -o ', CMYC_heatmap_dir, x, '.mat.gz'))
  
  # make heatmaps
  system(paste0('plotHeatmap -m ', CMYC_heatmap_dir, x, '.mat.gz --colorList "white,blue" --heatmapHeight 16 -out ', CMYC_heatmap_dir,x,'.png --sortUsing max --sortRegions descend --zMin 0 --zMax 4')     )
  
  # make profiles
  system(paste0('plotProfile -m ', CMYC_heatmap_dir, x,'.mat.gz -out ', CMYC_heatmap_dir,x,'.profile.png --colors blue --regionsLabel MYC --plotHeight 10 --plotWidth 15'))
}





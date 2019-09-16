# script to generate heatmaps
# COGN415 all histone marks with SE_filtered list
# top 5K peaks for MYCN amp/non-amp

mycn_bw <- '/mnt/isilon/maris_lab/target_nbl_ngs/ChipSeqAnalysis/mycn-myc-full-dataset/peakfiles/bigwig/'
cl <- c('COGN415-NMYC-20171205','Kelly-MYCN-20150914', 'LAN5-NMYC-20171205', 'NB1643-NMYC-20171205', 'NB69-NMYC-20171205', 'NGP-MYCN-20150914')

top_5K_peaks_dir <- '~/KP/heatmaps/2019-09-11/MYCN/top_5K_peaks/'
heatmap_dir <- '~/KP/heatmaps/2019-09-11/MYCN/'

#.. For MYCN-----

for(x in cl){
  print(paste0('Processing...', x))
  # compute matrix
  system(paste0('computeMatrix reference-point -S ',mycn_bw,x,'.macs2.SPMR.IP.sorted.bw -R ',top_5K_peaks_dir,x,'.macs2.SPMR.filtered.narrowPeak.top5K.peaks.bed -a 4000 -b 4000 -o ', heatmap_dir, x, '.mat.gz'))
  
  # make heatmaps
  system(paste0('plotHeatmap -m ', heatmap_dir, x, '.mat.gz --colorList "white,red" --heatmapHeight 25 --heatmapWidth 3 -out ',heatmap_dir,x,'.png --sortUsing max --sortRegions descend --zMin 0 --zMax 4'))
  
  # make profiles
  system(paste0('plotProfile -m ', heatmap_dir, x,'.mat.gz -out ',heatmap_dir,x,'.profile.png --colors red --regionsLabel Super_Enhancers --plotHeight 10 --plotWidth 15'))
}


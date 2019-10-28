# for Enh
cl <-  c('LAN5', 'NB69')
for(i in cl) {
  
  if(i == 'LAN5'){
    
    print("Computing matrix for LAN5...")
    system(paste0("/home/patelk26/miniconda3/bin/computeMatrix reference-point -S bigwigs/",i,"-MYCN.bw bigwigs/",i,"-H3K27Ac.bw bigwigs/",i,"-H3K27me3.bw bigwigs/",i,"-H3K4me1.bw bigwigs/",i,"-H3K4me3.bw bigwigs/",i,"-ATAC.bw -R files/",i,"_Enh.bed  -a 10000 -b 10000 --sortUsing max -o ",Sys.Date(),'_',i,"-Enh.mat.gz --outFileSortedRegions ", Sys.Date(),'_',i,"-Enh_regionsSorted.bed --outFileNameMatrix ",Sys.Date(),'_',i,"_Enh.matrix.tab"))
    
    
    # +/- 10kb for Enh only
    # make heatmap
    print(paste0("Generating heatmap for ",i))
    system(paste0("/home/patelk26/miniconda3/bin/plotHeatmap -m 2019-10-22_",i,"-Enh.mat.gz --colorList 'white,red' 'white,#ff7400' 'white,#004000' 'white,#00007F' 'white,#6F326F' 'white,#fc03d3' --sortRegions descend --sortUsing max --zMin 0 --zMax 10 --regionsLabel Enhancer --heatmapHeight 16 -out ",Sys.Date(),'_',i,"_Enh.pdf"))
    # make profile
    print(paste0("Generating profile plot for ",i))
    system(paste0("/home/patelk26/miniconda3/bin/plotProfile -m 2019-10-22_",i,"-Enh.mat.gz -out ", Sys.Date(),'_',i,"_Enh.profile.pdf --perGroup  --colors red orange green blue purple magenta --regionsLabel Enhancer --plotHeight 10 --plotWidth 15"))
    
    
  }
  if(i == 'NB69') {
    
    print("Computing matrix for NB69...")
    system(paste0("/home/patelk26/miniconda3/bin/computeMatrix reference-point -S bigwigs/",i,"-MYCN.bw bigwigs/",i,"-CMYC.bw bigwigs/",i,"-H3K27Ac.bw bigwigs/",i,"-H3K27me3.bw bigwigs/",i,"-H3K4me1.bw bigwigs/",i,"-H3K4me3.bw bigwigs/",i,"-ATAC.bw -R files/",i,"_Enh.bed  -a 10000 -b 10000 --sortUsing max -o ",Sys.Date(),'_',i,"-Enh.mat.gz --outFileSortedRegions ", Sys.Date(),'_',i,"-Enh_regionsSorted.bed --outFileNameMatrix ",Sys.Date(),'_',i,"_Enh.matrix.tab"))
    
    # +/- 10kb
    # make heatmap
    print(paste0("Generating heatmap for ",i))
    system(paste0("/home/patelk26/miniconda3/bin/plotHeatmap -m 2019-10-22_",i,"-Enh.mat.gz --colorList 'white,red' 'white,#04baf0' 'white,#ff7400' 'white,#004000' 'white,#00007F' 'white,#6F326F' 'white,#fc03d3' --sortRegions descend --sortUsing max --zMin 0 --zMax 10 --regionsLabel Enhancer --heatmapHeight 16 -out ",Sys.Date(),'_',i,"_Enh.pdf"))
    # make profile
    print(paste0("Generating profile plot for ",i))
    system(paste0("/home/patelk26/miniconda3/bin/plotProfile -m 2019-10-22_",i,"-Enh.mat.gz -out ",Sys.Date(),'_',i,"_Enh.profile.pdf --perGroup  --colors red lightblue orange green blue purple magenta --regionsLabel Enhancer --plotHeight 10 --plotWidth 15"))
    
  }
  
}

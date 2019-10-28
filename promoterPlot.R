cl <-  c('LAN5', 'NB69')

# promoters
for(i in cl){
  print(paste0("Plotting heatmaps for promoters: ",i))
  
  if(i == 'LAN5'){
    print("Computing matrix for LAN5...")
    system(paste0("/home/patelk26/miniconda3/bin/computeMatrix reference-point -S bigwigs/",i,"-MYCN.bw  bigwigs/",i,"-H3K27Ac.bw bigwigs/",i,"-H3K27me3.bw bigwigs/",i,"-H3K4me1.bw bigwigs/",i,"-H3K4me3.bw bigwigs/",i,"-ATAC.bw -R files/",i,"_promoter.bed  -a 4000 -b 4000 --sortRegions keep -o ",Sys.Date(),'_',i,"-promoter.mat.gz --outFileSortedRegions ", Sys.Date(),'_',i,"-promoter_regionsUnSorted.bed --outFileNameMatrix ",Sys.Date(),'_',i,"-promoter.unsorted.matrix.tab"))
    
    # # make heatmap
    # print(paste0("Generating heatmap for ",i))
    # system(paste0("/home/patelk26/miniconda3/bin/plotHeatmap -m 2019-10-22_",i,"-promoter.mat.gz --colorList 'white,red' 'white,#ff7400' 'white,#004000' 'white,#00007F' 'white,#6F326F' 'white,#fc03d3' --sortRegions descend --sortUsing max --zMin 0 --zMax 4 --regionsLabel Promoter --heatmapHeight 16 --missingDataColor white -out ", Sys.Date(),'_',i,"-promoter.pdf"))
    # # make profile
    # print(paste0("Generating profile plot for ",i))
    # system(paste0("/home/patelk26/miniconda3/bin/plotProfile -m 2019-10-22_",i,"-promoter.mat.gz -out ",Sys.Date(),'_',i,"-promoter.profile.pdf --perGroup  --colors red orange green blue purple magenta --regionsLabel Promoter --plotHeight 10 --plotWidth 15"))
    # 
    }
  if(i == 'NB69') {
    
    print("Computing matrix for NB69...")
    system(paste0("/home/patelk26/miniconda3/bin/computeMatrix reference-point -S bigwigs/",i,"-MYCN.bw bigwigs/",i,"-CMYC.bw bigwigs/",i,"-H3K27Ac.bw bigwigs/",i,"-H3K27me3.bw bigwigs/",i,"-H3K4me1.bw bigwigs/",i,"-H3K4me3.bw bigwigs/",i,"-ATAC.bw -R files/",i,"_promoter.bed --sortRegions keep  -a 4000 -b 4000 --sortUsing max -o ",Sys.Date(),'_',i,"-promoter.mat.gz --outFileSortedRegions ", Sys.Date(),'_',i,"-promoter_regionsUnSorted.bed --outFileNameMatrix ",Sys.Date(),'_',i,"-promoter.unsorted.matrix.tab"))
    
    # # make heatmap
    # print(paste0("Generating heatmap for ",i))
    # system(paste0("/home/patelk26/miniconda3/bin/plotHeatmap -m 2019-10-22_",i,"-promoter.mat.gz --colorList 'white,red' 'white,#04baf0' 'white,#ff7400' 'white,#004000' 'white,#00007F' 'white,#6F326F' 'white,#fc03d3' --sortRegions descend --sortUsing max --zMin 0 --zMax 4 --regionsLabel Promoter --heatmapHeight 16 --missingDataColor white -out ", Sys.Date(),'_',i,"-promoter.pdf"))
    # # make profile
    # print(paste0("Generating profile plot for ",i))
    # system(paste0("/home/patelk26/miniconda3/bin/plotProfile -m 2019-10-22_",i,"-promoter.mat.gz -out ", Sys.Date(),'_',i,"-promoter.profile.pdf --perGroup  --colors red lightblue orange green blue purple magenta --regionsLabel Promoter --plotHeight 10 --plotWidth 15"))
    # 
    }
  
}
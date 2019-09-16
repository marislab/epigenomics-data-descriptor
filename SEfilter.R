# script to keep LILY's SE present in 2+ lines for MYCN amp/non-amp
library(data.table)
library(tidyverse)

amp_cl <- c('COGN415-H3K27Ac', 'KELLY-H3K27Ac', 'LAN5-H3K27Ac', 'NB1643-H3K27Ac','NGP-H3K27Ac', 'SKNBE2C-H3K27Ac')
nonAmp_cl <- c('NB69-H3K27Ac', 'NBLS-H3K27Ac', 'SKNAS-H3K27Ac', 'SKNFI-H3K27Ac')

lily_dir <- '~/KP/enhancer-rank-list-H3K27Ac/LILY/result/SE_scores_anno/'
out_dir <-  '~/KP/enhancer-rank-list-H3K27Ac/LILY/result/SE_scores_anno/SE_filtered/'

#..MYCN AMPLIFIED LINES ----
for(x in amp_cl){
  lily_file <- paste0(x,'_lily')
  print(lily_file)
  # reading the file
  assign(lily_file, read.delim(paste0(lily_dir,x, '.scores.bed.SE.scores.bed.sorted.SE.bed.anno.txt'), header = T))
  # data subsetting
  name <-  paste0(lily_file,'.subset')
  assign(name, as.data.frame(get(lily_file)[,c(1,16)]))
  # renaming columns
  setnames(get(name), names(get(name))[1], paste0(x,'_PeakID'))
  setnames(get(name), names(get(name))[2], paste0('Gene.Name'))
}

# merging all dataframes
merge.df <- Reduce(function(x, y) merge(x, y, by = c('Gene.Name'), all = TRUE), list(`COGN415-H3K27Ac_lily.subset`, `KELLY-H3K27Ac_lily.subset`, `LAN5-H3K27Ac_lily.subset`,
                                                                   `NB1643-H3K27Ac_lily.subset`, `NGP-H3K27Ac_lily.subset`, `SKNBE2C-H3K27Ac_lily.subset`))

# filter and keep all SE's present in 2+ MYCN amp lines (having peakIDs)
merge.df$na_count <- apply(merge.df, 1, function(x) sum(is.na(x)))

# filter to keep rows with na_count 4 or less
mycn_amp_SE <- setDT(merge.df)[merge.df$na_count <= 4, 'Gene.Name']
mycn_amp_SE <-  unique(mycn_amp_SE$Gene.Name)

write.table(mycn_amp_SE, file = paste0(out_dir,Sys.Date(),'_mycnAmp_SE_LILY_calls_filtered.txt'), quote = F, col.names = F, row.names = F, sep = '\t')




#..MYCN NON-AMPLIFIED LINES ----
for(x in nonAmp_cl){
  lily_file <- paste0(x,'_lily')
  print(lily_file)
  # reading the file
  assign(lily_file, read.delim(paste0(lily_dir,x, '.scores.bed.SE.scores.bed.sorted.SE.bed.anno.txt'), header = T))
  # data subsetting
  name <-  paste0(lily_file,'.subset')
  assign(name, as.data.frame(get(lily_file)[,c(1,16)]))
  # renaming columns
  setnames(get(name), names(get(name))[1], paste0(x,'_PeakID'))
  setnames(get(name), names(get(name))[2], paste0('Gene.Name'))
}

# merging all dataframes
merge.df1 <- Reduce(function(x, y) merge(x, y, by = 'Gene.Name', all = TRUE), list(`NB69-H3K27Ac_lily.subset`, `NBLS-H3K27Ac_lily.subset`, `SKNAS-H3K27Ac_lily.subset`,
                                                                                   `SKNFI-H3K27Ac_lily.subset`))

# filter and keep all SE's present in 2+ MYCN Non-amp lines (having peakIDs)
merge.df1$na_count <- apply(merge.df1, 1, function(x) sum(is.na(x)))

# filter to keep rows with na_count 2 or less
mycn_Nonamp_SE <- setDT(merge.df1)[merge.df1$na_count <= 2, 'Gene.Name']
mycn_Nonamp_SE <-  unique(mycn_Nonamp_SE$Gene.Name)

write.table(mycn_Nonamp_SE, file = paste0(out_dir,Sys.Date(),'_mycnNonAmp_SE_LILY_calls_filtered.txt'), quote = F, col.names = F, row.names = F, sep = '\t')

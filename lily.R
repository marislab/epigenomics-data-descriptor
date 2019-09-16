# Rscript to run LILY on all H3K27Ac cell lines (amp/non-amp)

cl <- c('COGN415-H3K27Ac', 'KELLY-H3K27Ac', 'LAN5-H3K27Ac', 'NB1643-H3K27Ac',
        'NB69-H3K27Ac', 'NBLS-H3K27Ac', 'NGP-H3K27Ac', 'SKNAS-H3K27Ac', 'SKNBE2C-H3K27Ac', 'SKNFI-H3K27Ac',
        'SKNSH-H3K27Ac')

dir <- '~/KP/enhancer-rank-list-H3K27Ac/LILY/data/' 
result.dir <-  '~/KP/enhancer-rank-list-H3K27Ac/LILY/result/'

for(x in cl){
  system(paste0('cat runLILY.R | R --slave --args ', dir, x, ' ',result.dir ,' 12500 2500 hg19_refseq.ucsc hg19.chrom.sizes'))
}

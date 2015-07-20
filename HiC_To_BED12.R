library(data.table)
dat = fread('GSE65126/GSM1587701_HiC_mouse_liver_rep1_50000.txt')
# dat = as.data.table(dat)
# setDT(dat)[, id:= .GRP, by = c('chrom1','start1','end1')]

########################################## transform HiC to BED12 ########################################## 
# transformHiCToBED12 <- function(dat.t)
# {
#   dat.t$expected_count = NULL
#   dat.t$observed_count = NULL
#   dat.t$end2 = NULL
#   dat.t$chrom2 = NULL
#   # blockSizes = sizes rep(50000,n)
#   dat.t$blockSizes = 50000 
#   tmp = data.table::data.table(dat.t)[,lapply(.SD, paste, collapse = ","),c('chrom1','start1','end1')]
#   # blockCount = n = number of blocks
#   tmp$blockCount = data.table(dat.t)[,lapply(.SD, length),c('chrom1','start1','end1')][[4]]
#   # blockStarts = A comma-separated list of block starts
#   setnames(tmp, c("start2","chrom1","start1","end1"), c("blockStarts","chrom","chromStart","chromEnd"))
#   tmp$name = '.'
#   tmp$score = 0
#   tmp$strand = '.'
#   tmp$thickStart = '.'
#   tmp$thickEnd = '.'
#   tmp$itemRgb = '.'
#   tmp$chrom = paste('chr',tmp$chrom,sep = '')
#   setcolorder(tmp, c("chrom", "chromStart", "chromEnd","name","score","strand",
#                      "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"))
#   return(tmp)
# }
########################################## transform HiC to BED12 ########################################## 

########################################## transform HiC to BED12 ########################################## 
transformHiCToBED12 <- function(tmp)
{
  tmp$start.1 = ifelse(test = tmp$start2 < tmp$start1, tmp$start2, tmp$start1)
  tmp$end.1 = ifelse(test = tmp$end2 < tmp$end1 , tmp$end2, tmp$end1)
  tmp$start.2 = ifelse(test = tmp$start2 > tmp$start1, tmp$start2, tmp$start1)
  tmp$end.2 = ifelse(test = tmp$end2 > tmp$end1 , tmp$end2, tmp$end1)
  tmp$blockStarts = paste(tmp$start.1-tmp$start.1,tmp$start.2-tmp$start.1,sep = ',')
  tmp$blockSizes = paste(tmp$end.1-tmp$start.1,tmp$end.2-tmp$start.2,sep=',')
  tmp$blockCount = 2
  tmp$chrom = paste('chr',tmp$chrom1,sep='')
  tmp$chromStart = tmp$start.1
  tmp$chromEnd = tmp$end.2
  tmp$score = 0
  tmp$strand = '.'
  tmp$thickStart = 0
  tmp$thickEnd = 0
  tmp$chrom1 = NULL
  tmp$start1 = NULL
  tmp$end1 = NULL
  tmp$expected_count = NULL
  tmp$observed_count= NULL
  tmp$chrom2 = NULL
  tmp$start2 = NULL
  tmp$end2 = NULL
  tmp$start.1 = NULL
  tmp$start.2 = NULL
  tmp$end.1 = NULL
  tmp$end.2 = NULL
  tmp$itemRgb = '255,0,0'
  tmp = unique(tmp)
  setDT(tmp)[, id:= .GRP, by = c('chrom','chromStart','chromEnd')]
  tmp$name = paste('item',tmp$id,sep='')
  tmp$id = NULL
  setcolorder(tmp, c("chrom", "chromStart", "chromEnd","name","score","strand",
                     "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"))
  return(tmp)
}
########################################## transform HiC to BED12 ########################################## 

# call function to transform HiC data to BED12 format
res = transformHiCToBED12(tmp = dat)
write.table(res,'results.bed',quote=F,row.names=F,col.names = F,sep='\t')

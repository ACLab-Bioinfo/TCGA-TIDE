suppressPackageStartupMessages({
  library(scater)
  library(tidyverse)
  library(stringr)
})

### readDATA
hcc_fpkm <- read.table(gzfile("./data/TCGA-LIHC.htseq_fpkm.tsv.gz"), header = TRUE, row.names = 1, sep = "\t")
#hcc_fpkm <- read.table("./data/TCGA-LIHC.htseq_fpkm.tsv", header = TRUE, row.names = 1, sep = "\t")
map <- read.table("./data/gencode.v22.annotation.gene.probeMap", header=T, row.names=1)
map <- map[rownames(hcc_fpkm),]

### Annotation
rownames(hcc_fpkm) <- uniquifyFeatureNames(rownames(map),map$gene)
### Extract only tumor samples
hcc_T_fpkm <- hcc_fpkm[,sapply(colnames(hcc_fpkm), function(x) unlist(strsplit(x, split="\\."))[4]) == "01A"]
                               
GenerateTIDEinput <- function(gene){
  target.gene <- gene
  dir.create(paste0("./input_TIDE/", target.gene))
  dir.create(paste0("./output_TIDE/", target.gene))
  
  order_data_TUMOR <- hcc_T_fpkm[,order(hcc_T_fpkm[target.gene,],decreasing = TRUE)]
  order_data_TUMOR.normalized <- sweep(order_data_TUMOR, 1, rowMeans(order_data_TUMOR))

  order_data.normalized <- order_data_TUMOR.normalized

  first_last_15 <- order_data.normalized[,c(1:(ncol(order_data.normalized)*15/100),round(ncol(order_data.normalized)-(ncol(order_data.normalized)*15/100)+1):ncol(order_data.normalized))]
  first_last_25 <- order_data.normalized[,c(1:(ncol(order_data.normalized)*25/100),round(ncol(order_data.normalized)-(ncol(order_data.normalized)*25/100)+1):ncol(order_data.normalized))]
  first_last_30 <- order_data.normalized[,c(1:(ncol(order_data.normalized)*30/100),round(ncol(order_data.normalized)-(ncol(order_data.normalized)*30/100)+2):ncol(order_data.normalized))]

  write.table(first_last_15, paste0('./input_TIDE/',target.gene,'/order_data_TUMOR.first_last_15_normalized.',target.gene,'.txt'), sep="\t", quote=F)
  write.table(first_last_25, paste0('./input_TIDE/',target.gene,'/order_data_TUMOR.first_last_25_normalized.',target.gene,'.txt'), sep="\t", quote=F)
  write.table(first_last_30, paste0('./input_TIDE/',target.gene,'/order_data_TUMOR.first_last_30_normalized.',target.gene,'txt'), sep="\t", quote=F)
}

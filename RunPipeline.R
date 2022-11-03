library(ggplot2)
source("./scripts/GenerateTIDEInput.R")
source("./scripts/DrawTIDE.R")
source("./scripts/ModifyTIDEresult.R")

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

gene_list <- read.table("./data/genelist.txt")[,1]

for(gene in gene_list){
  GenerateTIDEinput(gene)
  system(paste0("bash ./scripts/RunTide.sh ", gene))
  ModifyTIDEresult(gene)
  DrawTIDE(gene)
}

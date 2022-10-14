library(ggplot2)
source("./scripts/GenerateTIDEInput.R")
source("./scripts/DrawTIDE.R")
source("./scripts/ModifyTIDEresult.R")

gene_list <- unlist(read.table("./data/genelist.txt", quote = F))

for(gene in gene_list){
  GenerateTIDEinput(gene)
  system("bash ./scripts/RunTide.sh")
  ModifyTIDEresult(gene)
  DrawTIDE(gene)
}


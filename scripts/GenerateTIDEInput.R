GenerateTIDEinput <- function(hcc_T_fpkm, gene){
  target.gene <- gene
  hcc_T_fpkm <- hcc_T_fpkm
  dir.create(paste0("./input_TIDE/", target.gene), recursive=TRUE)
  dir.create(paste0("./output_TIDE/", target.gene), recursive=TRUE)

  order_data_TUMOR <- hcc_T_fpkm[,order(hcc_T_fpkm[target.gene,],decreasing = TRUE)]
  order_data_TUMOR.normalized <- sweep(order_data_TUMOR, 1, rowMeans(order_data_TUMOR))

  order_data.normalized <- order_data_TUMOR.normalized

  first_last_15 <- order_data.normalized[,c(1:(ncol(order_data.normalized)*15/100),round(ncol(order_data.normalized)-(ncol(order_data.normalized)*15/100)+1):ncol(order_data.normalized))]
  first_last_25 <- order_data.normalized[,c(1:(ncol(order_data.normalized)*25/100),round(ncol(order_data.normalized)-(ncol(order_data.normalized)*25/100)+1):ncol(order_data.normalized))]
  first_last_30 <- order_data.normalized[,c(1:(ncol(order_data.normalized)*30/100),round(ncol(order_data.normalized)-(ncol(order_data.normalized)*30/100)+2):ncol(order_data.normalized))]

  write.table(first_last_15, paste0('./input_TIDE/',target.gene,'/order_data_TUMOR.first_last_15_normalized.',target.gene,'.txt'), sep="\t", quote=F)
  write.table(first_last_25, paste0('./input_TIDE/',target.gene,'/order_data_TUMOR.first_last_25_normalized.',target.gene,'.txt'), sep="\t", quote=F)
  write.table(first_last_30, paste0('./input_TIDE/',target.gene,'/order_data_TUMOR.first_last_30_normalized.',target.gene,'.txt'), sep="\t", quote=F)
}

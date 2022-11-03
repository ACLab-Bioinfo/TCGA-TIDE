suppressPackageStartupMessages({
	library(tidyverse)
})
ModifyTIDEresult <- function(gene){

  target.gene <- gene

  order_data_TUMOR <- hcc_T_fpkm[,order(hcc_T_fpkm[target.gene,],decreasing = TRUE)]
  order_data_TUMOR.normalized <- sweep(order_data_TUMOR, 1, rowMeans(order_data_TUMOR))

  first_last_15_TIDE.result <- read.table(paste0("./output_TIDE/", gene, "/TIDE.output.first_last_15_normalized.", gene, ".txt"),  sep="\t",  header = TRUE, row.names = 1)
  first_last_25_TIDE.result <- read.table(paste0("./output_TIDE/", gene, "/TIDE.output.first_last_25_normalized.", gene, ".txt"),  sep="\t",  header = TRUE, row.names = 1)
  first_last_30_TIDE.result <- read.table(paste0("./output_TIDE/", gene, "/TIDE.output.first_last_30_normalized.", gene, ".txt"),  sep="\t",  header = TRUE, row.names = 1)

  first_last_15_TIDE.result = first_last_15_TIDE.result %>%
    add_column(gene.expr = t(order_data_TUMOR.normalized[target.gene,rownames(first_last_15_TIDE.result)]), .after = 1)
  first_last_15_TIDE.result = first_last_15_TIDE.result %>%
    add_column(expr_high_or_low = ifelse(first_last_15_TIDE.result$gene.expr > 0, "high", "low"), .after = 1)
  write.table(first_last_15_TIDE.result, paste0("./output_TIDE/", gene, "/Modified.TIDE.result.first_last_15.with.",target.gene,".expr.txt" ), sep = ",", row.names = FALSE )

  first_last_25_TIDE.result = first_last_25_TIDE.result %>%
    add_column(gene.expr = t(order_data_TUMOR.normalized[target.gene,rownames(first_last_25_TIDE.result)]), .after = 1)
  first_last_25_TIDE.result = first_last_25_TIDE.result %>%
    add_column(expr_high_or_low = ifelse(first_last_25_TIDE.result$gene.expr > 0, "high", "low"), .after = 1)
  write.table(first_last_25_TIDE.result, paste0("./output_TIDE/", gene, "/Modified.TIDE.result.first_last_25.with.",target.gene,".expr.txt" ), sep = ",", row.names = FALSE )

  first_last_30_TIDE.result = first_last_30_TIDE.result %>%
    add_column(gene.expr = t(order_data_TUMOR.normalized[target.gene,rownames(first_last_30_TIDE.result)]), .after = 1)
  first_last_30_TIDE.result = first_last_30_TIDE.result %>%
    add_column(expr_high_or_low = ifelse(first_last_30_TIDE.result$gene.expr > 0, "high", "low"), .after = 1)
  write.table(first_last_30_TIDE.result, paste0("./output_TIDE/", gene, "/Modified.TIDE.result.first_last_30.with.",target.gene,".expr.txt" ), sep = ",", row.names = FALSE )

}

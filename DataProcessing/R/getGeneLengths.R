if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("EDASeq")

# EDAseq ------------------------------------------------------------------
library (EDASeq)
ensembl_list <- c("ENSG00000000003","ENSG00000000419","ENSG00000000457","ENSG00000000460")
getGeneLengthAndGCContent(ensembl_list, "hsa")


ensembl_list <- read.csv(file = 'Chiara_2023_genes.csv')$genes
isid <- grepl("^[Ee][Nn][Ss][A-Za-z]{0,3}[Gg][0-9]+", ensembl_list)
table(isid)
# res <- getGeneLengthAndGCContent(ensembl_list[1:20], "hsa")
res <- getGeneLengthAndGCContent(ensembl_list, "hsa")
tail(res)

write.csv(res, "Chiara_2023_geneLengths_170423.csv", row.names=TRUE)

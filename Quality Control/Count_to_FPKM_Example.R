setwd("~/Desktop/")
library(devtools)
library(biomaRt) #BiocManager::install("biomaRt")
require(gridExtra)
library(grid)
library(Biobase)
library(sva)
library(ggpubr)
#install.packages(c("devtools"))
#source("http://www.bioconductor.org/biocLite.R")
#BiocManager::install(c("Biobase","sva","bladderbatch","snpStats"))
#library(EMA)
library(tidyverse)
library(tidyr)
require(dplyr)
library(rrcov)
#install_github("lazappi/RNAtools")
library(data.table)
library(vroom)
library(tidyr)
library(tibble)
library(grid)
library(vroom)
library(data.table)


counts_to_fpkm <- function(counts, lengths) {
  exp(log(counts) + log(1e9) - log(lengths) - log(sum(counts)) )
}

#FPKM <- apply(count_df, 2, function(x) counts_to_fpkm(x, lengths ))

calculate_FPKM<- function(count_df) { #  change in the code 
  # ensembl_gene_id: based on the input gene format ensembl_gene_id/ hgnc_symbol
  # Name: colname of the gene ids
  df <- count_df
  geneSymbols <-  df$Name

  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl") # hgnc_symbol
  gene_coords=getBM(attributes=c("ensembl_gene_id", "start_position","end_position"), 
                    filters="ensembl_gene_id", values=geneSymbols, mart=human)
  gene_coords$size=gene_coords$end_position - gene_coords$start_position
  # Taking the maximum length for dublicated lengths
  gene_coords %>% group_by(ensembl_gene_id) %>% summarise(size= max(size)) -> gene_coords
  gene_coords <- gene_coords[gene_coords$ensembl_gene_id %in% df$Name,]
  # Taking the maximim of count for dublicated genes
  df %>% group_by(Name) %>% summarize(across(everything(), list(max)) )-> df
  
  df_knw_len <- df[df$Name %in% gene_coords$ensembl_gene_id,]
  df_len <- left_join(gene_coords,df_knw_len,by=c('ensembl_gene_id'="Name"))
  df_len <- df_len[!is.na(df_len$size),]

  df_fpkm <- apply(df_len[,3:ncol(df_len)], 2, function(x) counts_to_fpkm(x, df_len$size))
  df_fpkm <- data.frame(df_fpkm)
  rownames(df_fpkm) <- df_len$ensembl_gene_id
  df_fpkm
}

df_counts <- read_csv('./data_counts.csv')
# Make sure that the first column name in the count ddataframe is called "Name" 
## and has the gene names ("ensembl_gene_id" for enemble gene IDs and "hgnc_symbol" for official gene symbols)
# If the input is official gene symbols (BRCA, etc) replace in the previous function "ensembl_gene_id" with "hgnc_symbol"

colnames(df_counts) <- "Name"
df_fpkm <- calculate_FPKM(df_fpkm)
df_fpkm$ENSG <- row.names(df_fpkm)

write_csv(df_fpkm,'./data_FPKM.csv')

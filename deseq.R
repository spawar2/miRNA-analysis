if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
setwd("/Users/yalegenomecenter/Desktop")
data <- read.csv("FInal_Analysis.csv")
metadata <- read.csv("metadata.csv")
mycounts <- as.data.frame(data)
metadata <- as.data.frame(metadata)
head(mycounts)

#Comparisons
mycounts <- as.data.frame(data[,c(1,8,9,10,11,12,13,14,16,17,18,19,20,21)])
metadata <- as.data.frame(metadata[c(6,7,8,9,10,11,12,13,14,15,16,17,18),])
names(mycounts)[-1]
metadata$id
names(mycounts)[-1]==metadata$id
all(names(mycounts)[-1]==metadata$id)
dds <- DESeqDataSetFromMatrix(countData=mycounts, 
                              colData=metadata, 
                              design=~treatment, 
                              tidy=TRUE)                         
dds <- DESeq(dds)
sizeFactors(dds)
dispersions(dds)
results(dds)

res <- results(dds, tidy=TRUE)
res <- tbl_df(res)
res %>% 
  filter(padj<0.05) %>% 
  write_csv("SAMP8_Ct_Cx8.csv_vs_SAMP8_UNC0642.csv")




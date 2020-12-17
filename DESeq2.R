library(DESeq2)
library(ggplot2)
library(ggrepel)
Its <- read.csv("WT_dAGR_total_peaks.csv")
colData <- read_excel("condition.xlsx")
temp_1 <- Its$peak
temp_2 <- colData$X__1
Its <- Its[,-1]
colData <- colData[,-1]
rownames(Its) <- temp_1
rownames(colData) <- temp_2
Its[1:6] <- lapply(Its[1:6], as.integer)
dds <- DESeqDataSetFromMatrix(countData = Its,colData = colData,design = ~ condition)
dds$condition <- relevel(dds$condition, ref = 'wild-type')
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file="WT_dAGR_DE.csv")



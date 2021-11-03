### DE_analysis

### DESeq2.R

library(DESeq2)

colData <- data.frame(row.names = colnames(expr),group_list = group_list)
dds <- DESeqDataSetFromMatrix(countData = expr,colData = colData,design = ~group_list)
dds2 <- DESeq(dds)

contr <- c("Treat","Control")

tmp <- results(dds2, contrast = c("group_list",contr))
DEG_DESeq2 <- as.data.frame(tmp[order(tmp$padj),])
head(DEG_DESeq2)

DEG_DESeq2 <- na.omit(DEG_DESeq2)

fc_cutoff <- 2
fdr <- 0.05

DEG_DESeq2$regulated <- "normal"

loc_up <- intersect(which(DEG_DESeq2$log2FoldChange>log2(fc_cutoff)),which(DEG_DESeq2$padj<fdr))
loc_down <- intersect(which(DEG_DESeq2$log2FoldChange<(-log2(fc_cutoff))),which(DEG_DESeq2$padj<fdr))

DEG_DESeq2$regulated[loc_up] <- "up"
DEG_DESeq2$regulated[loc_down] <- "down"

table(DEG_DESeq2$regulated)

dir.create("./DESeq2_result")
library(openxlsx)
write.xlsx(DEG_DESeq2,file = "./DESeq2_result/DEG_DESeq2_all.xlsx")
save(DEG_DESeq2,file = "./DESeq2_result/DEG_DEGSeq2_all.rdata")

# DEG_DESeq2 <- DEG_DESeq2[,c(1,2,4,7)]
# save(DEG_DESeq2,file = "./DESeq2_result/DEG_DEGSeq2_nrDEG.rdata")

head(DEG_DESeq2)

suppressPackageStartupMessages(library(SummarizedExperiment))
library(fishpond)
library(DESeq2)
types <- c("txp","oracle","gene","tss")
for (t in types) {
  load(file=paste0("data/se_",t,".rda"))
  y <- se; rm(se)
  assays(y) <- assays(y)[1:3] # counts, TPM, length
  y <- labelKeep(y)
  y <- y[mcols(y)$keep,]
  y$sample <- factor(y$sample)
  dds <- DESeqDataSet(y, ~sample + allele)
  dds <- DESeq(dds, test="LRT", reduced=~sample, fitType="mean", quiet=TRUE)
  res <- results(dds, independentFiltering=FALSE)
  df <- data.frame(qvalue=res$padj, row.names=rownames(res))
  write.table(df, file=paste0("res/deseq2_",t,".tsv"), sep="\t",quote=FALSE)
}

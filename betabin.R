suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
library(fishpond)

types <- c("txp","oracle","gene","tss")
for (t in types) {
  load(file=paste0("data/se_",t,".rda"))
  y <- wide
  y <- labelKeep(y)
  y <- y[mcols(y)$keep,] 
  set.seed(1)
  y <- swish(y, x="allele", pair="sample")
  mcols(y)$keep <- NULL
  write.table(mcols(y), file=paste0("res/bb_",t,".tsv"), sep="\t",quote=FALSE)
}

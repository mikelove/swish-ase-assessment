suppressPackageStartupMessages(library(SummarizedExperiment))
library(fishpond)
library(tximeta)

# import tx2gene-------------------------------------------------------------
dir <- "quants"
samps <- list.files(dir)
files <- file.path(dir, samps, "quant.sf")

t2g <- read.table("tx2gene.oracle.tsv", header=TRUE)
table(table(sub("_[MP]","",t2g[,1]))) # check that we have 2x txps
se <- tximeta(files, skipMeta=TRUE, txOut=FALSE, tx2gene=t2g)

assayNames(se)

# wide format---------------------------------------------------------------------
ntxp <- nrow(se)/2
n <- ncol(se)
cts <- assay(se, "counts")

subset <- seq(1, length(rownames(se)), 2)
txp_nms_m <- rownames(se)[subset] 
txp_nms_p <- sub("_M","_P",txp_nms_m)
txp_nms <- sub("_M","",txp_nms_m)
cts_wide <- cbind(cts[txp_nms_m,], cts[txp_nms_p,])
rownames(cts_wide) <- txp_nms
coldata_wide <- data.frame(sample=factor(c(1:n, 1:n)),
                           allele=factor(rep(c("a1","a2"), each=n)))
rownames(coldata_wide) <- with(coldata_wide, paste0(sample, "-", allele))
wide <- SummarizedExperiment(assays=list(counts=cts_wide),
                             colData=coldata_wide)

for (j in 1:20) {
  idx <- paste0("infRep",j)
  ir_cts <- assay(se, idx)
  ir_cts_wide <- cbind(ir_cts[txp_nms_m,], ir_cts[txp_nms_p,])
  rownames(ir_cts_wide) <- txp_nms
  colnames(ir_cts_wide) <- rownames(coldata_wide)
  assay(wide, idx) <- ir_cts_wide
}

assayNames(wide)

save(wide, file="se.oracle.rda")

y <- wide
y <- labelKeep(y)
y <- y[mcols(y)$keep,] # keep the label doesn't make difference in the iCOBRA
set.seed(1)
y <- swish(y, x="allele", pair="sample", nperms=64)
mcols(y)$log10mean <- log10(rowMeans(assay(y)))
table(mcols(y)$qvalue < .05)
hist(mcols(y)$pvalue, col="grey")

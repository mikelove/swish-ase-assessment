suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
library(devtools)
load_all("../fishpond/fishpond")

# import tx2gene-------------------------------------------------------------
dir <- "boot_quants"
samps <- list.files(dir)
files <- file.path(dir, samps, "quant.sf")

t2g <- read.table("tx2gene.oracle.tsv", header=TRUE)
table(table(sub("_[MP]","",t2g[,1]))) # check that we have 2x txps
t2g$txp <- sub("_.*","",t2g$txp)
t2g$gene <- sub("_.*","",t2g$gene)
t2g <- t2g[!duplicated(t2g$txp),]

s <- paste0("samp",seq_along(files))
coldata <- data.frame(files=files, names=s, sample=s)
wide <- importAllelicCounts(coldata, a1="P", a2="M",
                            format="wide",
                            tx2gene=t2g,
                            ignoreAfterBar=TRUE)

assayNames(wide)

save(wide, file="se.oracle.rda")

load("/proj/milovelab/wu/bulk-ase/ase-analysis/Bootstrap_Analysis/wide.oracle.b+.rda")

y <- wide
y <- labelKeep(y)
y <- y[mcols(y)$keep,] 
set.seed(1)
y <- swish(y, x="allele", pair="sample")
table(mcols(y)$qvalue < .01)
hist(mcols(y)$pvalue, col="grey")

sig_txps <- rownames(y)[mcols(y)$qvalue < 0.01]
table(grepl("-2$", sig_txps))

grep("-2$", sig_txps)
sig_txps[grep("-2$", sig_txps)]
table(grepl("-2$", sig_txps))

hist(mcols(y)[sig_txps[grep("-2$", sig_txps)],]$log2FC, breaks=20, border="white", col="grey")

zzz <- mcols(y)[sig_txps[grep("-2$", sig_txps)],]
zzz <- zzz[order(zzz$pvalue),]
head(zzz,10)

txps[txps$gene_id == "FBgn0003372"]

y$s <- c(1:10,1:10)
y <- computeInfRV(y)
plotInfReps(y, idx="FBgn0001234-2", x="s", cov="allele", legend=TRUE, legendPos="bottom", xlab="sample")
plotInfReps(y, idx="FBgn0003372-2", x="s", cov="allele", legend=TRUE, legendPos="bottom", xlab="sample")
plotInfReps(y, idx="FBgn0003373-2", x="s", cov="allele", legend=TRUE, legendPos="bottom", xlab="sample")
plotInfReps(y, idx="FBgn0003943-2", x="s", cov="allele", legend=TRUE, legendPos="bottom", xlab="sample")

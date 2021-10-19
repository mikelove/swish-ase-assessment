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

y <- wide
y <- labelKeep(y)
y <- y[mcols(y)$keep,] 
set.seed(1)
y <- swish(y, x="allele", pair="sample")
table(mcols(y)$qvalue < .1)
hist(mcols(y)$pvalue, col="grey")

head(mcols(y)[mcols(y)$pvalue < 0.2,])

y$s <- factor(c(1:10,1:10))
y <- computeInfRV(y)
plotInfReps(y, idx="FBgn0000008-1.67", x="allele", cov="s", legend=TRUE)

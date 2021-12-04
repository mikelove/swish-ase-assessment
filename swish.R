suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
library(devtools)
load_all("../../fishpond/fishpond")

# import tx2gene
dir <- "../ase-sim/quants"
samps <- list.files(dir)
files <- file.path(dir, samps, "quant.sf")

if (FALSE) {
  load("../ase-sim/granges.rda")
  mcols(txps)$txp_groups <- paste0(txps$gene_id, "-", round(txps$abundance,2))
  t2g <- data.frame(txp=c(names(txps)), gene=c(txps$txp_groups))
  write.table(t2g, file="t2g.oracle.tsv", quote=FALSE, row.names=FALSE)
}

t2g <- read.table("t2g.oracle.tsv", header=TRUE)

s <- paste0("samp",seq_along(files))
coldata <- data.frame(files=files, names=s, sample=s)
wide <- importAllelicCounts(coldata, a1="P", a2="M",
                            format="wide",
                            tx2gene=t2g,
                            ignoreAfterBar=TRUE)

# no summarization
wide <- importAllelicCounts(coldata, a1="P", a2="M",
                            format="wide",
                            ignoreAfterBar=TRUE)

assayNames(wide)

#save(wide, file="se.oracle.rda")
save(wide, file="se.txp.rda")
# load object
#load("/proj/milovelab/wu/bulk-ase/ase-analysis/Bootstrap_Analysis/wide.txp.b+.rda")
#load("/proj/milovelab/wu/bulk-ase/ase-analysis/Bootstrap_Analysis/wide.oracle.b+.rda")
#load("/proj/milovelab/wu/bulk-ase/ase-analysis/Bootstrap_Analysis/wide.oracle.b.rda")

y <- wide
y <- labelKeep(y)
y <- y[mcols(y)$keep,] 
set.seed(1)
y <- swish(y, x="allele", pair="sample")

table(mcols(y)$qvalue < .01)
hist(mcols(y)$pvalue, col="grey")

# oracle
sig_txps <- rownames(y)[mcols(y)$qvalue < 0.01]
table(grepl("-2$", sig_txps))
table(ai=!grepl("-2$", rownames(y)), sig=mcols(y)$qvalue < 0.01)
prop.table(table(ai=!grepl("-2$", rownames(y)), sig=mcols(y)$qvalue < 0.01), 1)

# txp-level
mcols(y)$abund <- mcols(txps)[rownames(y),"abundance"]

sig <- mcols(y)$qvalue < 0.1
ai <- mcols(y)$abund != 2
table(ai, sig)


###

# check structure in plot

load("se.txp.rda")
y <- wide
y <- labelKeep(y)
y <- y[mcols(y)$keep,]
y <- computeInfRV(y)
#y2 <- y[1:2000,]
out <- swish(y, x="allele", pair="sample", returnNull=TRUE)
sigma <- sqrt(rowVars(out$nulls))
mcols(y)$mean <- rowMeans(assay(y))
mcols(y)$abund <- mcols(txps)[rownames(y),"abundance"]
plot(log10(mcols(y)$meanInfRV+.1), log10(sigma), cex=.2, col=ifelse(mcols(y2)$abund == 2,1,2))
plot(log10(mcols(y)$meanInfRV+.1), log10(sigma), cex=.2)
idx <- rowSums(assay(y) == 0) == 1
plot(log10(mcols(y)$meanInfRV+.1), log10(sigma), cex=1, col=ifelse(idx, 2, 1), pch=20)
plot(log10(mcols(y)$meanInfRV+.1), log10(mcols(y)$mean), cex=1, col=ifelse(idx, 2, 1), pch=20)
idx2 <- rowMedians(abs(assay(y)[,1:10] - assay(y)[,11:20])) < 1
plot(log10(mcols(y)$meanInfRV+.1), log10(mcols(y)$mean), cex=1, col=ifelse(idx2, 2, 1), pch=20)
plot(log10(mcols(y)$meanInfRV+.1), log10(sigma), cex=1, col=ifelse(idx2, 2, 1), pch=20)

idx3 <- rownames(y)[idx2]

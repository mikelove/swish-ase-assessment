suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
library(fishpond)

# import tx2gene
dir <- "../ase-sim/quants"
samps <- list.files(dir)
files <- file.path(dir, samps, "quant.sf")

if (FALSE) {
  load("../ase-sim/granges.rda")
  mcols(txps)$txp_groups <- paste0(txps$gene_id, "-", round(txps$abundance,2))
  t2g <- data.frame(txp=names(txps), gene=txps$txp_groups)
  write.table(t2g, file="t2g.oracle.tsv", quote=FALSE, row.names=FALSE)
  mcols(txps)$tss_groups <- paste0(txps$gene_id, "-", round(txps$tss,2))
  t2g <- data.frame(txp=names(txps), gene=txps$tss_groups)
  write.table(t2g, file="t2g.tss.tsv", quote=FALSE, row.names=FALSE)
  t2g <- data.frame(txp=names(txps), gene=txps$gene_id)
  write.table(t2g, file="t2g.gene.tsv", quote=FALSE, row.names=FALSE)
}

t2g <- read.table("t2g.oracle.tsv", header=TRUE)
t2g <- read.table("t2g.tss.tsv", header=TRUE)
t2g <- read.table("t2g.gene.tsv", header=TRUE)

s <- paste0("samp",seq_along(files))
coldata <- data.frame(files=files, names=s, sample=s)
wide <- importAllelicCounts(coldata, a1="P", a2="M",
                            format="wide",
                            tx2gene=t2g,
                            ignoreAfterBar=TRUE)

# no summarization = txp
wide <- importAllelicCounts(coldata, a1="P", a2="M",
                            format="wide",
                            ignoreAfterBar=TRUE)

#suffix <- "oracle"
#suffix <- "txp"
#suffix <- "tss"
suffix <- "gene"

save(wide, file=paste0("data/se_",suffix,".rda"))

load(file=paste0("data/se_",suffix,".rda"))

y <- wide
y <- labelKeep(y)
y <- y[mcols(y)$keep,] 
set.seed(1)
y <- swish(y, x="allele", pair="sample")

mcols(y)$keep <- NULL
write.table(mcols(y), file=paste0("res/",suffix,".tsv"), sep="\t",quote=FALSE)

###

truth <- mcols(txps)[,c("gene_id","tss","abundance","txp_groups","tss_groups")]
#write.table(truth, file="truth.tsv", sep="\t",quote=FALSE)

###


suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
library(devtools)
load_all("../fishpond/fishpond")

# import tx2gene
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

# load object
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
fp <- sub("-2","",sig_txps[grep("-2$", sig_txps)])
table(grepl("-2$", sig_txps))

# check repeat status
library(rtracklayer)
rep <- import("simpleRepeat_dm6_Aug2014.bed.gz")
seqlevelsStyle(rep) <- "NCBI"
rep <- rep[seqnames(rep) %in% c("2L","2R","3L","3R","4","X")]
rep <- keepStandardChromosomes(rep, pruning.mode="coarse")
library(AnnotationHub)
ah <- AnnotationHub()
release <- 100
q <- query(ah, c("EnsDb",release,"Drosophila"))
library(ensembldb)
edb <- q[[1]]
g <- exonsBy(edb, by="gene")
g <- keepStandardChromosomes(g, pruning.mode="coarse")
g_fp <- g[fp]
rep_size <- 50
big_rep <- rep[width(rep) >= rep_size]
table(overlapsAny(g_fp, big_rep, minoverlap=rep_size))
table(overlapsAny(g, big_rep, minoverlap=rep_size))
g_fp_extra <- g_fp[!overlapsAny(g_fp, big_rep, minoverlap=rep_size)]
names(g_fp_extra)

g0 <- genes(edb)
dat <- cbind(as.data.frame(g0[names(g_fp_extra)]),
             as.data.frame(mcols(y)[paste0(names(g_fp_extra),"-2"),]))
dat[dat$pvalue < .0001,]
table(dat[dat$pvalue < .0001,]$gene_biotype == "pseudogene")

###

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

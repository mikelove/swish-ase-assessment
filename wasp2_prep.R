library(readr)
library(dplyr)
dir <- "../ase-sim/wasp2_results"
samps <- list.files(dir)
files <- list.files(file.path(dir, samps[1]))
x <- read_tsv(file.path(dir, samps[1], files[1]), show_col_types=FALSE)
genes <- x$genes
for (s in samps) {
  x <- read_tsv(file.path(dir, s, files[1]), show_col_types=FALSE)
  genes <- union(genes, setdiff(x$genes, genes))
}
pvals <- data.frame(row.names=genes)
for (f in files) {
  for (s in samps) {
    x <- read_tsv(file.path(dir, s, f), show_col_types=FALSE)
    stopifnot(all(x$genes %in% genes))
    pvals[s] <- 1
    pvals[x$genes,s] <- x$pval
  }
  filename <- file.path("wasp2", sub("as_results_gene_", "", f))
  write.table(pvals, filename, sep="\t", quote=FALSE)
}

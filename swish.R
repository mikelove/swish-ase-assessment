# Dec 20 2021

suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
library(fishpond)

############################################
## read in the sample names from 'quants' ##
############################################

dir <- "../ase-sim/quants"
samps <- list.files(dir)
files <- file.path(dir, samps, "quant.sf")

###############################################################
## define different levels of aggregation, using granges.rda ##
###############################################################

if (FALSE) {
  load("../ase-sim/granges.rda")
  mcols(txps)$oracle_groups <- paste0(txps$gene_id, "-", round(txps$abundance,2))
  mcols(txps)$tss_groups <- paste0(txps$gene_id, "-", round(txps$tss,2))
  mcols(txps)$gene_groups <- txps$gene_id
  for (x in c("oracle","tss","gene")) {
    t2g <- data.frame(txp=names(txps),
                      gene=mcols(txps)[paste0(x, "_groups")])
    # these are saved as, e.g. t2g_oracle.tsv, etc.
    write.table(t2g, file=paste0("t2g_",x,".tsv"), quote=FALSE, row.names=FALSE)
  }
}

###############################################
## import counts at different level, save SE ##
###############################################

types <- c("txp","oracle","gene","tss")
for (t in types) {
  s <- paste0("samp",seq_along(files))
  coldata <- data.frame(files=files, names=s, sample=s)
  if (t == "txp") {
    # no summarization = txp
    wide <- importAllelicCounts(coldata, a1="P", a2="M",
                                format="wide",
                                ignoreAfterBar=TRUE)
  } else {
    t2g <- read.table(paste0("t2g_",t,".tsv"), header=TRUE)
    wide <- importAllelicCounts(coldata, a1="P", a2="M",
                                format="wide",
                                tx2gene=t2g,
                                ignoreAfterBar=TRUE)
  }
  # save the SE to e.g. 'data/se_oracle.rda'
  save(wide, file=paste0("data/se_",t,".rda"))
}

#######################################
## run swish and write results table ##
#######################################

types <- c("txp","oracle","gene","tss")
for (t in types) {
  load(file=paste0("data/se_",t,".rda"))
  y <- wide
  y <- labelKeep(y)
  y <- y[mcols(y)$keep,] 
  set.seed(1)
  y <- swish(y, x="allele", pair="sample")
  mcols(y)$keep <- NULL
  # write results to e.g. 'res/oracle.tsv'
  write.table(mcols(y), file=paste0("res/",t,".tsv"), sep="\t",quote=FALSE)
}

###############################################################
## write a table that can be used to define truth for iCOBRA ##
###############################################################

truth <- mcols(txps)[,c("gene_id","tss","abundance","oracle_groups","tss_groups","isoAI","geneAI")]
write.table(truth, file="truth.tsv", sep="\t",quote=FALSE)

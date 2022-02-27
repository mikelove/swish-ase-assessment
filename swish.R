# Feb 27 2022

suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
devtools::load_all("../../fishpond/fishpond")

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
  library(plyranges)
  load("../ase-sim/granges.rda")
  txps <- txps %>%
    mutate(oracle_groups = paste0(gene_id, "-", round(abundance,2)),
           tss_groups = paste0(gene_id, "-", round(tss,2)),
           gene_groups = gene_id)
  save(txps, file="granges_with_groups.rda")
  ## for (x in c("oracle","tss","gene")) {
  ##   t2g <- data.frame(txp=names(txps),
  ##                     gene=mcols(txps)[paste0(x, "_groups")])
  ##   # these are saved as, e.g. t2g_oracle.tsv, etc.
  ##   write.table(t2g, file=paste0("t2g_",x,".tsv"), quote=FALSE, row.names=FALSE)
  ## }
}

###############################################
## import counts at different level, save SE ##
###############################################

load("granges_with_groups.rda")

types <- c("txp","oracle","gene","tss")
for (t in types) {
  s <- paste0("samp",seq_along(files))
  coldata <- data.frame(files=files, names=s, sample=s)
  if (t == "txp") {
    # no summarization = txp
    se <- importAllelicCounts(coldata, a1="P", a2="M",
                              format="wide",

                              ignoreAfterBar=TRUE)
    rowRanges(se) <- txps[rownames(se)]
  } else {
    #t2g <- read.table(paste0("t2g_",t,".tsv"), header=TRUE)
    mcols(txps)$group_id <- mcols(txps)[[paste0(t,"_groups")]]
    se <- importAllelicCounts(coldata, a1="P", a2="M",
                              format="wide",
                              tx2gene=txps,
                              ignoreAfterBar=TRUE)
  }
  # save the SE to e.g. 'data/se_oracle.rda'
  save(se, file=paste0("data/se_",t,".rda"))
}

#######################################
## run swish and write results table ##
#######################################

types <- c("txp","oracle","gene","tss")
for (t in types) {
  load(file=paste0("data/se_",t,".rda"))
  y <- se
  y <- labelKeep(y)
  y <- y[mcols(y)$keep,] 
  set.seed(1)
  y <- swish(y, x="allele", pair="sample")
  mcols(y)$keep <- NULL
  # write results to e.g. 'res/oracle.tsv'
  write.table(mcols(y), file=paste0("res/",t,".tsv"), sep="\t",quote=FALSE)
  save(y, file=paste0("res/",t,"_se.rda"))
}

###############################################################
## write a table that can be used to define truth for iCOBRA ##
###############################################################

columns <- c("gene_id","tss","abundance","oracle_groups",
             "tss_groups","isoAI","geneAI")
truth <- mcols(txps)[,columns]
write.table(truth, file="truth.tsv", sep="\t",quote=FALSE)

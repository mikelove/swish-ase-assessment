# Jan 4 2022

library(iCOBRA)
library(dplyr)
library(tibble)

cols <- palette.colors()[1:6]
types <- c("txp","gene","tss","oracle","mmdiff")
names(cols)[1:5] <- types
names(cols)[6] <- "truth"

# Motivation for this script:

# For iCOBRA need to make a 'truth' table and a 'padj' table.
# This is complicated by the fact that we have different levels
# of resolution of analysis. A simpliciation is to just
# evaluate everything at transcript level. We can propagate
# the q-values from inner nodes to the leaves and evaluate
# transcript level sensitivity and FDR. We do this with
# inner_join()

# this TSV is made in swish.R
truth <- read.delim("truth.tsv")
truth <- truth %>%
  mutate(gene_groups = gene_id,
         status = ifelse(abundance == 2, 0, 1),
         group = case_when(isoAI ~ "iso", geneAI ~ "gene", TRUE ~ "null"))
# iCOBRA wants a data.frame with rownames
# but we need rownames as a column for dplyr's inner_join()
truth_tb <- truth %>% rownames_to_column("txp") %>% tibble()

# make an empty table which we will add columns to
padj <- data.frame(row.names=rownames(truth))

# loop of diferrent levels of analysis:
for (t in setdiff(types, "mmdiff")) {
  res <- read.delim(paste0("res/",t,".tsv")) # read swish results
  if (t == "txp") {
    # directly add the qvalues from txp-level
    padj$txp <- 1
    padj[rownames(res),"txp"] <- res$qvalue
  } else {
    # first join inner-node anlaysis to the txp-level table
    res_tb <- res %>% select(qvalue) %>%
      rownames_to_column(paste0(t,"_groups")) %>% # this column will allow joining
      tibble() %>%
      inner_join(truth_tb) # join with the truth table which has txp ID
    padj[[t]] <- 1 # pre-populate the column with 1's
    padj[res_tb$txp, t] <- res_tb$qvalue # add the q-values for matching rows
  }
}

# add mmdiff results
mmdiff <- read.table("../ase-sim/mmseq/mmdiff_results.txt", header=TRUE)
padj$mmdiff <- 1
padj[mmdiff$feature_id,"mmdiff"] <- 1 - mmdiff$posterior_probability

# ok now 'padj' is done and we can just run iCOBRA directly:

cd <- COBRAData(padj=padj, truth=truth)
cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                     "fdrtprcurve","overlap"),
                            splv="group", # breaks across isoAI / geneAI / null
                            thrs=c(.01,.05,.1),
                            thr_venn=.05)
cplot <- prepare_data_for_plot(cp, colorscheme=cols)
xrng <- c(0,.15)
yrng <- c(0,1)
plot_fdrtprcurve(cplot,
                 plottype="points",
                 xaxisrange=xrng,
                 yaxisrange=yrng,
                 title="txp-level AI testing")

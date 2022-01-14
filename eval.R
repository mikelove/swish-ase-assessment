# Jan 4 2022

library(iCOBRA)
library(dplyr)
library(tibble)

cols <- palette.colors(9)[c(1:4,6:8)]
types <- c("txp","gene","tss","oracle","mmdiff","mmdiff_gene")
names(cols)[1:(length(cols)-1)] <- types
names(cols)[length(cols)] <- "truth"

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
         mmdiff_gene_groups = gene_id,
         status = ifelse(abundance == 2, 0, 1),
         AI = case_when(isoAI ~ "iso", geneAI ~ "gene", TRUE ~ "null"))
# iCOBRA wants a data.frame with rownames
# but we need rownames as a column for dplyr's inner_join()
truth_tb <- truth %>% rownames_to_column("txp") %>% tibble()

# make an empty table which we will add columns to
padj <- data.frame(row.names=rownames(truth))

# loop of diferrent levels of analysis:
for (t in types) {
  if (!grepl("mmdiff", t)) {
    res <- read.delim(paste0("res/",t,".tsv")) # read swish results
  } else {
    res <- read.table(paste0("../ase-sim/mmseq/",t,"_results.txt"), header=TRUE)
    rownames(res) <- res$feature_id
    res$qvalue <- 1 - res$posterior_probability
  }
  if (t %in% c("txp","mmdiff")) {
    # directly add the qvalues from txp-level
    padj[[t]] <- 1
    padj[rownames(res),t] <- res$qvalue
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

# ok now 'padj' is done and we can just run iCOBRA directly:

cd <- COBRAData(padj=padj, truth=truth)
cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                     "fdrtprcurve","overlap"),
                            splv="AI", # breaks across isoAI / geneAI / null
                            thrs=c(.01,.05,.1),
                            thr_venn=.05)
cplot <- prepare_data_for_plot(cp, colorscheme=cols)
plot_tpr(cplot, pointsize=2)

cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                     "fdrtprcurve","overlap"),
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

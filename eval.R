library(iCOBRA)

cols <- palette.colors()[1:5]
types <- c("txp","gene","tss","oracle")
names(cols)[1:4] <- types
names(cols)[5] <- "truth"

library(dplyr)
library(tibble)

# For iCOBRA need to make a truth table and a padj table
# here complicated by the fact that we have different levels
# of resolution of analysis. A simpliciation is to just
# evaluate everything at transcript level. We can propagate
# the q-values from inner nodes to the leaves and evaluate
# transcript level sensitivity and FDR
truth <- read.delim("truth.tsv")
truth <- truth %>%
  mutate(gene_groups = gene_id,
         status = ifelse(abundance == 2, 0, 1),
         group = case_when(isoAI ~ "iso", geneAI ~ "gene", TRUE ~ "null"))
truth_tb <- truth %>% rownames_to_column("txp") %>% tibble()

padj <- data.frame(row.names=rownames(truth))

for (t in types) {
  res <- read.delim(paste0("res/",t,".tsv"))
  #res <- read.delim(paste0("res/bb_",t,".tsv"))
  #res <- read.delim(paste0("res/deseq2_",t,".tsv"))
  if (t == "txp") {
    # just directly add the qvalues from txp-level
    padj$txp <- 1
    padj[rownames(res),"txp"] <- res$qvalue
  } else {
    # first join inner-node anlaysis to the txp-level table
    res_tb <- res %>% select(qvalue) %>%
      rownames_to_column(paste0(t,"_groups")) %>% # this will allow joining
      tibble() %>%
      inner_join(truth_tb) # join with the truth table which has txp ID
    padj[[t]] <- 1 # pre-populate with 1
    padj[res_tb$txp, t] <- res_tb$qvalue # add the q-values when we have them
  }
}

cd <- COBRAData(padj=padj, truth=truth)
cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr","fdrtprcurve","overlap"),
                            splv="group",
                            thrs=c(.01,.05,.1),
                            thr_venn=.05)

cp <- prepare_data_for_plot(cp, colorscheme=cols)

# simple plot
plot_tpr(cp, pointsize=2.5)

# TPR over FDR
yrng <- c(0,1)
xrng <- c(0,.15)
#xrng <- c(0,.2)
plot_fdrtprcurve(cp,
                 plottype="points",
                 xaxisrange=xrng,
                 yaxisrange=yrng,
                 title="txp-level AI testing")

if (FALSE) {
  # FDR and number
  plot_fdrnbrcurve(cp, xaxisrange=xrng) +
    ggplot2::ylim(0,4500)
  
  # upset plot
  plot_upset(cp, order.by="freq", nintersects=10)
}

library(iCOBRA)

cols <- palette()[2:6]
meths <- c("txp","gene","tss","oracle")
names(cols)[1:4] <- meths
names(cols)[5] <- "truth"

library(dplyr)
library(tibble)

truth <- read.delim("truth.tsv")
truth <- truth %>% mutate(status = ifelse(abundance == 2, 0, 1))
truth_tb <- truth %>% rownames_to_column("txp") %>% tibble()

padj <- data.frame(row.names=rownames(truth))

match_col <- c(gene="gene_id", tss="tss_groups", oracle="txp_groups")

for (m in meths) {
  res <- read.delim(paste0("res/",m,".tsv"))
  if (m == "txp") {
    # just directly add the qvalues
    padj[[m]] <- 1
    padj[rownames(res),m] <- res$qvalue
  } else {
    # first join with the txp-level information
    res_tb <- res %>% select(qvalue) %>%
      rownames_to_column(match_col[m]) %>%
      tibble() %>%
      inner_join(truth_tb)
    padj[[m]] <- 1
    padj[res_tb$txp,m] <- res_tb$qvalue
  }
}

cd <- COBRAData(padj=padj, truth=truth)
cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr","fdrtprcurve","overlap"),
                            thrs=c(.01,.05,.1),
                            thr_venn=.05)

cobraplot <- prepare_data_for_plot(cp,
                                   colorscheme=cols,
                                   facetted=TRUE,
                                   incloverall=TRUE)

# simple plot
plot_tpr(cobraplot, pointsize=2.5)

# TPR over FDR
yrng <- c(0,1)
xrng <- c(0,.5)
plot_fdrtprcurve(cobraplot,
                 plottype="points",
                 pointsize=2.5,
                 xaxisrange=xrng,
                 yaxisrange=yrng,
                 title="txp-level AI testing")

# FDR and number
plot_fdrnbrcurve(prepare_data_for_plot(cp, colorscheme=cols),
                 xaxisrange=xrng) +
  ggplot2::ylim(0,4500)

# upset plot
plot_upset(prepare_data_for_plot(cp, colorscheme=cols),
           order.by="freq", nintersects=10)

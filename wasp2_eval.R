single <- read.delim("wasp2/single.tsv")
linear <- read.delim("wasp2/linear.tsv")

load("../ase-sim/granges.rda")
library(GenomicRanges)
mapdf <- data.frame(symbol=g$symbol, id=g$gene_id)
table(rownames(single) %in% mapdf$symbol)
table(rownames(linear) %in% mapdf$symbol)
single <- single[rownames(single) %in% mapdf$symbol,]
linear <- linear[rownames(linear) %in% mapdf$symbol,]

padj$mmdiff <- NULL
padj$mmdiff_gene <- NULL

rownames(single) <- mapdf$id[match(rownames(single), mapdf$symbol)]
rownames(linear) <- mapdf$id[match(rownames(linear), mapdf$symbol)]

fisherP <- function(p) {
  pchisq(-2 * sum(log(p)), 2*length(p), lower.tail=FALSE)
}

single_padj <- tibble(gene_id=rownames(single),
                      wasp2_s=p.adjust(apply(single, 1, fisherP), method="BH"))
linear_padj <- tibble(gene_id=rownames(linear),
                      wasp2_l=p.adjust(apply(linear, 1, fisherP), method="BH"))

dat <- truth_tb %>%
  select(txp, gene_id) %>%
  left_join(single_padj) %>%
  left_join(linear_padj)

table(is.na(dat$wasp2_s))
table(is.na(dat$wasp2_l))

library(tidyr)
dat <- dat %>% replace_na(list(wasp2_s=1, wasp2_l=1))

padj$wasp2_s <- 1
padj$wasp2_l <- 1
padj[dat$txp,"wasp2_s"] <- dat$wasp2_s
padj[dat$txp,"wasp2_l"] <- dat$wasp2_l

names(cols)[names(cols) %in% c("mmdiff","mmdiff_gene")] <- c("wasp2_s","wasp2_l")
cols <- cols[c(1,3,2,4:7)]

cd <- COBRAData(padj=padj, truth=truth)
cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                     "fdrtprcurve","overlap"),
                            splv="AI", # breaks across isoAI / geneAI / null
                            thrs=c(.01,.05,.1),
                            thr_venn=.05)
cplot <- prepare_data_for_plot(cp, colorscheme=cols)
cplot <- reorder_levels(cplot, names(cols))
pdf(file="figs/sim_with_wasp2_tpr.pdf")
plot_tpr(cplot, pointsize=2)
dev.off()

cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                     "fdrtprcurve","overlap"),
                            thrs=c(.01,.05,.1),
                            thr_venn=.05)
cplot <- prepare_data_for_plot(cp, colorscheme=cols)
cplot <- reorder_levels(cplot, names(cols))
xrng <- c(0,.3)
yrng <- c(0,1)
pdf(file="figs/sim_with_wasp2_fdrtpr.pdf")
plot_fdrtprcurve(cplot,
                 plottype="points",
                 xaxisrange=xrng,
                 yaxisrange=yrng)
dev.off()

pdf(file="figs/sim_with_wasp2_upset.pdf")
plot_upset(cplot, order.by="freq", decreasing=TRUE, nintersects=12)
dev.off()

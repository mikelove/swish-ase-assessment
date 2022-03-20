library(tibble)
library(dplyr)
library(ggplot2)

wasp <- read.table("../ase-sim/wasp_cht/cht_results.txt", header=TRUE)
wasp <- wasp %>% mutate(ratio=ALT.AS.READ.COUNT/REF.AS.READ.COUNT) %>%
  select(chr=TEST.SNP.CHROM, start=REGION.START, end=REGION.END,
         snp=TEST.SNP.POS, ratio,
         total=TOTAL.AS.READ.COUNT, pvalue=P.VALUE) %>%
  tibble()
# remove duplicate rows
wasp <- wasp %>% mutate(id = paste(chr, snp, start, end, sep="-")) %>%
  filter(!duplicated(id))
sum(duplicated(wasp$snp))

# plot
## wasp %>% filter(total > 50 & total < 1e5) %>%
##   ggplot(aes(total, log2(ratio), col=-log10(pvalue+1e-15))) +
##   geom_point() + coord_cartesian(ylim=c(-.5, .5))

load("../ase-sim/data/drosophila_test_target.rda")

test_target <- wasp_test_target %>% mutate(id = paste(chr, pos, tstart, tend, sep="-")) %>%
  select(id, gene_id) %>% tibble()
# this adds gene_id
wasp <- wasp %>% left_join(test_target)

# all matches, some duplicate tests
table(is.na(wasp$gene_id))
sum(duplicated(wasp$id))
wasp <- wasp %>% filter(!duplicated(id))
sum(duplicated(wasp$gene_id))

# cutoff at min p
min_p <- min(wasp$pvalue[wasp$pvalue > 0])
wasp <- wasp %>%
  mutate(pvalue=ifelse(pvalue == 0, min_p, pvalue))

load("../ase-sim/granges.rda")

# check on overlaps for FPs
library(plyranges)
gene_status <- g %>% select(gene_id, isoAI, geneAI, .drop_ranges=TRUE) %>%
  as.data.frame() %>% tibble()
wasp <- wasp %>% left_join(gene_status)
fp_gene_ids <- wasp %>% filter(pvalue < .1 & !isoAI & !geneAI) %>% pull(gene_id)
fp_genes <- g[fp_gene_ids] %>% select(isoAI, geneAI)
ai_genes <- g %>% filter(isoAI | geneAI) %>%
  select(isoAI, geneAI)

overlap_fps <- fp_genes %>% join_overlap_inner(ai_genes, maxgap=100) %>% names()
length(overlap_fps)

# re-do locfdr calculation
redo <- TRUE
if (redo) {
  library(locfdr)
  z <- qnorm(wasp$pvalue, lower.tail=FALSE)
  fit <- locfdr(z) # check plot to examine left tail
  fdr <- fit$fdr
  fdr[z < 0] <- 1
  wasp$qvalue <- fdr
  # wasp %>% ggplot(aes(qvalue, fdr)) + geom_point()
} else {
  wasp$qvalue <- p.adjust(wasp$pvalue, method="BH")
}

# propagate q-values to txps
t2g <- txps %>% select(tx_id, gene_id, .drop_ranges=TRUE) %>%
  as.data.frame %>% tibble()
wasp_qval <- wasp %>% select(gene_id, qvalue) %>%
    mutate(qvalue = ifelse(gene_id %in% overlap_fps, 1, qvalue)) %>%
  left_join(t2g) %>%
  select(qvalue, tx_id)

### need to run eval.R script to build 'padj' and 'truth' ###

padj$mmdiff <- NULL
padj$mmdiff_gene <- NULL

padj$wasp <- 1
padj[wasp_qval$tx_id,"wasp"] <- wasp_qval$qvalue

cols <- palette.colors()[c(1,3,6,7,5,8)]
types <- c("txp","gene","tss","oracle","wasp","truth")
names(cols) <- types

library(iCOBRA)

cd <- COBRAData(padj=padj, truth=truth)

cp <- calculate_performance(
  cd, binary_truth="status", aspect="tpr", splv="AI",
  thrs=c(.01,.05,.1), thr_venn=.05)
cplot <- prepare_data_for_plot(cp, colorscheme=cols)
cplot <- reorder_levels(cplot, names(cols))
pdf("figs/sim_with_wasp_tpr.pdf")
plot_tpr(cplot, pointsize=2.5)
dev.off()

cp <- calculate_performance(
  cd, binary_truth="status",
  aspect=c("tpr","fdr","fdrtpr","fdrnbr",
           "fdrtprcurve","overlap"),
  thrs=c(.01,.05,.1), thr_venn=.05)
cplot <- prepare_data_for_plot(cp, colorscheme=cols)
cplot <- reorder_levels(cplot, names(cols))
xrng <- c(0,.15)
yrng <- c(0,1)
pdf("figs/sim_with_wasp_fdrtpr_no_locfdr.pdf")
plot_fdrtprcurve(cplot,
                 plottype="points",
                 xaxisrange=xrng,
                 yaxisrange=yrng)
dev.off()

pdf(file="figs/sim_with_wasp_upset.pdf")
plot_upset(cplot, order.by="freq", decreasing=TRUE, nintersects=8)
dev.off()

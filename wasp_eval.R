cols <- palette.colors()[1:6]
types <- c("txp","gene","tss","oracle")
names(cols)[1:4] <- types
names(cols)[5] <- "truth"
names(cols)[6] <- "wasp"

library(dplyr)
library(ggplot2)
wasp <- read.table("../ase-sim/wasp_cht/cht_results.txt", header=TRUE)
wasp <- wasp %>% mutate(
                   start=as.numeric(sapply(strsplit(REGION.START,";"), head, 1)),
                   end=as.numeric(sapply(strsplit(REGION.END,";"), tail, 1)),
                   ratio=ALT.AS.READ.COUNT/REF.AS.READ.COUNT) %>%
  select(seqnames=TEST.SNP.CHROM, start, end,
         snp=TEST.SNP.POS, ratio,
         total=TOTAL.AS.READ.COUNT, pvalue=P.VALUE) %>%
  tibble()
# remove duplicate rows
wasp <- wasp %>% mutate(id = paste(seqnames, start, end, snp, sep="-")) %>%
  filter(!duplicated(id))
sum(duplicated(wasp$snp))

# plot
## wasp %>% filter(total > 50 & total < 1e5) %>%
##   ggplot(aes(total, log2(ratio), col=-log10(pvalue+1e-15))) +
##   geom_point() + coord_cartesian(ylim=c(-.5, .5))

# load truth to identify gene from WASP SNP-level results
library(plyranges)
load("../ase-sim/granges.rda")
txps <- txps %>%
  select(tx_id, gene_id, snp_loc, isoAI, geneAI) %>%
  mutate(snp = min(snp_loc))
strand(txps) <- "*"
# ~1 minute
system.time({
  genes <- txps %>%
    group_by(gene_id) %>%
    reduce_ranges(snp=min(snp), isoAI=any(isoAI), geneAI=any(geneAI))
})
genes_tb <- genes %>% as.data.frame() %>%
  mutate(id = paste(seqnames, start, end, snp, sep="-")) %>%
  tibble()

# this adds gene_id and other variables
wasp <- wasp %>% left_join(genes_tb)

# all matches, no duplicates
table(is.na(wasp$gene_id))
sum(duplicated(wasp$gene_id))

# create q-value
min_p <- min(wasp$pvalue[wasp$pvalue > 0])
wasp <- wasp %>%
  mutate(pvalue=ifelse(pvalue == 0, min_p, pvalue), 
         qvalue=p.adjust(pvalue, method="BH"))

# create locfdr
z <- qnorm(wasp$pvalue, lower.tail=FALSE)
library(locfdr)
fit <- locfdr(z)
fdr <- fit$fdr
fdr[z < 0] <- 1
wasp$locfdr <- fdr
wasp %>% ggplot(aes(qvalue, fdr)) + geom_point()

# check on overlaps for FPs
fp_gene_ids <- wasp %>% filter(pvalue < .1 & !isoAI & !geneAI) %>% pull(gene_id)
fp_genes <- g[fp_gene_ids] %>% select(isoAI, geneAI)
ai_genes <- genes %>% filter(isoAI | geneAI) %>% select(isoAI, geneAI)

overlap_fps <- fp_genes %>% join_overlap_inner(ai_genes, maxgap=100) %>% names()
length(overlap_fps)

wasp %>%
  mutate(pvalue = ifelse(gene_id %in% overlap_fps, NA, pvalue)) %>%
  mutate(AI=factor(isoAI | geneAI, levels=c("TRUE","FALSE"))) %>% 
  ggplot(aes(x=pvalue, fill=AI)) + geom_histogram() +
  ggtitle("WASP gene-level p-values")

# propagate locfdr / q-values to txps
t2g <- txps %>% select(tx_id, gene_id, .drop_ranges=TRUE) %>%
  as.data.frame %>% tibble()
wasp_qval <- wasp %>% select(gene_id, qvalue=locfdr) %>%
    mutate(qvalue = ifelse(gene_id %in% overlap_fps, 1, qvalue)) %>%
  left_join(t2g) %>%
  select(qvalue, tx_id)

# need to run eval.R script to build 'padj' and 'truth'

padj$wasp <- 1
padj[wasp_qval$tx_id,"wasp"] <- wasp_qval$qvalue

cd <- COBRAData(padj=padj, truth=truth)
cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                     "fdrtprcurve","overlap"),
                            thrs=c(.01,.05,.1),
                            thr_venn=.05)

cplot <- prepare_data_for_plot(cp, colorscheme=cols)

# simple plot
plot_tpr(cplot, pointsize=2.5)

xrng <- c(0,.15)
yrng <- c(0,1)
plot_fdrtprcurve(cplot,
                 plottype="points",
                 xaxisrange=xrng,
                 yaxisrange=yrng,
                 title="txp-level AI testing")

# FDR and number
plot_fdrnbrcurve(cplot, xaxisrange=xrng) +
  ggplot2::ylim(0,4500)

# upset plot
plot_upset(cplot, order.by="freq", nintersects=10)

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

wasp %>% as.data.frame() %>%
  mutate(AI=factor(isoAI | geneAI, levels=c("TRUE","FALSE"))) %>% 
  ggplot(aes(x=pvalue, fill=AI)) + geom_histogram() +
  ggtitle("WASP gene-level p-values")

wasp_qval <- wasp %>% select(gene_id, pvalue, .drop_ranges=TRUE) %>%
  as.data.frame() %>% tibble()

###

  distinct(gene_id, pvalue) %>%
  mutate(qvalue=p.adjust(pvalue, method="BH"))
wasp_qval <- wasp_qval %>% filter(!duplicated(gene_id)) %>%
  inner_join(dat) %>% select(qvalue, txp)

# need to run eval.R script to build 'padj'

padj$wasp <- 1
padj[wasp_qval$txp,"wasp"] <- wasp_qval$qvalue

cpnf <- calculate_performance(cd,
                              binary_truth="status",
                              aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                       "fdrtprcurve","overlap"),
                              thrs=c(.01,.05,.1),
                              thr_venn=.05)

cplotnf <- prepare_data_for_plot(cpnf, colorscheme=cols)

# simple plot
plot_tpr(cplot, pointsize=2.5)

# FDR and number
plot_fdrnbrcurve(cplotnf, xaxisrange=xrng) +
  ggplot2::ylim(0,4500)

# upset plot
plot_upset(cplotnf, order.by="freq", nintersects=10)

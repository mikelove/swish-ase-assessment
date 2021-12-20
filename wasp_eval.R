cols <- palette.colors()[1:6]
types <- c("txp","gene","tss","oracle")
names(cols)[1:4] <- types
names(cols)[5] <- "truth"
names(cols)[6] <- "wasp"

library(dplyr)
library(ggplot2)
wasp <- read.table("../ase-sim/wasp_cht/cht_results.txt", header=TRUE)
wasp <- wasp %>% mutate(ratio=ALT.AS.READ.COUNT/REF.AS.READ.COUNT,
                        total=TOTAL.AS.READ.COUNT,
                        pvalue=P.VALUE,
                        chr_loc=paste0(TEST.SNP.CHROM,"-",TEST.SNP.POS)) %>%
  tibble()
wasp %>% filter(total > 50 & total < 1e5) %>%
  ggplot(aes(total, log2(ratio), col=-log10(pvalue+1e-15))) +
  geom_point() + ylim(-.5, .5)
library(GenomicRanges)
load("../ase-sim/granges.rda")
dat <- data.frame(chr=seqnames(txps),
                  gene_id=txps$gene_id,
                  txp=names(txps),
                  loc=sapply(txps$snp_loc, min))
dat <- dat %>% mutate(chr_loc=paste0(chr, "-", loc))
wasp_qval <- wasp %>% inner_join(dat, by="chr_loc") %>%
  distinct(gene_id, pvalue) %>%
  mutate(qvalue=p.adjust(pvalue, method="BH"))
# fix this...
wasp_qval <- wasp_qval %>% filter(!duplicated(gene_id)) %>%
  inner_join(dat) %>% select(qvalue, txp)
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

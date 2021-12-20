library(iCOBRA)

cols <- palette.colors()[1:6]
types <- c("txp","gene","tss","oracle")
names(cols)[1:4] <- types
names(cols)[5] <- "truth"
names(cols)[6] <- "wasp"

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

wasp <- TRUE
if (wasp) {
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
}

cd <- COBRAData(padj=padj, truth=truth)
cp <- calculate_performance(cd,
                            binary_truth="status",
                            aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                     "fdrtprcurve","overlap"),
                            splv="group",
                            thrs=c(.01,.05,.1),
                            thr_venn=.05)
cpnf <- calculate_performance(cd,
                              binary_truth="status",
                              aspect=c("tpr","fdr","fdrtpr","fdrnbr",
                                       "fdrtprcurve","overlap"),
                              thrs=c(.01,.05,.1),
                              thr_venn=.05)

cplot <- prepare_data_for_plot(cp, colorscheme=cols)
cplotnf <- prepare_data_for_plot(cpnf, colorscheme=cols)

# simple plot
plot_tpr(cplot, pointsize=2.5)

# TPR over FDR
yrng <- c(0,1)
#xrng <- c(0,.15)
xrng <- c(0,.5)
plot_fdrtprcurve(cplotnf,
                 plottype="points",
                 xaxisrange=xrng,
                 yaxisrange=yrng,
                 title="txp-level AI testing")

if (FALSE) {
  # FDR and number
  plot_fdrnbrcurve(cplotnf, xaxisrange=xrng) +
    ggplot2::ylim(0,4500)
  
  # upset plot
  plot_upset(cplotnf, order.by="freq", nintersects=10)
}

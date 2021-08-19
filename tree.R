dir <- "/proj/milovelab/wu/bulk-ase/ase-analysis/exp4_s10/tree_analysis/"
dir0 <- "/proj/milovelab/wu/bulk-ase/ase-analysis/exp4_s10/"
res <- read.csv(file.path(dir,"tree_swish_res.csv"))
cov <- read.csv(file.path(dir,"tree_infRV.csv"))
load(file.path(dir0,"tree_dat.rda"))

library(dplyr)
dat$gene <- sub(".*_","",dat$txp)
dat$abund <- sub(".*\\|.*_(.*)_FBgn.*","\\1",dat$txp)
dat %>% filter(Group_4 != "unknown") %>%
  summarize(uniq = n_distinct(gene))
length(unique(dat$gene))

dat %>% filter(Group_4 != "unknown") %>%
  group_by(gene) %>%
  summarize(ai = any(abund != 2)) %>%
  summarize(sum(ai))

gene <- dat[dat$Group_2 == "69693",]
txps <- unique(sub("_.*","",gene$txp))


res %>% filter(X %in% c(txps, unique(gene$Group_3))) %>%
  select(X, log2FC, pvalue)


library(tibble)
library(dplyr)
library(tidyr)

tss <- read.delim("../res/tss.tsv") %>%
  rownames_to_column("tss_groups") %>%
  tibble()

txp <- read.delim("../res/txp.tsv") %>%
  rownames_to_column("tx_id") %>%
  tibble()

load("../granges_with_groups.rda")
suppressPackageStartupMessages(library(GenomicRanges))
tab <- mcols(txps)[,c("tx_id","gene_id","tss_groups")] %>%
  as.data.frame() %>%
  tibble()

tss_thin <- tss %>% select(tss_groups, tss_q = qvalue)

dat <- txp %>%
  inner_join(tab) %>%
  inner_join(tss_thin) %>%
  group_by(tss_groups) %>%
  filter(pvalue == min(pvalue))
         
library(ggplot2)
library(ggrepel)

dat2 <- dat %>%
    filter(tss_q < 1e-3 & pvalue > .1)

dat %>%
    ggplot(aes(-log10(qvalue), -log10(tss_q), color=log10mean)) +
    geom_point() +
    geom_label_repel(data=dat2, aes(-log10(qvalue), -log10(tss_q), label=tss_groups))



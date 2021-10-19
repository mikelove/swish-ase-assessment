dir <- "/proj/milovelab/wu/bulk-ase/ase-analysis/exp4_s10/tree_analysis/"
dir0 <- "/proj/milovelab/wu/bulk-ase/ase-analysis/exp4_s10/"
res <- read.csv(file.path(dir,"tree_swish_res.csv"))
cov <- read.csv(file.path(dir,"tree_infRV.csv"))
load(file.path(dir0,"tree_dat.rda"))

head(dat)
head(cov)
head(res)

library(dplyr)
dat <- dat[grep("_M", dat$txp),]
#dat$gene <- sub(".*_","",dat$txp)
#dat$txp0 <- sub("_M\\|.*","",dat$txp)
#dat$abund <- sub(".*\\|.*_(.*)_FBgn.*","\\1",dat$txp)
names(dat)[12:14] <- c("txp0","gene","abund")

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

###

# Noor's question:
# what is reduction in inferential uncertainty for all the genes with AI
# that dont have groups and see where they lie w.r.t threshold

dat_m <- merge(dat, cov[,1:2], by.x="txp0", by.y="X") %>%
  arrange(gene, abund)

ai_genes <- dat_m %>% group_by(gene) %>%
  summarize(ai = any(abund != 2)) %>%
  filter(ai) %>% pull(gene)
length(ai_genes) # 1000

dat_ai <- dat_m %>% filter(gene %in% ai_genes) %>%
  mutate(sub_gene_agg = Group_4 != "unknown")

dat_ai %>% ggplot(aes(sub_gene_agg, meanInfRV)) +
  geom_boxplot() +
  scale_y_log10() +
  ggtitle("Transcripts in AI genes")

dat_ai %>% group_by(gene) %>%
  summarize(geneInfRV = mean(meanInfRV), sub_gene_agg=any(sub_gene_agg)) %>%
  ggplot(aes(sub_gene_agg, geneInfRV)) +
  geom_boxplot() +
  scale_y_log10() +
  ggtitle("Per gene (mean InfRV over txps) for AI genes")

gene_ai <- dat_ai %>% group_by(gene) %>%
  summarize(geneInfRV = mean(meanInfRV), sub_gene_agg=any(sub_gene_agg))

agg_genes <- dat %>% filter(gene %in% ai_genes) %>%
  arrange(gene, abund) %>%
  group_by(gene) %>%
  summarize(gene_agg = all(Group_2 == Group_2[1])) %>%
  filter(gene_agg) %>% pull(gene)

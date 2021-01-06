load("txi.rda")
dim(txi$counts)
library(DESeq2)
m <- nrow(txi$counts)/2
n <- 6 # total count of samples
cts <- cbind(txi$counts[1:m,], txi$counts[(m+1):nrow(txi$counts),])
cts <- round(cts)
coldata <- data.frame(condition=factor(rep(c("M","P"),each=n), levels=c("M","P")),
                      id=factor(rep(1:n,2)))
dds <- DESeqDataSetFromMatrix(cts, coldata, ~0 + id + condition)

keep <- rowSums(counts(dds) >= 10) >= 6
table(keep) # this will look better when lib size ~5e6

dds <- dds[keep,]
dds <- dds[1:1000,] # for speed

dds <- DESeq(dds, fitType="mean")

plotDispEsts(dds)
abline(h=1/100, col="purple")
dispersionFunction(dds)
table(mcols(dds)$dispGeneEst < .01)

plot(mcols(dds)$baseMean, mcols(dds)$dispGeneEst, log="xy", pch=20)
abline(h=1/100, col="purple")
identify(mcols(dds)$baseMean, mcols(dds)$dispGeneEst)

idx <- 358
dat <- plotCounts(dds, gene=idx, intgroup=c("condition","id"), returnData=TRUE)
library(ggplot2)
ggplot(dat, aes(condition, count, group=id)) + geom_point() + geom_line()
rownames(dds)[idx]

txpNm <- "FBtr0305084" # suspicious txp
(txpNms <- grep(txpNm, rownames(txi$counts), value=TRUE))
susp_txps <- grep(txpNm, rownames(txi$counts))
round(txi$variance[susp_txps,] / (txi$counts[susp_txps,] + 1), 1)

load("sim_counts_matrix.rda")
dat2 <- dat
dat2$count <- as.vector(t(counts_matrix[txpNms,]))
dat <- rbind(dat, dat2)
dat$type <- factor(rep(c("salmon","true"),each=12))
dat$id <- as.numeric(dat$id)
dat$id[13:24] <- rep(7:12, 2)
dat$id <- factor(dat$id)
ggplot(dat, aes(condition, count, group=id, col=type)) +
  geom_point() +
  geom_line()


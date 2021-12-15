table(mcols(y)$qvalue < .01)
hist(mcols(y)$pvalue, col="grey")

# oracle
alpha <- .01
sig_txps <- rownames(y)[mcols(y)$qvalue < alpha]
table(grepl("-2$", sig_txps))
table(ai=!grepl("-2$", rownames(y)), sig=mcols(y)$qvalue < alpha)
prop.table(table(ai=!grepl("-2$", rownames(y)), sig=mcols(y)$qvalue < alpha), 2)

# txp-level
mcols(y)$abund <- mcols(txps)[rownames(y),"abundance"]

sig <- mcols(y)$qvalue < 0.1
ai <- mcols(y)$abund != 2
table(ai, sig)

suppressPackageStartupMessages(library(SummarizedExperiment))
library(fishpond)
library(apeglm)
library(locfdr)
types <- c("txp","oracle","gene","tss")
for (t in types) {
  load(file=paste0("data/se_",t,".rda"))
  y <- wide; rm(wide)
  assays(y) <- assays(y)["counts"]
  y <- labelKeep(y)
  y <- y[mcols(y)$keep,]
  n <- ncol(y)/2
  ase_cts <- round(assay(y)[,(n+1):(2*n)])
  tot_cts <- round(assay(y)[,1:n] + ase_cts)
  idx <- tot_cts == 0
  tot_cts[idx] <- 1
  set.seed(1)
  ase_cts[idx] <- rbinom(sum(idx), 1, .5)
  theta.hat <- 10
  x <- matrix(rep(1,n), ncol=1)
  param <- cbind(theta.hat, tot_cts)
  mle <- rowMeans((ase_cts+1) / (tot_cts+2))
  beta_hat <- log(mle / (1 - mle))
  theta_hat <- bbEstDisp(success=ase_cts, size=tot_cts, x=x,
                         beta=beta_hat, minDisp=1, maxDisp=500)
  param <- cbind(theta_hat, tot_cts)
  fit <- apeglm(Y=ase_cts, x=x, log.lik=NULL, param=param,
                no.shrink=TRUE, log.link=FALSE, method="betabinCR")
  z <- fit$map[,1]/fit$sd[,1]
  fit <- locfdr(z, df=n-1)
  df <- data.frame(qvalue=fit$fdr, row.names=rownames(y))
  write.table(df, file=paste0("res/bb_",t,".tsv"), sep="\t",quote=FALSE)
}

suppressPackageStartupMessages(library(SummarizedExperiment))
library(fishpond)
library(apeglm)
library(locfdr)
types <- c("txp","oracle","gene","tss")
for (t in types) {
  load(file=paste0("data/se_",t,".rda"))
  y <- wide
  assays(y) <- assays(y)["counts"]
  rm(wide)
  y <- labelKeep(y)
  y <- y[mcols(y)$keep,]
  n <- ncol(y)/2
  cts <- assay(y)
  
  theta.hat <- 10
  x <- matrix(rep(1,n),ncol=1)
  param <- cbind(theta.hat, cts)
  mle <- rowMeans(ase.cts/cts)
  beta.hat <- log(mle / (1 - mle))
  theta.hat <- bbEstDisp(success=ase.cts, size=cts,
                         x=x, beta=beta.hat,
                         minDisp=1, maxDisp=50)
  param <- cbind(theta.hat, cts)
  fit <- apeglm(Y=ase.cts, x=x, log.lik=NULL, param=param,
                no.shrink=TRUE, log.link=FALSE, method="betabinCR")
  mle <- cbind(fit$map[,1], fit$sd[,1])
  z <- mle[,1]/mle[,2]
  fit <- locfdr(z)
  write.table(..., file=paste0("res/bb_",t,".tsv"), sep="\t",quote=FALSE)
}

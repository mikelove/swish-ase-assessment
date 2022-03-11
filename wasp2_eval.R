library(readr)
library(dplyr)
dir <- "../ase-sim/wasp2_results"
samps <- list.files(dir)
files <- list.files(file.path(dir, samps[1]))
x <- read_tsv(file.path(dir, samps[1], files[1]))
hist(x$pval)
x <- read_tsv(file.path(dir, samps[1], files[2]))
hist(x$pval)

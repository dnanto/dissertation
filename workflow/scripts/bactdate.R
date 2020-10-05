#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(lubridate)
  library(BactDating)
})

path.tree <- snakemake@input[["tree"]]
path.meta <- snakemake@input[["meta"]]

meta <- read_tsv(path.meta, col_types = cols(.default = "c"))

pre <- file.path(dirname(path.tree), str_split_fixed(basename(path.tree), "\\.", 2)[, 1])
gub <- loadGubbins(pre)

tip.date <- decimal_date(ymd(with(meta, date[match(gub$tip.label, accver)])))

p <- snakemake@params
set.seed(p$seed)
res <- bactdate(gub, tip.date, model = p$mod, nbIts = p$nbi, thin = p$thn, useRec = T)

qs::qsave(res, snakemake@output[["res"]])

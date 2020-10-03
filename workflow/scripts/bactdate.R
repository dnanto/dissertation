#!/usr/bin/env RScript --vanilla

library(stringr, warn.conflicts = F, quietly = T, verbose = F)
library(lubridate, warn.conflicts = F, quietly = T, verbose = F)
library(BactDating, warn.conflicts = F, quietly = T, verbose = F)


path <- snakemake@input[[1]]
pre <- file.path(dirname(path), str_split_fixed(basename(path), "\\.", 2)[, 1])

gub <- loadGubbins(pre)

tip.date <- decimal_date(ymd(str_match(gub$tip.label, "(\\d{4}-\\d{2}-\\d{2})$")[, 2]))

res <- with(
  snakemake@params,
  bactdate(gub, tip.date, model = mod, nbIts = nbi, thin = thn, useRec = T)
)

qs::qsave(res, snakemake@output[[1]])

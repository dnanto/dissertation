#!/usr/bin/env RScript --vanilla

rmarkdown::render(
  snakemake@params$rmd,
  rmarkdown::html_fragment(),
  snakemake@output[[1]],
  params = tail(snakemake@input, -4)
)

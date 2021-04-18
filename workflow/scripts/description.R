#!/usr/bin/env RScript --vanilla

rmarkdown::render(
  snakemake@params$rmd,
  rmarkdown::html_fragment(),
  snakemake@output[[1]],
  params = list(
    tsv = snakemake@input[[1]],
    gff = snakemake@input[[2]],
    ess = snakemake@params$ess,
    typ = snakemake@params$typ
  ),
  knit_root_dir = snakemake@params$root
)

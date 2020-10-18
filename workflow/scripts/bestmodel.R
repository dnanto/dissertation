#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  library(qs)
  library(coda)
  library(rjson)
  library(dplyr)
  library(readr)
  library(stringr)
  library(treeio)
  library(BactDating)
})


node_json <- function(res, file) {
  root <- min(res$tree$edge[, 1])
  iden <- c(root, res$tree$edge[, 2])
  label <- c(res$tree$tip.label, res$tree$node.label)[iden]
  data <- data.frame(
    branch_length = c(0, res$tree$edge.length),
    clock_length = c(0, res$tree$edge.length),
    mutation_length = c(0, res$tree$subs),
    num_date = c(leafDates(res$tree), nodeDates(res$tree))[iden],
    num_date_confidence = I(res$CI[iden, ])
  )
  toJSON(list(nodes = split(data, label)), indent = 1)
}

path <- do.call(c, snakemake@input["models"])
qlen <- ncol(ape::read.dna(snakemake@input$qry, format = "fasta", as.character = T))
ncpu <- snakemake@params[[1]]

data.frame(
  path = path,
  species = basename(dirname(dirname(path))),
  gene = basename(dirname(path)),
  model = str_split(basename(tools::file_path_sans_ext(path)), "-", n = 2, simplify = T)[, 2],
  bind_rows(
    lapply(path, function(ele) {
      res <- qs::qread(ele, nthreads = ncpu)
      # 50% burn-in
      record <- res$record
      record <- with(res, {
        record[max(1, round(nrow(record) * 0.5)):nrow(record), c("mu", "sigma", "alpha")]
      })
      med <- setNames(apply(record, 2, median), c("median.mu", "median.sigma", "median.alpha"))
      ess <- setNames(effectiveSize(as.mcmc(record)), c("ess.mu", "ess.sigma", "ess.alpha"))
      c(dic = res$dic, ssy = med[1] / qlen, med, ess)
    })
  ),
  stringsAsFactors = F
) %>%
  write_tsv(snakemake@output[["tsv"]]) %>%
  group_by(species, gene) %>%
  group_walk(function(x, y) {
    path <- x[which.min(x$dic), ]$path
    res <- qread(path, nthreads = ncpu)
    write(node_json(res), snakemake@output[["json"]])
    write.tree(res$tree, snakemake@output[["tree"]])
  })

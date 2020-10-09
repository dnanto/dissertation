#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  library(qs)
  library(rjson)
  library(dplyr)
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

path <- do.call(c, snakemake@input)
ncpu <- snakemake@params[[1]]

data.frame(
  path = path,
  species = basename(dirname(dirname(path))),
  gene = basename(dirname(path)),
  model = str_split(basename(tools::file_path_sans_ext(path)), "-", n = 2, simplify = T)[, 2],
  dic = sapply(path, function(ele) qread(ele, nthreads = ncpu)$dic),
  stringsAsFactors = F
) %>%
  group_by(species, gene) %>%
  group_walk(function(x, y) {
    path <- x[which.min(x$dic), ]$path
    res <- qread(path, nthreads = ncpu)
    write(node_json(res), snakemake@output[["json"]])
    write.tree(res$tree, snakemake@output[["tree"]])
  })

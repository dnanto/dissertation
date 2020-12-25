#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
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
  rjson::toJSON(list(nodes = split(data, label)), indent = 1)
}

path.in.meta <- snakemake@input[["meta"]]
path.in.tree <- snakemake@input[["tree"]]
path.out.tsv <- snakemake@output[["tsv"]]
path.out.tree <- snakemake@output[["tree"]]
path.out.json <- snakemake@output[["json"]]

# params <- list(root = dirname(path.in.tree), nbIts = 10000, thin = 10, burn = 0.5, threads = 16)
params <- snakemake@params

threads <- snakemake@threads

set.seed(params$seed)

meta <- read_tsv(path.in.meta, col_types = cols(.default = "c"))

# ##gff-version 3
# ##sequence-region SEQUENCE 1 34892
len <- (
  file.path(dirname(path.in.tree), "gub.recombination_predictions.gff") %>%
    readLines(n = 2) %>%
    head(n = 2) %>%
    last() %>%
    str_split_fixed(pattern = " ", 4) %>%
    last() %>%
    as.integer()
)

pre <- file.path(dirname(path.in.tree), str_split_fixed(basename(path.in.tree), "\\.", 2)[, 1])
phy <- loadGubbins(pre)
phy$tip.date <- decimal_date(ymd(with(meta, date[match(phy$tip.label, accver)])))
phy <- initRoot(phy, phy$tip.date, useRec = T)

path <-(
parallel::mclapply(params$model, function(ele) {
  path <- file.path(params$root, paste0("bac-", ele, ".qs"))
  run <- bactdate(phy, phy$tip.date, nbIts = params$nbIts, thin = params$thin, model = ele, useRec = T)
  qs::qsave(run, path)
  ape::write.tree(run$tree, file.path(params$root, paste0("bac-", ele, ".tree")))
  path
}, mc.cores = threads) %>%
  lapply(function(path) {
    run <- qs::qread(path, nthreads = threads)
    key <- c("likelihood", "mu", "sigma", "alpha", "prior")
    rec <- with(run, record[max(1, round(nrow(record) * params$burn)):nrow(record), key])
    est <-
      apply(rec[, key], 2, summary) %>%
      t() %>%
      apply(1, function(row) {
        sprintf("%.2e [%.2e;%.2e]", row["Median"], row["1st Qu."], row["3rd Qu."])
      }) %>%
      setNames(paste0("est.", key))
    ess <- setNames(coda::effectiveSize(rec), paste0("ess.", key))
    mod <- (
      basename(path) %>%
        tools::file_path_sans_ext() %>%
        str_remove("^bac-")
    )
    c(
      mod = mod, path = path,
      dic = run$dic, rootprob = run$rootprob,
      rate = median(rec[, "mu"]) / len,
      est, ess
    )
  }) %>%
  bind_rows() %>%
  mutate(
    across(c("dic", "rootprob", "rate"), as.numeric),
    across(starts_with("ess."), as.numeric)
  ) %>%
  arrange(dic) %>%
  select(path, everything()) %>%
  write_tsv(path.out.tsv) %>%
  slice_min(n = 1, order_by = dic) %>%
  pull(path)
)

file.copy(str_replace(path, ".qs$", ".tree"), path.out.tree)
qs::qread(path, nthreads = threads) %>%
  node_json() %>%
  write(file = path.out.json)

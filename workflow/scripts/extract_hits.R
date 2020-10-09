#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

parse_8C <- function(path) {
  col_types <- cols(
    `% identity` = "d", `alignment length` = "i",
    mismatches = "i", `gap opens` = "i",
    `q. start` = "i", `q. end` = "i",
    `s. start` = "i", `s. end` = "i",
    evalue = "d", `bit score` = "d",
    .default = "c"
  )
  header <- read_lines(path, n_max = 6)
  qlen <- as.integer(tail(strsplit(header[3], " ")[[1]], 2)[1])
  fields <- strsplit(substr(header[5], 11, nchar(header[5])), ", ")[[1]]
  cbind(
    read_tsv(path, col_names = fields, col_types = col_types, comment = "#"),
    `query length` = qlen
  )
}

df.meta <- read_tsv(
  snakemake@input[["meta"]],
  col_types = cols(.default = "c"),
  na = "?"
)

parse_8C(snakemake@input[["hits"]]) %>%
  bind_rows() %>%
  mutate(
    qpidcov = abs(`q. end` - `q. start`) / `query length` * `% identity`,
    strand = if_else(`q. start` < `q. end`, "plus", "minus")
  ) %>%
  group_by(`subject id`) %>%
  slice_max(qpidcov) %>%
  ungroup() %>%
  left_join(df.meta, by = c(`subject id` = "accver")) ->
df.qpidcov

with(
  df.qpidcov,
  {
    writeLines(
      paste0(strand, "\t", `subject id`, ":", `s. start`, "-", `s. end`),
      snakemake@output[["tsv"]]
    )
  }
)

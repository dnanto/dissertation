#!/usr/bin/env RScript --vanilla

library(readr, warn.conflicts = F, quietly = T, verbose = F)
library(dplyr, warn.conflicts = F, quietly = T, verbose = F)


args <- commandArgs(trailingOnly = T)
f <- if (args[1] == "-") file("stdin") else args[1] # file
d <- args[2] # delimiter
c1 <- args[3] # column to convert to
f1 <- args[4] # format to convert to
c2 <- args[5] # column to convert from
f2 <- args[6:length(args)] # formats to convert from


df <- read_delim(f, d, col_types = cols(.default = "c"))

lubridate::parse_date_time(df[[c2]], f2, quiet = T) %>%
  strftime(f1) %>%
  as.data.frame() %>%
  setNames(c1) %>%
  mutate_all(coalesce, "?") %>%
  bind_cols(df, .) %>%
  format_delim(d) %>%
  cat()

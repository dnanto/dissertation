#!/usr/bin/env RScript --vanilla

library(readr, warn.conflicts = F, quietly = T, verbose = F)
library(dplyr, warn.conflicts = F, quietly = T, verbose = F)
library(tidyr, warn.conflicts = F, quietly = T, verbose = F)
library(tibble, warn.conflicts = F, quietly = T, verbose = F)
library(stringr, warn.conflicts = F, quietly = T, verbose = F)


args <- commandArgs(trailingOnly = T)
f <- if (args[1] == "-") file("stdin") else args[1] # file
d <- args[2] # delimiter
k <- args[3] # column
p <- args[4] # pattern
n <- if (length(args) <= 4) parallel::detectCores() else args[5] # cores


df <- read_delim(f, d, col_types = cols(.default = "c"))
country <- unique(df[[k]])

# subset and order by feature code priority
codes <- c("country_code", paste0("admin", 1:4, "_code"))
levels <- as.vector(do.call(c, neogeonames::akfc))
geoname <- select(neogeonames::geoname, feature_code, all_of(codes), latitude, longitude)
geoname <- geoname[order(ordered(geoname$feature_code, levels = levels)), ]

# coordinates
coor <- (
  parallel::mclapply(country, neogeonames::adminify_delim, delim = "[:,]", mc.cores = n) %>%
    lapply(`[[`, "ac") %>%
    bind_rows() %>%
    setNames(codes) %>%
    apply(1, function(row) {
      Filter(Negate(is.na), row) %>%
        t() %>%
        as.data.frame() %>%
        left_join(geoname, by = intersect(names(.), codes)) %>%
        slice_head(n = 1)
    }) %>%
    bind_rows()
)

# location
select(coor, all_of(codes)) %>%
  mutate(across(.fns = ~ replace_na(., "?"))) %>%
  apply(1, paste, collapse = ".") %>%
  str_remove("(\\.\\?)+$") %>%
  enframe(name = NULL, value = "location") %>%
  bind_cols(coor) %>%
  mutate(country = country) %>%
  left_join(df, ., by = "country") %>%
  select(all_of(names(df)), location, latitude, longitude, feature_code, all_of(codes)) %>%
  # output
  mutate_all(as.character) %>%
  mutate_all(coalesce, "?") %>%
  format_delim(d) %>%
  cat()

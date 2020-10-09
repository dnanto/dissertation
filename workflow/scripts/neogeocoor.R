#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

codes <- c("country_code", paste0("admin", 1:4, "_code"))

df.meta <- read_tsv(
  snakemake@input[[1]],
  col_types = cols(.default = "c"),
  na = ""
)

# subset and order by feature code priority
levels <- as.vector(do.call(c, neogeonames::akfc))
geoname <- select(
  neogeonames::geoname,
  feature_code, codes, latitude, longitude
)
geoname <- geoname[order(ordered(geoname$feature_code, levels = levels)), ]

separate(df.meta, "location", codes, "\\.", extra = "drop", fill = "right", remove = F) %>%
  select(location, all_of(codes)) %>%
  distinct() %>%
  apply(1, function(row) {
    Filter(Negate(is.na), row) %>%
      .[. != "?"] %>%
      t() %>%
      data.frame() %>%
      left_join(geoname, by = tail(names(.), -1)) %>%
      slice_head(n = 1)
  }) %>%
  bind_rows() %>%
  mutate(type = "location") %>%
  select(type, location, latitude, longitude) %>%
  write_tsv(snakemake@output[[1]], col_names = F)

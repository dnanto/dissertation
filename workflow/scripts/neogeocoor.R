#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

df.meta <- read_tsv(
  snakemake@input[[1]],
  col_types = cols(.default = "c"),
  na = "?"
)

levels <- as.vector(do.call(c, neogeonames::akfc))
geoname <- select(
  neogeonames::geoname,
  feature_code, country_code, paste0("admin", 1:4, "_code"), latitude, longitude
)
geoname <- geoname[order(ordered(geoname$feature_code, levels = levels)), ]

codes1 <- c("country", paste0("div", 1:4))
codes2 <- c("country_code", paste0("admin", 1:4, "_code"))

lapply(1:5, function(idx) {
  bind_rows(
    select(df.meta, codes1[1:idx]) %>%
      filter(!is.na(.[idx])) %>%
      distinct() %>%
      apply(1, function(row) {
        if (!all(row==F)) {
          df <- data.frame(t(row[(i <- which(!is.na(row)))]))
          left_join(setNames(df, codes2[i]), geoname, by = codes2[i])[1, ] %>%
            replace(is.na(.), "?") %>%
            mutate(
              type = codes1[idx],
              name = apply(.[codes2[1:idx]], 1, paste0, collapse = ".")
            )
        }
      })
  )
}) %>%
  bind_rows() %>%
  select(type, name, latitude, longitude) %>%
  write_tsv(snakemake@output[[1]], col_names = F)

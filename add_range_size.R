# add_range_size.R --- add each species' range size (number of 0.25-degree grid
# cells in its range, the denominator of every "fraction of range exposed") to
# the app's species metadata, from the pipeline's species attributes table.
# Run after regenerating app data:  Rscript add_range_size.R
suppressPackageStartupMessages({library(qs2); library(dplyr)})
ext  <- file.path("inst", "extdata")
attr_file <- "/Users/cory.merow/Documents/SDMs/Exposure_2025/V8/metadata_V8/spAttributes_v8.qs"
meta <- qs_read(file.path(ext, "allExpForShiny_species.qs"))
rs <- suppressWarnings(qs_read(attr_file)) %>% distinct(spName, rangeSize) %>% filter(!is.na(rangeSize))
stopifnot(!any(duplicated(rs$spName)))
meta$rangeSize <- NULL
meta <- meta %>% left_join(rs, by = "spName")
cat(sprintf("species: %d, with rangeSize: %d, missing: %d\n", nrow(meta), sum(!is.na(meta$rangeSize)), sum(is.na(meta$rangeSize))))
qs_save(meta, file.path(ext, "allExpForShiny_species.qs"))
cat("written: inst/extdata/allExpForShiny_species.qs (rangeSize = range cells on the app grid)\n")

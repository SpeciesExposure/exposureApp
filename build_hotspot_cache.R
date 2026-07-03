#!/usr/bin/env Rscript
# Build precomputed hotspot raster data for recent years.
# Run once (or after each pipeline run that updates year shards).
# Output: exposureApp/inst/extdata/allExpForShiny_hotspot_cache_v1.qs
#
# Usage:
#   Rscript exposureApp/build_hotspot_cache.R
#   # or to rebuild specific years:
#   Rscript exposureApp/build_hotspot_cache.R 2020 2021 2022 2023 2024 2025

suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(terra)
})

args <- commandArgs(trailingOnly = TRUE)

# locate extdata relative to this script or cwd
script_dir <- tryCatch(
  dirname(normalizePath(commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))][[1]]
                        |> sub("--file=", "", x = .))),
  error = function(e) "."
)
EXTDATA <- file.path(script_dir, "inst/extdata")
if (!dir.exists(EXTDATA)) EXTDATA <- "exposureApp/inst/extdata"
if (!dir.exists(EXTDATA)) EXTDATA <- file.path(getwd(), "exposureApp/inst/extdata")
if (!dir.exists(EXTDATA))
  stop("Cannot locate exposureApp/inst/extdata. Run from project root or exposureApp/.")

MANIFEST   <- file.path(EXTDATA, "allExpForShiny_manifest_v1.qs")
TEMPLATE   <- file.path(EXTDATA, "landTemplate.tif")
OUT_CACHE  <- file.path(EXTDATA, "allExpForShiny_hotspot_cache_v1.qs")

stopifnot(file.exists(MANIFEST), file.exists(TEMPLATE))

m   <- qs2::qs_read(MANIFEST)
tpl <- terra::rast(TEMPLATE)
n_cells <- terra::ncell(tpl)

year_paths <- setNames(
  vapply(names(m$year_files),
         function(yr) normalizePath(file.path(EXTDATA, m$year_files[[yr]]), mustWork = FALSE),
         character(1)),
  names(m$year_files)
)
all_years <- as.integer(names(year_paths))
vars_avail <- m$vars_avail
vars_avail <- vars_avail[!is.na(vars_avail) & nzchar(vars_avail)]
# The pipeline manifest leaves vars_avail empty; the app itself falls back to
# the variables present in the data. Do the same here so the hotspot counts are
# not silently filtered to zero rows.
if (!length(vars_avail)) {
  tc_path <- file.path(EXTDATA, "allExpForShiny_trends_cache_v1.qs")
  if (file.exists(tc_path)) {
    tc <- tryCatch(qs2::qs_read(tc_path), error = function(e) NULL)
    if (!is.null(tc) && !is.null(tc$sp_trend_df))
      vars_avail <- sort(unique(as.character(tc$sp_trend_df$var)))
  }
  if (!length(vars_avail)) {
    fp1 <- normalizePath(file.path(EXTDATA, m$year_files[[length(m$year_files)]]), mustWork = FALSE)
    if (file.exists(fp1)) vars_avail <- sort(unique(as.character(qs2::qs_read(fp1)$var)))
  }
  cat(sprintf("manifest vars_avail empty; using %d vars from data\n", length(vars_avail)))
}

# years to precompute: use args if given, else most recent 5 + last 3
if (length(args)) {
  target_years <- as.integer(args)
} else {
  target_years <- sort(unique(c(
    tail(sort(all_years), 5),       # 5 most recent
    2023L, 2024L, 2025L             # always include these if present
  )))
}
target_years <- target_years[target_years %in% all_years]
cat(sprintf("Building hotspot cache for %d years: %s\n",
            length(target_years), paste(target_years, collapse = ", ")))

# load existing cache so we can update incrementally
existing <- if (file.exists(OUT_CACHE)) {
  tryCatch(qs2::qs_read(OUT_CACHE), error = function(e) list())
} else list()

cache <- existing
for (yr in target_years) {
  key <- as.character(yr)
  fp  <- year_paths[key]
  if (!file.exists(fp)) { cat(sprintf("  skip %d (shard missing)\n", yr)); next }

  cat(sprintf("  computing year %d ...", yr))
  d <- qs2::qs_read(fp)

  counts <- d %>%
    filter(!is.na(cell), cell >= 1L, cell <= n_cells,
           var %in% vars_avail) %>%
    group_by(cell) %>%
    summarize(n = n_distinct(spName), .groups = "drop")

  if (!nrow(counts)) {
    cat(" (no rows)\n")
    next
  }

  # store as compact integer vector (NA outside exposed cells)
  cell_vals            <- rep(NA_integer_, n_cells)
  cell_vals[counts$cell] <- as.integer(counts$n)

  cache[[key]] <- list(
    year       = yr,
    cell_vals  = cell_vals,   # integer vector length ncell
    minN       = min(counts$n),
    maxN       = max(counts$n),
    n_cells    = n_cells,
    vars_used  = vars_avail,
    built_at   = Sys.time()
  )
  cat(sprintf(" %d exposed cells, max=%d\n", nrow(counts), max(counts$n)))
}

qs2::qs_save(cache, OUT_CACHE)
cat(sprintf("Saved hotspot cache (%d years) -> %s\n", length(cache), OUT_CACHE))

suppressPackageStartupMessages(suppressWarnings({
  library(shiny)
  library(qs2)
  library(dplyr)
  library(terra)
  library(leaflet)
  library(leaflet.extras)
  library(jsonlite)
  library(DT)
  library(raster)   # needed for addRasterImage
}))

`%||%` <- function(a, b) if (!is.null(a)) a else b
APP_VERBOSE <- identical(tolower(Sys.getenv("EXPOSURE_APP_VERBOSE", "false")), "true")
log_msg <- function(...) if (isTRUE(APP_VERBOSE)) message(...)

read_qs <- function(path) {
  qs2::qs_read(path)
}

write_qs <- function(object, path) {
  qs2::qs_save(object, path)
}

# ── 1. Locate Int directory ──────────────────────────────────────────────────
find_int_dir <- function() {
  env <- Sys.getenv("INT_DIR", "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env))
  # Search one and two levels up from the app directory
  for (parent in c(normalizePath(file.path(getwd(), "..")),
                   normalizePath(file.path(getwd(), "../..")))) {
    cand <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
    hit  <- grepl("Int_", basename(cand))
    if (any(hit)) return(cand[which(hit)[1]])
  }
  stop("Cannot find Int_* directory. Set INT_DIR env var before launching.")
}

intDir <- tryCatch(find_int_dir(), error = function(e) stop(e$message))
allcells_path <- file.path(intDir, "allExpForShiny.qs")
template_path <- file.path(intDir, "landTemplate.tif")
manifest_path <- file.path(intDir, "allExpForShiny_manifest_v1.qs")
species_meta_path <- file.path(intDir, "allExpForShiny_species.qs")
trends_cache_path <- file.path(intDir, "allExpForShiny_trends_cache_v1.qs")
species_year_cells_path  <- file.path(intDir, "allExpForShiny_species_year_cells_v1.qs")
range_cells_mb_path      <- file.path(intDir, "range_cells_mb.qs")
hotspot_cache_path       <- file.path(intDir, "allExpForShiny_hotspot_cache_v1.qs")
stopifnot(file.exists(template_path))

log_msg("Loading template raster ...")
tpl    <- rast(template_path)

to_title_case <- function(x) {
  x <- as.character(x)
  x <- trimws(tolower(x))
  x <- gsub("(^|[[:space:]-])([[:alpha:]])", "\\1\\U\\2", x, perl = TRUE)
  x[nchar(x) == 0] <- NA_character_
  x
}

empty_app_df <- function() {
  data.frame(
    spName = character(),
    cell = integer(),
    var = character(),
    year = integer(),
    group = character(),
    orderName = character(),
    familyName = character(),
    redlistCategory = character(),
    propExposed = numeric(),
    stringsAsFactors = FALSE
  )
}

normalize_app_df <- function(d) {
  if (is.null(d) || !nrow(d)) return(empty_app_df())
  if (!("spName" %in% names(d))) d$spName <- NA_character_
  if (!("cell" %in% names(d))) d$cell <- NA_integer_
  if (!("var" %in% names(d))) d$var <- NA_character_
  if (!("year" %in% names(d))) d$year <- NA_integer_
  if (!("group" %in% names(d))) d$group <- NA_character_
  if (!("orderName" %in% names(d))) d$orderName <- NA_character_
  if (!("familyName" %in% names(d))) d$familyName <- NA_character_
  if (!("redlistCategory" %in% names(d))) d$redlistCategory <- NA_character_
  if (!("propExposed" %in% names(d))) d$propExposed <- NA_real_
  d <- d[, c("spName", "cell", "var", "year", "group", "orderName", "familyName", "redlistCategory", "propExposed")]
  d$spName <- as.factor(as.character(d$spName))
  d$var <- as.factor(as.character(d$var))
  d$cell <- suppressWarnings(as.integer(d$cell))
  d$year <- suppressWarnings(as.integer(d$year))
  d$group <- as.character(d$group)
  d$orderName <- to_title_case(d$orderName)
  d$familyName <- to_title_case(d$familyName)
  d$redlistCategory <- as.character(d$redlistCategory)
  d$propExposed <- suppressWarnings(as.numeric(d$propExposed))
  d
}

load_species_history <- NULL
load_year_df <- NULL
cell_trend_df <- NULL
sp_trend_df <- NULL
# species_year_cells_df (~439 MB) is used only in single-species mode, so it is
# loaded lazily on first use rather than at startup — this lets the default
# hotspot landing page appear ~1.3s sooner. .syc_holder caches the loaded table;
# get_species_year_cells() reads it once on demand. .syc_holder$loaded guards
# against re-reading when the file is missing/empty (result stays NULL).
.syc_holder <- new.env(parent = emptyenv())
.syc_holder$df <- NULL
.syc_holder$loaded <- FALSE
get_species_year_cells <- function() {
  if (!.syc_holder$loaded) {
    .syc_holder$loaded <- TRUE
    if (exists("species_year_cells_path") && file.exists(species_year_cells_path)) {
      log_msg("Lazy-loading species_year_cells (single-species mode) ...")
      .syc_holder$df <- tryCatch(read_qs(species_year_cells_path), error = function(e) NULL)
    }
  }
  .syc_holder$df
}
species_lookup <- empty_app_df()[0, c("spName", "orderName", "familyName"), drop = FALSE]
avail_cells_master <- integer(0)
app_df <- NULL

manifest_obj <- if (file.exists(manifest_path)) {
  tryCatch(read_qs(manifest_path), error = function(e) NULL)
} else NULL

manifest_year_files <- if (!is.null(manifest_obj) && is.list(manifest_obj)) manifest_obj$year_files %||% list() else list()
manifest_year_paths <- vapply(manifest_year_files,
                              function(x) normalizePath(file.path(intDir, x), mustWork = FALSE),
                              character(1))

use_sharded_app_data <- !is.null(manifest_obj) && is.list(manifest_obj) &&
  identical(manifest_obj$version, 1L) &&
  file.exists(species_meta_path) &&
  length(manifest_year_paths) > 0 &&
  any(file.exists(manifest_year_paths))

if (use_sharded_app_data) {
  log_msg("Loading sharded app data manifest ...")
  species_meta <- tryCatch(read_qs(species_meta_path), error = function(e) NULL)
  if (is.null(species_meta) || !nrow(species_meta)) {
    species_meta <- empty_app_df()[0, c("spName", "group", "orderName", "familyName", "redlistCategory"), drop = FALSE]
  }
  species_meta <- species_meta %>%
    mutate(
      spName = as.character(spName),
      group = as.character(group),
      orderName = to_title_case(orderName),
      familyName = to_title_case(familyName),
      redlistCategory = as.character(redlistCategory)
    ) %>%
    distinct(spName, group, orderName, familyName, redlistCategory)

  year_files <- as.list(manifest_year_paths[file.exists(manifest_year_paths)])
  years_avail <- sort(as.integer(names(year_files)))
  year_cache <- new.env(parent = emptyenv())

  load_year_df <- function(year) {
    key <- as.character(year)
    if (exists(key, envir = year_cache, inherits = FALSE)) {
      return(get(key, envir = year_cache, inherits = FALSE))
    }
    fp <- year_files[[key]]
    d <- if (!is.null(fp) && file.exists(fp)) normalize_app_df(read_qs(fp)) else empty_app_df()
    assign(key, d, envir = year_cache)
    d
  }

  load_species_history <- function(sp, sel_vars) {
    sel_vars <- as.character(sel_vars %||% character(0))
    if (!length(sel_vars) || !nzchar(sp %||% "") || !(sp %in% species_avail)) return(empty_app_df())
    bind_rows(lapply(years_avail, function(yr) {
      d <- load_year_df(yr)
      d %>% filter(spName == sp, var %in% sel_vars)
    }))
  }

  if (file.exists(trends_cache_path)) {
    tc <- tryCatch(read_qs(trends_cache_path), error = function(e) NULL)
    if (!is.null(tc) && is.list(tc)) {
      cell_trend_df <- tc$cell_trend_df
      sp_trend_df <- tc$sp_trend_df
      log_msg("Loaded precomputed trend tables: ", basename(trends_cache_path))
    }
    rm(tc)  # release the raw cache list (~124 MB); derived tables are kept
  }
  # species_year_cells_df is now lazy-loaded on first single-species use
  # (see get_species_year_cells); not read here at startup.
  if (is.null(cell_trend_df) || is.null(sp_trend_df)) {
    log_msg("Trend cache missing; rebuilding from sharded app data ...")
    all_year_dfs <- bind_rows(lapply(years_avail, load_year_df))
    cell_trend_df <- all_year_dfs %>%
      distinct(spName, cell, year) %>%
      count(cell, year, name = "n_sp")
    sp_trend_df <- all_year_dfs %>%
      group_by(spName, var, year) %>%
      summarize(mean_prop = mean(propExposed, na.rm = TRUE), .groups = "drop")
  }

  log_msg(sprintf("cell_trend_df: %d rows; sp_trend_df: %d rows",
                  nrow(cell_trend_df), nrow(sp_trend_df)))

  species_avail <- sort(unique(species_meta$spName))
  vars_raw <- as.character(manifest_obj$vars_avail %||% unique(as.character(sp_trend_df$var)))
  vars_raw <- na.omit(vars_raw)
  vars_raw <- vars_raw[nzchar(vars_raw)]
  vars_avail <- vars_raw[order(ifelse(grepl("^temp", vars_raw), 0L, 1L), vars_raw)]
  orders_avail <- sort(na.omit(unique(species_meta$orderName)))
  families_avail <- sort(na.omit(unique(species_meta$familyName)))
  default_year <- if (2025 %in% years_avail) 2025 else max(years_avail)
  species_lookup <- species_meta %>%
    dplyr::select(spName, orderName, familyName) %>%
    distinct()
  avail_cells_master <- as.integer(manifest_obj$avail_cells_all %||% unique(cell_trend_df$cell))
  log_msg(sprintf("Ready. %d years, %d species (sharded app data).",
                  length(years_avail), length(species_avail)))
} else {
  stopifnot(file.exists(allcells_path))
  index_path <- file.path(intDir, "allExpForShiny_year_index_v1.qs")
  log_msg("Loading exposure data: ", basename(allcells_path), " ...")
  all_df <- read_qs(allcells_path)

  app_df <- all_df[, intersect(c("spName", "cell", "var", "year", "group",
                                 "orderName", "familyName", "redlistCategory", "propExposed"), names(all_df))]
  app_df <- normalize_app_df(app_df)
  log_msg(sprintf("app_df: %d rows, %.0f MB (full all_df was %.0f MB)",
                  nrow(app_df),
                  as.numeric(object.size(app_df)) / 1e6,
                  as.numeric(object.size(all_df)) / 1e6))

  year_row_idx <- NULL
  index_mtime <- tryCatch(file.info(allcells_path)$mtime, error = function(e) NA)
  if (file.exists(index_path)) {
    idx_obj <- tryCatch(read_qs(index_path), error = function(e) NULL)
    if (!is.null(idx_obj) && is.list(idx_obj) &&
        identical(idx_obj$version, 1L) &&
        identical(idx_obj$nrows, nrow(app_df)) &&
        identical(as.character(idx_obj$data_mtime), as.character(index_mtime)) &&
        !is.null(idx_obj$year_row_idx)) {
      year_row_idx <- idx_obj$year_row_idx
      log_msg("Loaded cached year index: ", basename(index_path))
    }
  }
  if (is.null(year_row_idx)) {
    log_msg("Building row indices by year ...")
    year_row_idx <- split(seq_len(nrow(app_df)), app_df$year)
    names(year_row_idx) <- as.character(names(year_row_idx))
    tryCatch({
      write_qs(list(
        version = 1L,
        nrows = nrow(app_df),
        data_mtime = index_mtime,
        year_row_idx = year_row_idx
      ), index_path)
      log_msg("Saved cached year index: ", basename(index_path))
    }, error = function(e) {
      log_msg("Year index cache write skipped: ", e$message)
    })
  }

  if (file.exists(trends_cache_path)) {
    tc <- tryCatch(read_qs(trends_cache_path), error = function(e) NULL)
    if (!is.null(tc) && is.list(tc) &&
        identical(tc$version, 1L) &&
        identical(tc$nrows, nrow(app_df)) &&
        identical(as.character(tc$data_mtime), as.character(index_mtime))) {
      cell_trend_df <- tc$cell_trend_df
      sp_trend_df   <- tc$sp_trend_df
      log_msg("Loaded cached trend tables: ", basename(trends_cache_path))
    }
    rm(tc)  # release the raw cache list (~124 MB); derived tables are kept
  }
  # species_year_cells_df lazy-loaded on first single-species use (get_species_year_cells)
  if (is.null(cell_trend_df)) {
    log_msg("Building trend tables ...")
    cell_trend_df <- app_df %>%
      distinct(spName, cell, year) %>%
      count(cell, year, name = "n_sp")
    sp_trend_df <- app_df %>%
      group_by(spName, var, year) %>%
      summarize(mean_prop = mean(propExposed, na.rm = TRUE), .groups = "drop")
    tryCatch({
      write_qs(list(
        version = 1L,
        nrows = nrow(app_df),
        data_mtime = index_mtime,
        cell_trend_df = cell_trend_df,
        sp_trend_df   = sp_trend_df
      ), trends_cache_path)
      log_msg("Saved trend tables cache: ", basename(trends_cache_path))
    }, error = function(e) log_msg("Trend cache write skipped: ", e$message))
  }
  log_msg(sprintf("cell_trend_df: %d rows; sp_trend_df: %d rows",
                  nrow(cell_trend_df), nrow(sp_trend_df)))

  rm(all_df); gc()

  load_year_df <- function(year) {
    idx <- year_row_idx[[as.character(year)]]
    if (is.null(idx) || !length(idx)) return(app_df[0, , drop = FALSE])
    app_df[idx, , drop = FALSE]
  }
  load_species_history <- function(sp, sel_vars) {
    sel_vars <- as.character(sel_vars %||% character(0))
    if (!length(sel_vars) || !nzchar(sp %||% "") || !(sp %in% species_avail)) return(app_df[0, , drop = FALSE])
    app_df %>% filter(spName == sp, var %in% sel_vars)
  }

  years_avail    <- sort(as.integer(names(year_row_idx)))
  species_avail  <- sort(unique(as.character(app_df$spName)))
  vars_raw       <- na.omit(unique(as.character(app_df$var)))
  vars_raw       <- vars_raw[nzchar(vars_raw)]
  vars_avail     <- vars_raw[order(ifelse(grepl("^temp", vars_raw), 0L, 1L), vars_raw)]
  orders_avail   <- sort(na.omit(unique(app_df$orderName)))
  families_avail <- sort(na.omit(unique(app_df$familyName)))
  default_year   <- if (2025 %in% years_avail) 2025 else max(years_avail)
  species_lookup <- app_df %>%
    dplyr::select(spName, orderName, familyName) %>%
    distinct()
  avail_cells_master <- sort(unique(app_df$cell[!is.na(app_df$cell) & app_df$cell >= 1 & app_df$cell <= ncell(tpl)]))
  log_msg(sprintf("Ready. %d rows across %d years, %d species.",
                  nrow(app_df), length(year_row_idx), length(levels(app_df$spName))))
}

# Precomputed hotspot rasters (all-vars, no-filter default) for fast initial render
hotspot_raster_cache <- local({
  if (!file.exists(hotspot_cache_path)) return(list())
  raw <- tryCatch(read_qs(hotspot_cache_path), error = function(e) NULL)
  if (is.null(raw) || !is.list(raw)) return(list())
  tpl_r <- tryCatch(raster::raster(template_path), error = function(e) NULL)
  if (is.null(tpl_r)) return(list())
  out <- lapply(raw, function(entry) {
    # NB: raster::setValues() returns the RasterLayer WITH values; wrapping it in
    # raster::raster() copies only the geometry and drops the values (all-NA),
    # which blanked the default hotspot map. Keep the setValues() result directly.
    tryCatch(raster::setValues(tpl_r, as.numeric(entry$cell_vals)),
             error = function(e) NULL)
  })
  out <- out[!vapply(out, is.null, logical(1))]
  log_msg(sprintf("Loaded precomputed hotspot rasters for years: %s",
                  paste(sort(names(out)), collapse = ", ")))
  out
})

# Full-range cell lookup for Mammals and Birds (for grey background overlay)
log_msg("Loading species range cells (Mammals + Birds) ...")
range_cells_mb <- tryCatch(
  if (file.exists(range_cells_mb_path)) qs2::qs_read(range_cells_mb_path) else list(),
  error = function(e) { log_msg("range_cells_mb load failed: ", e$message); list() }
)
log_msg(sprintf("  range_cells_mb: %d species", length(range_cells_mb)))

# Per-variable colours — keyed to the actual var strings in the data
VAR_COLORS <- c(
  temp__12_up      = "#d73027",   # annual hot
  temp__12_lo      = "#4575b4",   # annual cold
  temp__3__max_up  = "#f46d43",   # seasonal hot extreme
  temp__3__min_lo  = "#313695",   # seasonal cold extreme
  precip__12_up    = "#1a9850",   # annual wet
  precip__12_lo    = "#9400D3",   # annual dry
  precip__3__max_up = "#33a02c",  # seasonal wet extreme
  precip__3__min_lo = "#6a0dad"   # seasonal dry extreme
)
var_col <- function(v) {
  col <- VAR_COLORS[as.character(v)]
  if (!is.na(col)) col else "#e41a1c"
}

VAR_LABELS <- c(
  temp__12_up       = "Temperature: annual upper extreme",
  temp__12_lo       = "Temperature: annual lower extreme",
  temp__3__max_up   = "Temperature: seasonal max upper extreme",
  temp__3__min_lo   = "Temperature: seasonal min lower extreme",
  precip__12_up     = "Precipitation: annual upper extreme",
  precip__12_lo     = "Precipitation: annual lower extreme",
  precip__3__max_up = "Precipitation: seasonal max upper extreme",
  precip__3__min_lo = "Precipitation: seasonal min lower extreme"
)
var_label <- function(v) {
  lab <- VAR_LABELS[as.character(v)]
  if (!is.na(lab)) lab else as.character(v)
}
var_choices <- setNames(as.list(vars_avail), vapply(vars_avail, var_label, character(1)))

hotspot_transform <- function(x, min_x, max_x, power = 1.8) {
  x <- as.numeric(x)
  if (!length(x) || !is.finite(min_x) || !is.finite(max_x) || max_x <= min_x) {
    return(rep(0, length(x)))
  }
  scaled <- (pmin(pmax(x, min_x), max_x) - min_x) / (max_x - min_x)
  scaled ^ power
}

species_display_label <- function(x) gsub("_", " ", as.character(x), fixed = TRUE)
normalize_species_value <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (!nzchar(x)) return("")
  cand <- gsub("[[:space:]]+", "_", x)
  if (cand %in% species_avail) cand else x
}
species_choice_df <- function(x) {
  data.frame(
    value = as.character(x),
    label = species_display_label(x),
    stringsAsFactors = FALSE
  )
}
species_choices <- species_choice_df(species_avail)

# ── helpers ───────────────────────────────────────────────────────────────────
get_species_year_cell_list <- function(sp, year, sel_vars) {
  sp <- normalize_species_value(sp)
  sel_vars <- as.character(sel_vars %||% character(0))
  sel_vars <- sel_vars[nzchar(sel_vars)]
  if (!length(sel_vars) || !nzchar(sp) || !(sp %in% species_avail)) return(list())

  syc <- get_species_year_cells()
  if (!is.null(syc) && nrow(syc)) {
    yr <- suppressWarnings(as.integer(year))
    rows <- syc %>%
      filter(year == yr, spName == sp, var %in% sel_vars)
    if (!nrow(rows)) return(list())
    out <- setNames(rows$cells, as.character(rows$var))
    return(out[intersect(sel_vars, names(out))])
  }

  d <- single_species_year_df()
  out <- list()
  for (v in intersect(sel_vars, as.character(unique(d$var)))) {
    cv <- unique(d$cell[as.character(d$var) == v])
    cv <- cv[!is.na(cv) & cv >= 1 & cv <= ncell(tpl)]
    if (length(cv)) out[[v]] <- cv
  }
  out
}

build_species_layers_from_cell_list <- function(cell_list) {
  out <- list()
  for (v in names(cell_list)) {
    cv <- unique(as.integer(cell_list[[v]]))
    cv <- cv[!is.na(cv) & cv >= 1 & cv <= ncell(tpl)]
    if (!length(cv)) next
    rv <- rep(NA_real_, ncell(tpl)); rv[cv] <- 1
    r_full <- setValues(tpl, rv)
    xy   <- terra::xyFromCell(tpl, cv)
    xbuf <- max((max(xy[, 1]) - min(xy[, 1])) * 0.5, 5)
    ybuf <- max((max(xy[, 2]) - min(xy[, 2])) * 0.5, 5)
    e    <- terra::ext(
      max(-180, min(xy[, 1]) - xbuf),
      min(180,  max(xy[, 1]) + xbuf),
      max(-61,  min(xy[, 2]) - ybuf),
      min(86,   max(xy[, 2]) + ybuf)
    )
    out[[v]] <- raster(terra::crop(r_full, e))
  }
  out
}

build_species_layers <- function(d, sel_vars) {
  out <- list()
  for (v in intersect(sel_vars, as.character(unique(d$var)))) {
    cv <- unique(d$cell[as.character(d$var) == v])
    cv <- cv[!is.na(cv) & cv >= 1 & cv <= ncell(tpl)]
    if (!length(cv)) next
    rv <- rep(NA_real_, ncell(tpl)); rv[cv] <- 1
    r_full <- setValues(tpl, rv)
    # Crop to species extent + buffer so the overlay PNG is tight around the
    # species range — prevents single-cell species from vanishing in a
    # world-scale PNG rendered at low resolution.
    xy   <- terra::xyFromCell(tpl, cv)
    xbuf <- max((max(xy[, 1]) - min(xy[, 1])) * 0.5, 5)
    ybuf <- max((max(xy[, 2]) - min(xy[, 2])) * 0.5, 5)
    e    <- terra::ext(
      max(-180, min(xy[, 1]) - xbuf),
      min(180,  max(xy[, 1]) + xbuf),
      max(-61,  min(xy[, 2]) - ybuf),
      min(86,   max(xy[, 2]) + ybuf)
    )
    out[[v]] <- raster(terra::crop(r_full, e))
  }
  out
}

build_hotspot_raster <- function(d) {
  counts <- d %>%
    filter(!is.na(cell), cell >= 1, cell <= ncell(tpl)) %>%
    group_by(cell) %>%
    summarize(n = n_distinct(spName), .groups = "drop")
  if (!nrow(counts)) return(NULL)
  rv <- rep(NA_real_, ncell(tpl)); rv[counts$cell] <- counts$n
  raster(setValues(tpl, rv))
}

# Per-cell exposed-species count for a single year, taken from the precomputed
# cell_trend_df (all variables, unfiltered). Returns a length-ncell(tpl) numeric
# vector; cells with no exposed species are 0. Used for the two-year change map.
hotspot_count_vec_for_year <- function(yr) {
  v <- rep(0, ncell(tpl))
  if (is.null(cell_trend_df) || !nrow(cell_trend_df)) return(v)
  yr <- suppressWarnings(as.integer(yr))
  d <- cell_trend_df[!is.na(cell_trend_df$year) & cell_trend_df$year == yr,
                     c("cell", "n_sp"), drop = FALSE]
  d <- d[!is.na(d$cell) & d$cell >= 1 & d$cell <= ncell(tpl), , drop = FALSE]
  if (nrow(d)) v[d$cell] <- as.numeric(d$n_sp)
  v
}

# Fast default-view hotspot raster straight from the in-memory cell_trend_df
# (all variables, no filters) — the same quantity build_hotspot_raster()
# recomputes from a disk shard, but ~35x faster because cell_trend_df is already
# the precomputed per-cell/per-year species count. Results are memoized per year
# so animation loops are instant after the first pass. Only valid for the DEFAULT
# view; filtered/threatened/DD views must use the exact per-row path.
.default_raster_memo <- new.env(parent = emptyenv())
default_hotspot_raster <- function(yr) {
  key <- as.character(yr)
  if (!is.null(.default_raster_memo[[key]])) return(.default_raster_memo[[key]])
  v <- hotspot_count_vec_for_year(yr)          # length ncell(tpl), 0 where none
  if (!any(v > 0)) { .default_raster_memo[[key]] <- NULL; return(NULL) }
  v[v == 0] <- NA_real_
  rr <- raster(setValues(tpl, v))
  .default_raster_memo[[key]] <- rr
  rr
}

# Diverging change raster: species-count(target) - species-count(baseline).
# Cells that are 0 in BOTH years stay NA (transparent) so only cells that
# gained or lost exposure are drawn.
build_change_raster <- function(baseline_year, target_year) {
  vb <- hotspot_count_vec_for_year(baseline_year)
  vt <- hotspot_count_vec_for_year(target_year)
  delta <- vt - vb
  delta[vb == 0 & vt == 0] <- NA_real_
  if (!any(is.finite(delta))) return(NULL)
  raster(setValues(tpl, delta))
}

safe_read_json <- function(url) {
  tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
}

BOUNDARY_CACHE_DIR <- file.path(intDir, "boundary_cache")
if (!dir.exists(BOUNDARY_CACHE_DIR)) dir.create(BOUNDARY_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

geoboundaries_country_index <- local({
  # Cache the parsed country index to disk so launches are instant and work
  # offline. The geoBoundaries ADM0 list is effectively static, so a cached copy
  # is safe to reuse indefinitely; delete the file to force a refresh.
  cache_file <- file.path(BOUNDARY_CACHE_DIR, "geoboundaries_country_index_v1.qs")
  if (file.exists(cache_file)) {
    cached <- tryCatch(read_qs(cache_file), error = function(e) NULL)
    if (is.data.frame(cached) && nrow(cached) > 0) {
      log_msg("Loaded cached geoBoundaries country index: ", basename(cache_file))
      return(cached)
    }
  }
  x <- safe_read_json("https://www.geoboundaries.org/api/current/gbOpen/ALL/ADM0/")
  if (is.null(x) || !is.data.frame(x)) return(data.frame())
  out <- x %>%
    dplyr::transmute(
      iso3 = as.character(boundaryISO),
      country = as.character(boundaryName)
    ) %>%
    filter(!is.na(iso3), nchar(iso3) == 3, !is.na(country), nzchar(country)) %>%
    distinct()
  out <- out[order(out$country, out$iso3), , drop = FALSE]
  if (nrow(out) > 0) {
    tryCatch({
      write_qs(out, cache_file)
      log_msg("Saved geoBoundaries country index cache: ", basename(cache_file))
    }, error = function(e) log_msg("Country index cache write skipped: ", e$message))
  }
  out
})

gadm_gpkg_path <- function(iso3) {
  iso3 <- toupper(trimws(as.character(iso3 %||% "")))
  if (!nzchar(iso3)) return(NA_character_)
  file.path(BOUNDARY_CACHE_DIR, sprintf("gadm41_%s.gpkg", iso3))
}

ensure_gadm_file <- function(iso3) {
  iso3 <- toupper(trimws(as.character(iso3 %||% "")))
  if (!nzchar(iso3)) stop("Missing ISO3 code")
  out <- gadm_gpkg_path(iso3)
  if (file.exists(out)) return(out)
  url <- sprintf("https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_%s.gpkg", iso3)
  tmp <- tempfile(fileext = ".gpkg")
  utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
  file.copy(tmp, out, overwrite = TRUE)
  out
}

read_gadm_level <- function(iso3, level) {
  stopifnot(level %in% c(0L, 1L, 2L))
  gpkg <- ensure_gadm_file(iso3)
  lyr <- sprintf("ADM_ADM_%d", as.integer(level))
  terra::vect(gpkg, layer = lyr)
}

spatvector_to_leaflet_geojson <- function(v) {
  tmp <- tempfile(fileext = ".geojson")
  terra::writeVector(v, tmp, filetype = "GeoJSON", overwrite = TRUE)
  paste(readLines(tmp, warn = FALSE), collapse = "\n")
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(
    tags$style(HTML(
    "body{font-size:13px}
     #click_panel{margin-top:12px;padding:10px;border:1px solid #ddd;border-radius:4px;background:#fafafa}
     #click_panel h4{margin-top:4px}
     #polygon_panel{margin-top:12px;padding:10px;border:1px solid #e67e00;border-radius:4px;background:#fffaf5}
     #polygon_panel h4{margin-top:4px}
    .sidebar-panel{font-size:1.1em}
    .sidebar .well{padding:10px 12px}
    .sidebar hr{margin:8px 0}
    .sidebar h5{margin:6px 0 4px 0}
    .sidebar .form-group{margin-bottom:6px}
    .sidebar .control-label{margin-bottom:2px}
    .sidebar .radio{margin-top:2px;margin-bottom:2px}
    .sidebar .radio-inline,.sidebar .checkbox-inline{padding-top:2px;padding-bottom:2px}
    .sidebar .checkbox{margin-top:3px;margin-bottom:3px}
    .sidebar .selectize-control{margin-bottom:4px}
    .sidebar .selectize-input{min-height:30px;padding-top:4px;padding-bottom:4px}
    .sidebar .btn{padding-top:4px;padding-bottom:4px}
    .sidebar p{margin:3px 0}
  .btn,.radio label,.checkbox label{min-height:34px}
  .radio-inline,.checkbox-inline{padding:6px 10px}
  .btn-sm{min-height:34px}
  .btn:focus-visible,.selectize-input:focus-within,input:focus-visible,select:focus-visible{outline:3px solid #1e90ff;outline-offset:1px}
     .tip{display:inline-block;margin-left:5px;width:15px;height:15px;line-height:15px;text-align:center;
          font-size:10px;font-weight:bold;color:#fff;background:#888;border-radius:50%;
          cursor:help;vertical-align:middle;position:relative}
  .tip:hover::after,.tip:focus::after,.tip.tip-open::after{content:attr(data-tip);position:absolute;left:20px;top:-4px;
          background:#333;color:#fff;padding:5px 8px;border-radius:4px;font-size:11px;
          font-weight:normal;white-space:normal;width:200px;z-index:9999;line-height:1.4}
  @media (pointer: coarse){
    .tip{width:18px;height:18px;line-height:18px;font-size:11px}
    .tip:hover::after,.tip:focus::after,.tip.tip-open::after{position:fixed;left:12px;right:12px;top:auto;bottom:12px;width:auto;max-width:none;font-size:12px}
  }
     #app-progress-wrap{display:none;margin:6px 0 2px 0}
     #app-progress-bar{height:8px;border-radius:4px;background:#3498db;
          transition:width 0.25s ease;width:0%}
    #app-progress-label{font-size:11px;color:#555;margin-top:2px}
    .leaflet-control .legend{background:#fff !important;opacity:1 !important}
    .map-empty-note{background:rgba(255,255,255,0.92);border:1px solid #d9534f;border-radius:6px;
         padding:10px 16px;font-size:14px;font-weight:600;color:#a94442;
         box-shadow:0 1px 4px rgba(0,0,0,0.25);text-align:center;max-width:320px}
    .summary-card{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin:6px 0 2px 0;
         padding:8px 12px;border:1px solid #d9e2ec;border-radius:6px;background:#f5f7fa;font-size:12px}
    .summary-card .stat{display:flex;flex-direction:column;line-height:1.15}
    .summary-card .stat .num{font-size:18px;font-weight:700;color:#2c3e50}
    .summary-card .stat .lab{font-size:11px;color:#607086;text-transform:uppercase;letter-spacing:0.03em}
    .summary-card .stat.wide{min-width:150px}
    .summary-card .stat.wide .num{font-size:13px;font-weight:600}
    "
    )),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('setProgress', function(msg) {
        var wrap  = document.getElementById('app-progress-wrap');
        var bar   = document.getElementById('app-progress-bar');
        var label = document.getElementById('app-progress-label');
        if (!wrap || !bar || !label) return;
        if (msg.pct < 0) {
          wrap.style.display = 'none';
          bar.style.width = '0%';
          label.textContent = '';
        } else {
          wrap.style.display = 'block';
          bar.style.width = Math.min(100, msg.pct) + '%';
          label.textContent = msg.label || '';
        }
      });
      document.addEventListener('DOMContentLoaded', function(){
        var tips = document.querySelectorAll('.tip');
        tips.forEach(function(t){
          t.setAttribute('tabindex', '0');
          t.setAttribute('role', 'button');
          t.setAttribute('aria-label', 'Help');
        });
        document.addEventListener('click', function(ev){
          tips.forEach(function(t){ if (t !== ev.target) t.classList.remove('tip-open'); });
          if (ev.target.classList && ev.target.classList.contains('tip')) {
            ev.target.classList.toggle('tip-open');
            ev.preventDefault();
          }
        });
      });
    "))
  ),

  tags$div(style = "display:flex;align-items:center;justify-content:space-between;padding:10px 15px 4px 15px",
    tags$h2("Species exposure to extreme climate", style = "margin:0"),
    tags$a(href = "https://speciesexposure.github.io/", target = "_blank",
           style = "font-size:13px;color:#3498db;text-decoration:none;white-space:nowrap",
           "speciesexposure.github.io ↗")
  ),
  sidebarLayout(
    sidebarPanel(width = 3, class = "sidebar-panel",
      radioButtons("mode", tags$span("Map mode", tags$span("?", class="tip", `data-tip`="Single species: map one species' exposed cells. Hotspot: map how many species are exposed per cell.")),
                   choices  = list("Single species" = "single",
                                   "Hotspot: # species/cell" = "hotspot"),
                   selected = "hotspot", inline = TRUE),
      actionButton("reset_all", "Reset all controls", width = "100%"),
      hr(),
      fluidRow(
        column(4, tags$label("Year", tags$span("?", class="tip", `data-tip`="Select the year to display exposure for."), style = "padding-top:7px;font-weight:600")),
        column(8, selectInput("year", NULL, choices = years_avail, selected = default_year))
      ),
      conditionalPanel("input.mode == 'hotspot'",
        fluidRow(
          column(6, actionButton("play_years", tags$span(icon("play"), " Play years"),
                                  class = "btn-sm", width = "100%")),
          column(6, actionButton("pause_years", tags$span(icon("pause"), " Pause"),
                                  class = "btn-sm", width = "100%"))
        ),
        div(style = "font-size:11px;color:#607086;margin:2px 0 2px 2px",
            tags$span("?", class="tip", `data-tip`="Animate the map forward through years automatically. It loops back to the first year at the end. Pause to stop."),
            " Auto-advance through years"),
        sliderInput("play_speed",
                    tags$span("Playback speed (sec/year)",
                              tags$span("?", class="tip", `data-tip`="Seconds each year is shown during playback. Lower is faster. Takes effect on the next frame.")),
                    min = 0.25, max = 2.5, value = 0.75, step = 0.25, ticks = FALSE),
        hr(),
        checkboxInput("change_mode",
                      tags$span("Show change vs baseline year",
                                tags$span("?", class="tip", `data-tip`="Map the change in number of exposed species per cell between a baseline year and the selected year. Blue = fewer species exposed than baseline; red = more.")),
                      value = FALSE),
        conditionalPanel("input.change_mode == true",
          fluidRow(
            column(4, tags$label("Baseline", style = "padding-top:7px;font-weight:600")),
            column(8, selectInput("baseline_year", NULL, choices = years_avail,
                                  selected = min(years_avail)))
          )
        )
      ),
      checkboxGroupInput("sel_vars", tags$span("Climate variables", tags$span("?", class="tip", `data-tip`="Choose which extreme climate variables to include. The map shows cells where at least one selected variable exceeds the species' historical threshold.")), choices = var_choices, selected = vars_avail),
      fluidRow(
        column(6, actionButton("vars_all", "Select all", width = "100%")),
        column(6, actionButton("vars_none", "Deselect all", width = "100%"))
      ),
      hr(),
      conditionalPanel("input.mode == 'single'",
        h5(tags$span("Species", tags$span("?", class="tip", `data-tip`="Type to search for a species. The map will show all grid cells where that species is exposed to the selected climate variable(s) in the chosen year."))),
        checkboxInput("flt_exposed_only", "Only show species with exposure in selected year/variable(s)", value = FALSE),
        selectizeInput("species", NULL,
                       choices  = NULL,
                       selected = character(0),
                       options  = list(maxItems = 1, placeholder = "Type species name...",
                                       valueField = "value",
                                       labelField = "label",
                                       searchField = c("label", "value")))
      ),
      conditionalPanel("input.mode == 'hotspot'",
        h5(tags$span("Filter species (hotspot mode)", tags$span("?", class="tip", `data-tip`="Restrict which species count toward the hotspot map. Leave all blank to include all species."))),
        fluidRow(
          column(6, checkboxInput("flt_threatened", "Only threatened (CR/EN/VU)", value = FALSE)),
          column(6, checkboxInput("flt_data_deficient", "Only Data Deficient (DD)", value = FALSE))
        ),
        div(style = "margin-bottom:4px;",
            tags$label("Taxonomic groups",
                       tags$span("?", class = "tip",
                         `data-tip` = "Restrict the hotspot map to one or more taxonomic groups. All four are included by default. Note: reptile and amphibian ranges are more coarsely represented, so their hotspots behave differently from birds and mammals."),
                       style = "font-weight:600"),
            checkboxGroupInput("flt_groups", NULL,
                               choices = c("Amphibians", "Birds", "Mammals", "Reptiles"),
                               selected = c("Amphibians", "Birds", "Mammals", "Reptiles"),
                               inline = TRUE)
        ),
        if (length(orders_avail))
          fluidRow(
            column(4, tags$label("Order", style = "padding-top:7px;font-weight:600")),
            column(8, selectizeInput("flt_order", NULL,
                           choices = c("", orders_avail), selected = "", multiple = TRUE,
                           options = list(placeholder = "All")))
          )
        else helpText("Order filter unavailable — no metadata loaded."),
        if (length(families_avail))
          fluidRow(
            column(4, tags$label("Family", style = "padding-top:7px;font-weight:600")),
            column(8, selectizeInput("flt_family", NULL,
                           choices = c("", families_avail), selected = "", multiple = TRUE,
                           options = list(placeholder = "All")))
          )
        else helpText("Family filter unavailable — no metadata loaded."),
        hr(),
        h5(tags$span("Political unit", tags$span("?", class="tip", `data-tip`="Select a country, and optionally a state or county. Click 'Load boundary' to highlight that region and list all exposed species within it."))),
        fluidRow(
          column(4, tags$label("Country", style = "padding-top:7px;font-weight:600")),
          column(8, selectizeInput("pol_country", NULL, choices = NULL,
                                   selected = "", multiple = FALSE,
                                   options = list(placeholder = "Select...")))
        ),
        fluidRow(
          column(4, tags$label("State", style = "padding-top:7px;font-weight:600")),
          column(8, selectizeInput("pol_state", NULL, choices = NULL,
                                   selected = "", multiple = FALSE,
                                   options = list(placeholder = "All")))
        ),
        fluidRow(
          column(4, tags$label("County", style = "padding-top:7px;font-weight:600")),
          column(8, selectizeInput("pol_county", NULL, choices = NULL,
                                   selected = "", multiple = FALSE,
                                   options = list(placeholder = "All")))
        ),
        fluidRow(
          column(6, actionButton("load_pol_unit", "Load boundary", class = "btn-sm", width = "100%")),
          column(6, actionButton("clear_polygon", "Clear",          class = "btn-sm", width = "100%"))
        ),
        hr(),
        h5(tags$span("Upload shapefile", tags$span("?", class="tip", `data-tip`="Select all components of a shapefile at once (.shp, .dbf, .shx, .prj). The boundary will be drawn on the map and exposed species listed below."))),
        tags$p(style = "font-size:11px;color:#666;margin-bottom:4px",
               "Select all shapefile components (.shp, .dbf, .shx, .prj) at once."),
        fileInput("shp_upload", NULL, multiple = TRUE,
                  accept = c(".shp", ".dbf", ".shx", ".prj", ".cpg"),
                  buttonLabel = "Browse…", placeholder = "No file selected"),
        hr(),
        h5(tags$span("Exposure threshold",
                     tags$span("?", class = "tip",
                       `data-tip` = "Count a species toward the hotspot map in a given year only if at least this fraction of its range is exposed that year. Applied per species per year. Default 10%. Raise it to down-weight wide-ranging species whose exposure is a small fraction of a large range. Note: species with <1% of their range exposed in a year are not included in the dataset, to keep it light — so 1% is the lowest available threshold."))),
        sliderInput("excl_threshold", "Minimum % of range exposed (per species-year)",
                    min = 1, max = 100, value = 10, step = 1, post = "%", ticks = FALSE)
      )
    ),
    mainPanel(width = 9,
      uiOutput("summary_card"),
      uiOutput("map_year_title"),
      leafletOutput("map", height = "778px"),
      div(style = "margin-top:6px;",
          downloadButton("download_map_below", "Download map raster (GeoTIFF)", class = "btn-sm"),
          downloadButton("download_summary", "Download summary stats (CSV)", class = "btn-sm")),
      conditionalPanel("input.mode == 'hotspot'",
        div(style = "margin-top:8px;",
            plotOutput("trend_ts", height = "150px"),
            div(style = "font-size:11px;color:#666;text-align:center;margin-top:-2px;",
                "Total exposed species per year at the current exposure threshold (all taxa; order/family/threatened/DD filters not applied to this trend). Vertical line = displayed year.")
        ),
        div(style = "display:flex;justify-content:flex-end;align-items:center;gap:18px;margin-top:4px;margin-bottom:2px;",
            div(style = "display:flex;align-items:center;",
                checkboxInput("fix_color_scale", "Fix colour scale", value = FALSE),
                tags$span("?", class = "tip",
                  `data-tip` = "When checked, the map's colour scale is locked to a fixed range instead of rescaling to each year, so colours are comparable across years (useful when stepping through or animating). The fixed range is the span of species counts in 2024 computed at the CURRENT exposure threshold (and taxon filters) — so if you change the threshold, the fixed scale is recomputed to match, since raising the threshold lowers the counts. Leave unchecked to rescale each year to its own min/max.")
            ),
            div(style = "width:260px;",
                sliderInput("hotspot_opacity", "Hotspot opacity", min = 0, max = 1, value = 0.95, step = 0.05)
            )
        )
      ),
      div(id = "app-progress-wrap",
          div(style = "background:#e0e0e0;border-radius:4px;overflow:hidden",
              div(id = "app-progress-bar")),
          div(id = "app-progress-label")
      ),
      uiOutput("status_msg"),
      uiOutput("exposure_msg"),
      conditionalPanel("input.mode == 'single'",
        plotOutput("species_range_trend", height = "200px"),
        h4("Selected species summary"),
        DTOutput("species_summary"),
        br(),
        downloadButton("download_filtered",  "Download species table"),
        downloadButton("download_template",  "Download raster template"),
        downloadButton("download_rscript",   "Download R script"),
        downloadButton("download_map_raster", "Download map raster")
      ),
      conditionalPanel("input.mode == 'hotspot'",
        conditionalPanel("output.polygon_active != 'yes'",
          div(id = "click_panel",
              h4(textOutput("click_title", inline = TRUE)),
              DTOutput("click_table"),
              br(),
              downloadButton("download_click", "Download clicked-cell table")
          ),
          h4("Hotspot diagnostics for clicked cell"),
          DTOutput("hotspot_diag"),
          h4("Clicked cell trend over time"),
          plotOutput("cell_trend", height = "220px"),
          hr()
        ),
        conditionalPanel("output.polygon_active == 'yes'",
          div(id = "polygon_panel",
              h4("Selected area summary"),
              uiOutput("polygon_summary"),
              plotOutput("polygon_trend", height = "180px"),
              div(style = "font-size:11px;color:#666;text-align:center;margin-top:-4px;",
                  "Total exposed species\u00d7cells per year within the selected area (a species is counted once per cell it occupies; all taxa, no threshold applied to this trend). Vertical line = displayed year."),
              h4("Species in selected area"),
              DTOutput("polygon_table"),
              br(),
              downloadButton("download_polygon",       "Download polygon table"),
              downloadButton("download_map_raster_hs", "Download map raster")
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  exposure_msg_rv <- reactiveVal(NULL)
  status_msg_rv <- reactiveVal("Ready")
  polygon_empty_reason_rv <- reactiveVal(NULL)

  set_status <- function(msg, done = FALSE) {
    stamp <- format(Sys.time(), "%H:%M:%S")
    status_msg_rv(paste0("[", stamp, "] ", msg))
    if (isTRUE(done)) showNotification(msg, type = "message", duration = 2)
  }

  show_progress <- function(pct, label = "") {
    session$sendCustomMessage("setProgress", list(pct = pct, label = label))
  }
  hide_progress <- function() {
    session$sendCustomMessage("setProgress", list(pct = -1, label = ""))
  }

  # Centered "no data" annotation drawn directly on the map, so an empty result
  # reads as an intentional state rather than a broken/blank map.
  show_map_empty <- function(msg = "No exposed cells for the current selection.") {
    leafletProxy("map") %>%
      removeControl("empty_note") %>%
      addControl(html = sprintf("<div class='map-empty-note'>%s</div>", msg),
                 position = "topright", layerId = "empty_note")
  }
  clear_map_empty <- function() {
    leafletProxy("map") %>% removeControl("empty_note")
  }

  ensure_prop_exposed <- function(d) {
    if (!"propExposed" %in% names(d)) d$propExposed <- NA_real_
    d$propExposed <- suppressWarnings(as.numeric(d$propExposed))
    d
  }

  # Normalize IUCN labels: the data stores them with underscores
  # (e.g. "Critically_Endangered", "Data_Deficient"), so map "_" -> " " before
  # matching. Without this, underscored categories like Critically_Endangered
  # slip past a space-based match.
  .iucn_norm <- function(x) tolower(trimws(gsub("_", " ", as.character(x))))
  is_threatened <- function(x) {
    .iucn_norm(x) %in% c("critically endangered", "endangered", "vulnerable", "cr", "en", "vu")
  }
  is_data_deficient <- function(x) {
    .iucn_norm(x) %in% c("data deficient", "dd")
  }

  ALL_GROUPS <- c("Amphibians", "Birds", "Mammals", "Reptiles")
  apply_hotspot_filters <- function(d, ord, fam, iuc, threatened = FALSE, data_deficient = FALSE,
                                    groups = ALL_GROUPS) {
    ord <- ord[nzchar(ord)]
    fam <- fam[nzchar(fam)]
    iuc <- iuc[nzchar(iuc)]
    groups <- groups[nzchar(groups)]
    # Only apply a group filter when it is a real subset (empty or all-four = no-op)
    if (length(groups) && !setequal(groups, ALL_GROUPS) && "group" %in% names(d))
      d <- d %>% filter(group %in% groups)
    if (length(ord)) d <- d %>% filter(orderName %in% ord)
    if (length(fam)) d <- d %>% filter(familyName %in% fam)
    if (length(iuc)) d <- d %>% filter(redlistCategory %in% iuc)
    if (isTRUE(threatened)) d <- d %>% filter(is_threatened(redlistCategory))
    if (isTRUE(data_deficient)) d <- d %>% filter(is_data_deficient(redlistCategory))
    d
  }

  is_valid_species_selection <- function(sp) {
    sp_norm <- normalize_species_value(sp)
    nzchar(sp_norm) && (sp_norm %in% species_avail)
  }

  updateSelectizeInput(session, "species",
                       choices = species_choices,
                       selected = species_avail[1],
                       server = TRUE)

  country_choices <- if (nrow(geoboundaries_country_index)) {
    setNames(as.list(geoboundaries_country_index$iso3),
             paste0(geoboundaries_country_index$country, " (", geoboundaries_country_index$iso3, ")"))
  } else {
    list("United States (USA)" = "USA")
  }
  updateSelectizeInput(session, "pol_country", choices = country_choices, selected = "USA", server = TRUE)
  updateSelectizeInput(session, "pol_state", choices = character(0), selected = "", server = TRUE)
  updateSelectizeInput(session, "pol_county", choices = character(0), selected = "", server = TRUE)

  # ── Shareable-link state: read/write key inputs via the URL query string ──
  # On load, parse ?...= parameters and apply them; thereafter, mirror the
  # current view back into the URL so it can be copied and shared. A short
  # restore window prevents the initial write from clobbering restored values.
  url_restoring <- reactiveVal(TRUE)
  # Apply saved URL state exactly once, on the first flush, then release the lock
  # so subsequent input changes are mirrored back into the URL. Tying this to the
  # session flush (rather than a wall-clock timer) keeps it within the session
  # lifecycle and avoids stale callbacks firing after the session closes.
  session$onFlushed(function() {
    qs <- tryCatch(parseQueryString(session$clientData$url_search), error = function(e) list())
    isolate({
      if (!is.null(qs$mode) && qs$mode %in% c("single", "hotspot"))
        updateRadioButtons(session, "mode", selected = qs$mode)
      if (!is.null(qs$year) && qs$year %in% as.character(years_avail))
        updateSelectInput(session, "year", selected = qs$year)
      if (!is.null(qs$vars)) {
        v <- strsplit(qs$vars, ",", fixed = TRUE)[[1]]
        v <- intersect(v, vars_avail)
        if (length(v)) updateCheckboxGroupInput(session, "sel_vars", selected = v)
      }
      if (!is.null(qs$species)) {
        sp <- normalize_species_value(qs$species)
        if (is_valid_species_selection(sp))
          updateSelectizeInput(session, "species", selected = sp, server = TRUE)
      }
      if (!is.null(qs$threat))
        updateCheckboxInput(session, "flt_threatened", value = identical(qs$threat, "1"))
      if (!is.null(qs$dd))
        updateCheckboxInput(session, "flt_data_deficient", value = identical(qs$dd, "1"))
      if (!is.null(qs$change))
        updateCheckboxInput(session, "change_mode", value = identical(qs$change, "1"))
      if (!is.null(qs$baseline) && qs$baseline %in% as.character(years_avail))
        updateSelectInput(session, "baseline_year", selected = qs$baseline)
      if (!is.null(qs$country) && nzchar(qs$country))
        updateSelectizeInput(session, "pol_country", selected = toupper(qs$country), server = TRUE)
      if (!is.null(qs$excl)) {
        ev <- suppressWarnings(as.numeric(qs$excl))
        if (is.finite(ev) && ev >= 1 && ev <= 100)
          updateSliderInput(session, "excl_threshold", value = round(ev))
      }
      if (!is.null(qs$groups)) {
        g <- intersect(strsplit(qs$groups, ",")[[1]], ALL_GROUPS)
        if (length(g)) updateCheckboxGroupInput(session, "flt_groups", selected = g)
      }
    })
    # Release the restore lock so later input changes mirror back into the URL
    url_restoring(FALSE)
  }, once = TRUE)

  observe({
    # Mirror current state into the URL (skip during the restore window)
    mode <- input$mode; yr <- input$year
    if (is.null(mode) || is.null(yr)) return()
    if (isTRUE(url_restoring())) return()
    q <- list(
      mode     = mode,
      year     = as.character(yr),
      vars     = paste(input$sel_vars %||% character(0), collapse = ","),
      threat   = if (isTRUE(input$flt_threatened)) "1" else "0",
      dd       = if (isTRUE(input$flt_data_deficient)) "1" else "0",
      change   = if (isTRUE(input$change_mode)) "1" else "0",
      baseline = as.character(input$baseline_year %||% min(years_avail)),
      excl     = as.character(input$excl_threshold %||% 10),
      groups   = paste(input$flt_groups %||% ALL_GROUPS, collapse = ",")
    )
    if (identical(mode, "single")) {
      sp <- normalize_species_value(input$species)
      if (nzchar(sp)) q$species <- sp
    }
    ct <- toupper(trimws(input$pol_country %||% ""))
    if (nzchar(ct)) q$country <- ct
    qstr <- paste0("?", paste(names(q), vapply(q, function(v) utils::URLencode(as.character(v), reserved = TRUE), character(1)),
                              sep = "=", collapse = "&"))
    updateQueryString(qstr, mode = "replace", session = session)
  })

  gadm_cache <- new.env(parent = emptyenv())
  shp_cache <- new.env(parent = emptyenv())
  intersection_cache <- new.env(parent = emptyenv())

  get_gadm_cached <- function(iso3, level) {
    key <- paste0(toupper(iso3), "_", as.integer(level))
    if (exists(key, envir = gadm_cache, inherits = FALSE)) return(get(key, envir = gadm_cache, inherits = FALSE))
    v <- read_gadm_level(iso3, level)
    assign(key, v, envir = gadm_cache)
    v
  }

  avail_cells_all <- sort(unique(avail_cells_master[!is.na(avail_cells_master) & avail_cells_master >= 1 & avail_cells_master <= ncell(tpl)]))
  avail_pts <- if (length(avail_cells_all)) {
    xy <- terra::xyFromCell(tpl, avail_cells_all)
    terra::vect(xy, type = "points", crs = terra::crs(tpl))
  } else NULL

  intersect_cached <- function(unit_key, unit_v) {
    if (exists(unit_key, envir = intersection_cache, inherits = FALSE)) {
      return(get(unit_key, envir = intersection_cache, inherits = FALSE))
    }
    if (is.null(avail_pts) || !length(avail_cells_all)) {
      out <- integer(0)
    } else {
      inside <- tryCatch(
        terra::is.related(avail_pts, unit_v, "intersects"),
        error = function(e) rep(FALSE, length(avail_cells_all))
      )
      out <- avail_cells_all[inside]
    }
    assign(unit_key, out, envir = intersection_cache)
    out
  }

  # Re-filter species list when checkbox or year/vars change
  observe({
    req(input$mode == "single")
    if (isTRUE(input$flt_exposed_only)) {
      d <- year_df() %>% filter(var %in% (input$sel_vars %||% character(0)))
      spp <- sort(unique(as.character(d$spName)))
      if (!length(spp)) spp <- species_avail
    } else {
      spp <- species_avail
    }
    cur <- normalize_species_value(isolate(input$species))
    updateSelectizeInput(session, "species",
                         choices = species_choice_df(spp),
                         selected = if (nzchar(cur %||% "") && cur %in% spp) cur else spp[1],
                         server = TRUE)
  })

  observeEvent(input$vars_all, {
    updateCheckboxGroupInput(session, "sel_vars", selected = vars_avail)
  })

  observeEvent(input$vars_none, {
    updateCheckboxGroupInput(session, "sel_vars", selected = character(0))
  })

  output$exposure_msg <- renderUI({
    msg <- exposure_msg_rv()
    if (!is.null(msg))
      div(style = "color:#c0392b; padding:6px 10px; margin-top:6px;",
          icon("exclamation-circle"), " ", msg)
  })

  output$status_msg <- renderUI({
    msg <- status_msg_rv()
    div(style = "color:#2c3e50; background:#f5f7fa; border:1px solid #d9e2ec; border-radius:4px; padding:6px 10px; margin-top:6px;",
        icon("spinner"), " ", msg)
  })

  # O(1) year lookup: picks the pre-split slice instead of scanning all rows
  year_df <- reactive({
    load_year_df(input$year)
  })

  # filtered_df: apply var + optional hotspot filters on the already-small year slice.
  # bindCache means identical inputs reuse the previous result instantly.
  filtered_df <- reactive({
    sel <- input$sel_vars
    d   <- year_df() %>% filter(var %in% sel)
    if (input$mode == "hotspot") {
      d <- apply_hotspot_filters(d,
                                 input$flt_order %||% character(0),
                                 input$flt_family %||% character(0),
                                 character(0),
                                 input$flt_threatened %||% FALSE,
                                 input$flt_data_deficient %||% FALSE,
                                 input$flt_groups %||% ALL_GROUPS)
      # Exposure threshold (per species-year): keep only species whose exposed
      # fraction of range that year (propExposed) is >= the slider value. The
      # shards already floor at 1%, so thr = 1% is a no-op that reproduces the
      # delivered map; higher thresholds down-weight wide-ranging species.
      thr <- (as.numeric(input$excl_threshold %||% 10)) / 100
      if (is.finite(thr) && thr > 0.01 && "propExposed" %in% names(d)) {
        d <- d %>% filter(is.na(.data$propExposed) | .data$propExposed >= thr)
      }
    }
    d
  }) |> bindCache(input$year, input$sel_vars, input$mode,
                  input$flt_order, input$flt_family, input$flt_threatened,
                  input$flt_data_deficient, input$excl_threshold, input$flt_groups)

  # Fixed colour scale reference: the per-cell species-count range for a fixed
  # reference year (2024) under the CURRENT variable/taxon/threshold settings.
  # Used when input$fix_color_scale is checked so colours are comparable across
  # years. Memoized per settings key so animation does not rebuild it each frame.
  FIXED_SCALE_YEAR <- 2024L
  .fixed_scale_memo <- new.env(parent = emptyenv())
  fixed_scale_range <- reactive({
    sel <- input$sel_vars %||% character(0)
    if (!length(sel)) return(NULL)
    thr <- (as.numeric(input$excl_threshold %||% 10)) / 100
    key <- paste(FIXED_SCALE_YEAR, paste(sort(sel), collapse = ","),
                 paste(sort(input$flt_order  %||% ""), collapse = ","),
                 paste(sort(input$flt_family %||% ""), collapse = ","),
                 isTRUE(input$flt_threatened), isTRUE(input$flt_data_deficient),
                 paste(sort(input$flt_groups %||% ALL_GROUPS), collapse = ","),
                 thr, sep = "|")
    if (!is.null(.fixed_scale_memo[[key]])) return(.fixed_scale_memo[[key]])
    d <- tryCatch(load_year_df(FIXED_SCALE_YEAR), error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    d <- d %>% filter(var %in% sel)
    d <- apply_hotspot_filters(d,
                               input$flt_order %||% character(0),
                               input$flt_family %||% character(0),
                               character(0),
                               input$flt_threatened %||% FALSE,
                               input$flt_data_deficient %||% FALSE,
                               input$flt_groups %||% ALL_GROUPS)
    if (is.finite(thr) && thr > 0.01 && "propExposed" %in% names(d))
      d <- d %>% filter(is.na(.data$propExposed) | .data$propExposed >= thr)
    rr <- build_hotspot_raster(d)
    if (is.null(rr)) return(NULL)
    v <- raster::values(rr); v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    rng <- c(min(v, na.rm = TRUE), max(v, na.rm = TRUE))
    .fixed_scale_memo[[key]] <- rng
    rng
  })

  single_species_year_df <- reactive({
    req(input$mode == "single", input$species)
    sel <- input$sel_vars %||% character(0)
    req(length(sel))
    sp_sel <- normalize_species_value(input$species)
    if (!is_valid_species_selection(sp_sel)) return(empty_app_df())
    d <- year_df()
    keep <- as.character(d$spName) == sp_sel & as.character(d$var) %in% sel
    d[keep, , drop = FALSE]
  }) |> bindCache(input$year, input$species, input$sel_vars)

  species_summary_df <- reactive({
    req(input$mode == "single", input$species)
    sel <- input$sel_vars
    req(length(sel))
    d <- load_species_history(normalize_species_value(input$species), sel)
    d <- ensure_prop_exposed(d)
    d %>%
      group_by(year) %>%
      summarize(
        Exposed_cells     = n_distinct(cell),
        Mean_prop_exposed = round(mean(.data$propExposed, na.rm = TRUE), 4),
        .groups = "drop"
      ) %>%
      arrange(desc(year))
  }) |> bindCache(input$species, input$sel_vars)

  output$species_summary <- renderDT({
    req(input$mode == "single")
    set_status("Loading species summary...")
    out <- datatable(species_summary_df(), rownames = FALSE,
                     options = list(pageLength = 8, dom = "tp"),
                     class = "compact stripe hover")
    set_status("Species summary complete", done = TRUE)
    out
  })

  output$map <- renderLeaflet({
    leaflet::leaflet() %>% leaflet::addProviderTiles(leaflet::providers$CartoDB.DarkMatter) %>%
      setView(lng = 0, lat = 15, zoom = 3) %>%
      leaflet.extras::addDrawToolbar(
        polylineOptions     = FALSE,
        circleOptions       = FALSE,
        markerOptions       = FALSE,
        circleMarkerOptions = FALSE,
        rectangleOptions    = FALSE,
        polygonOptions      = leaflet.extras::drawPolygonOptions(shapeOptions = leaflet.extras::drawShapeOptions(fillOpacity = 0.1, color = "#e67e00", weight = 2)),
        editOptions         = leaflet.extras::editToolbarOptions(edit = FALSE, remove = TRUE)
      )
  })

  observe({
    set_status(if (identical(input$mode, "hotspot")) "Rendering hotspot map..." else "Rendering species map...")
    show_progress(5, if (identical(input$mode, "hotspot")) "Rendering hotspot map..." else "Rendering species map...")
    withProgress(
      message = if (identical(input$mode, "hotspot")) "Rendering hotspot map..." else "Rendering species map...",
      value = 0,
      {
        tryCatch({
          proxy <- leafletProxy("map")
          clear_map_empty()
          sel <- input$sel_vars %||% character(0)
          incProgress(0.2)

          show_progress(20, "Checking data...")
          if (input$mode == "single") {
            req(input$species)
            sp_sel <- normalize_species_value(input$species)
            if (!is_valid_species_selection(sp_sel)) {
              exposure_msg_rv(paste0("Species not found: '", input$species, "'."))
              proxy %>% clearImages() %>% clearControls() %>% setView(lng = 0, lat = 15, zoom = 3)
              map_raster_rv(NULL)
              hide_progress()
              set_status("Species render complete (species not found)", done = TRUE)
              return()
            }
            cell_list <- get_species_year_cell_list(sp_sel, input$year, sel)
            if (!length(cell_list)) {
              exposure_msg_rv(paste0("No exposed cells for '", species_display_label(sp_sel),
                                     "' in ", input$year,
                                     " with the selected variable(s)."))
              proxy %>% clearImages() %>% clearControls() %>% setView(lng = 0, lat = 15, zoom = 3)
              show_map_empty(sprintf("No exposed cells for %s in %s.",
                                     species_display_label(sp_sel), input$year))
              hide_progress()
              set_status("Species render complete (no exposed cells)", done = TRUE)
              return()
            }
            layers <- build_species_layers_from_cell_list(cell_list)
            if (!length(layers)) {
              exposure_msg_rv(paste0("No exposed cells for '", species_display_label(sp_sel),
                                     "' in ", input$year,
                                     " with the selected variable(s)."))
              proxy %>% clearImages() %>% clearControls() %>% setView(lng = 0, lat = 15, zoom = 3)
              show_map_empty(sprintf("No exposed cells for %s in %s.",
                                     species_display_label(sp_sel), input$year))
              map_raster_rv(NULL)
              hide_progress()
              set_status("Species render complete (no exposed cells)", done = TRUE)
              return()
            }
            exposure_msg_rv(NULL)
            proxy %>% clearImages() %>% clearControls()
            drawn <- character(0)

            # --- full range background (Mammals + Birds only) ---
            rc <- range_cells_mb[[sp_sel]]
            if (!is.null(rc) && length(rc)) {
              rc <- rc[rc >= 1L & rc <= ncell(tpl)]
              if (length(rc)) {
                rv_bg <- rep(NA_real_, ncell(tpl)); rv_bg[rc] <- 1
                r_bg  <- raster(setValues(tpl, rv_bg))
                xy_bg <- terra::xyFromCell(tpl, rc)
                xbuf  <- max((max(xy_bg[,1]) - min(xy_bg[,1])) * 0.5, 5)
                ybuf  <- max((max(xy_bg[,2]) - min(xy_bg[,2])) * 0.5, 5)
                e_bg  <- terra::ext(
                  max(-180, min(xy_bg[,1]) - xbuf), min(180, max(xy_bg[,1]) + xbuf),
                  max(-61,  min(xy_bg[,2]) - ybuf), min(86,  max(xy_bg[,2]) + ybuf)
                )
                r_bg_crop <- raster(terra::crop(setValues(tpl, rv_bg), e_bg))
                proxy %>% addRasterImage(r_bg_crop, colors = "#a8e063",
                                         opacity = 0.35, layerId = "range_bg",
                                         method = "ngb")
              }
            }

            # store cropped terra rasters for download
            terra_layers <- list()
            for (nm in names(layers)) {
              proxy %>% addRasterImage(layers[[nm]], colors = var_col(nm),
                                       opacity = 0.7, layerId = paste0("lyr_", nm),
                                       method = "ngb")
              drawn <- c(drawn, nm)
              # keep a terra version (full-extent, single layer)
              cv2 <- unique(as.integer(cell_list[[nm]]))
              cv2 <- cv2[!is.na(cv2) & cv2 >= 1 & cv2 <= ncell(tpl)]
              rv2 <- rep(NA_real_, ncell(tpl)); rv2[cv2] <- 1
              terra_layers[[nm]] <- setValues(tpl, rv2)
            }
            if (length(terra_layers) > 0) {
              names(terra_layers) <- make.names(names(terra_layers))
              map_raster_rv(tryCatch(do.call(c, terra_layers), error = function(e) terra_layers[[1]]))
            }
            proxy %>% addLegend(position = "bottomright",
                                colors = vapply(drawn, var_col, character(1)),
                                labels = vapply(drawn, var_label, character(1)),
                                title = "Climate variable", layerId = "legend_main")
            incProgress(0.8)
            show_progress(90, "Fitting map bounds...")
            all_cells <- unique(as.integer(unlist(cell_list, use.names = FALSE)))
            all_cells <- all_cells[!is.na(all_cells) & all_cells >= 1 & all_cells <= ncell(tpl)]
            if (length(all_cells)) {
              xy      <- terra::xyFromCell(tpl, all_cells)
              lng_min <- min(xy[, 1]); lng_max <- max(xy[, 1])
              lat_min <- min(xy[, 2]); lat_max <- max(xy[, 2])
              lng_pad <- max((lng_max - lng_min) * 0.25, 2)
              lat_pad <- max((lat_max - lat_min) * 0.25, 2)
              proxy %>% fitBounds(lng_min - lng_pad, lat_min - lat_pad,
                                  lng_max + lng_pad, lat_max + lat_pad)
            }
            hide_progress()
            set_status("Species map complete", done = TRUE)
          } else if (isTRUE(input$change_mode)) {
            # ── Two-year change map: delta of exposed-species count per cell ──
            baseline_yr <- suppressWarnings(as.integer(input$baseline_year %||% min(years_avail)))
            target_yr   <- suppressWarnings(as.integer(input$year))
            if (identical(baseline_yr, target_yr)) {
              exposure_msg_rv("Baseline and selected year are the same — pick different years to see change.")
              proxy %>% clearImages() %>% clearControls()
              map_raster_rv(NULL)
              show_map_empty("Baseline equals the selected year — pick different years to see change.")
              hide_progress()
              set_status("Change map (baseline equals target)", done = TRUE)
              return()
            }
            rr <- build_change_raster(baseline_yr, target_yr)
            if (is.null(rr)) {
              proxy %>% clearImages() %>% clearControls()
              map_raster_rv(NULL)
              show_map_empty(sprintf("No cells changed between %d and %d.", baseline_yr, target_yr))
              hide_progress()
              set_status("Change map complete (no change)", done = TRUE)
              return()
            }
            exposure_msg_rv(NULL)
            vals <- raster::values(rr); vals <- vals[is.finite(vals)]
            if (!length(vals)) {
              proxy %>% clearImages() %>% clearControls()
              map_raster_rv(NULL)
              show_map_empty(sprintf("No cells changed between %d and %d.", baseline_yr, target_yr))
              hide_progress()
              set_status("Change map complete (no change)", done = TRUE)
              return()
            }
            # Symmetric diverging palette centred on 0
            lim <- max(abs(vals), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1
            div_cols <- colorRampPalette(c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"))(255)
            pal_change <- colorNumeric(div_cols, domain = c(-lim, lim), na.color = "transparent")
            proxy %>% clearImages() %>% clearControls()
            leg_breaks <- pretty(c(-lim, lim), n = 5)
            leg_breaks <- leg_breaks[abs(leg_breaks) <= lim]
            proxy %>%
              addRasterImage(rr, colors = pal_change, opacity = input$hotspot_opacity %||% 0.95,
                             layerId = "hotspot", method = "ngb") %>%
              addLegend(position = "bottomright",
                        colors  = pal_change(leg_breaks),
                        labels  = ifelse(leg_breaks > 0, paste0("+", leg_breaks), as.character(leg_breaks)),
                        opacity = input$hotspot_opacity %||% 0.95,
                        title   = sprintf("&Delta; species<br>%d &rarr; %d", baseline_yr, target_yr),
                        layerId = "legend_main")
            map_raster_rv(terra::rast(rr))
            hide_progress()
            set_status(sprintf("Change map complete (%d \u2192 %d)", baseline_yr, target_yr), done = TRUE)
          } else {
            # Use precomputed raster when defaults are active (all vars, no filters,
            # and the exposure threshold is at the 1% floor). Above 1% the count
            # depends on per-species propExposed, which the precomputed/fast paths
            # (cell_trend_df) do not carry, so route to the exact per-row path.
            excl_at_floor <- ((as.numeric(input$excl_threshold %||% 10)) / 100) <= 0.01
            all_groups_selected <- setequal(input$flt_groups %||% ALL_GROUPS, ALL_GROUPS)
            is_default_hotspot <- (
              excl_at_floor &&
              setequal(sel, vars_avail) &&
              all_groups_selected &&
              !isTRUE(input$flt_threatened) &&
              !isTRUE(input$flt_data_deficient) &&
              !nzchar(paste(input$flt_order   %||% "", collapse = "")) &&
              !nzchar(paste(input$flt_family  %||% "", collapse = ""))
            )
            yr_key <- as.character(input$year)
            if (is_default_hotspot && !is.null(hotspot_raster_cache[[yr_key]])) {
              rr <- hotspot_raster_cache[[yr_key]]
              log_msg(sprintf("Using precomputed hotspot raster for year %s", yr_key))
            } else if (is_default_hotspot) {
              # Fast path: default view for a non-precomputed year, built from the
              # in-memory cell_trend_df instead of a disk read + re-aggregation.
              rr <- default_hotspot_raster(input$year)
              if (is.null(rr)) {
                proxy %>% clearImages() %>% clearControls()
                map_raster_rv(NULL)
                show_map_empty("No exposed cells for the current selection.")
                hide_progress()
                set_status("Hotspot map complete (no exposed cells)", done = TRUE)
                return()
              }
              log_msg(sprintf("Using fast cell_trend_df raster for year %s", yr_key))
            } else {
              d <- filtered_df()
              if (!length(sel) || nrow(d) == 0) {
                exposure_msg_rv("No exposed cells for the selected year/filter/variable combination.")
                proxy %>% clearImages() %>% clearControls()
                show_map_empty("No exposed cells for the selected year, variables, or filters.")
                hide_progress()
                set_status("Map update complete (no exposed cells)", done = TRUE)
                return()
              }
              rr <- build_hotspot_raster(d)
            }
            if (is.null(rr)) {
              proxy %>% clearImages() %>% clearControls()
              map_raster_rv(NULL)
              show_map_empty("No exposed cells for the current selection.")
              hide_progress()
              set_status("Hotspot map complete (no exposed cells)", done = TRUE)
              return()
            }
            exposure_msg_rv(NULL)
            vals <- raster::values(rr)
            vals <- vals[is.finite(vals)]
            if (!length(vals)) {
              proxy %>% clearImages() %>% clearControls()
              map_raster_rv(NULL)
              show_map_empty("No exposed cells for the current selection.")
              hide_progress()
              set_status("Hotspot map complete (no exposed cells)", done = TRUE)
              return()
            }
            minN <- suppressWarnings(min(vals, na.rm = TRUE))
            maxN <- suppressWarnings(max(vals, na.rm = TRUE))
            if (!is.finite(minN) || !is.finite(maxN) || maxN == 0) {
              proxy %>% clearImages() %>% clearControls()
              map_raster_rv(NULL)
              show_map_empty("No exposed cells for the current selection.")
              hide_progress()
              set_status("Hotspot map complete (no exposed cells)", done = TRUE)
              return()
            }
            # Fixed colour scale: lock limits to the 2024 reference range so
            # colours are comparable across years (values outside are clamped by
            # hotspot_transform's pmin/pmax). Fall back to this year's range if
            # the reference is unavailable.
            if (isTRUE(input$fix_color_scale)) {
              fr <- fixed_scale_range()
              if (!is.null(fr) && is.finite(fr[1]) && is.finite(fr[2]) && fr[2] > fr[1]) {
                minN <- fr[1]; maxN <- fr[2]
              }
            }
            if (minN == maxN) {
              minN <- minN - 0.5
              maxN <- maxN + 0.5
            }
            n_cols <- max(2, min(100, as.integer(maxN - minN + 1)))
            cm_cols <- colorRampPalette(c("steelblue4", "steelblue1", "gold", "red1", "red4"))(n_cols)
            pal_base <- colorNumeric(cm_cols, domain = c(0, 1), na.color = "transparent")
            pal <- function(x) pal_base(hotspot_transform(x, minN, maxN))
            proxy %>% clearImages() %>% clearControls()
            # build explicit legend colours (addLegend requires a proper colorNumeric
            # palette object; using colors+labels avoids that restriction)
            n_leg        <- min(7L, length(vals))
            leg_breaks   <- seq(minN, maxN, length.out = n_leg)
            leg_cols     <- pal(leg_breaks)
            leg_labels   <- formatC(round(leg_breaks), format = "d", big.mark = ",")
            leg_title <- sprintf("# Species exposed<br><span style='font-weight:normal;font-size:10px'>&ge; %d%% of range%s</span>",
                                 as.integer(as.numeric(input$excl_threshold %||% 10)),
                                 if (isTRUE(input$fix_color_scale)) "<br>fixed scale (2024)" else "")
            proxy %>%
              addRasterImage(rr, colors = pal, opacity = input$hotspot_opacity %||% 0.95, layerId = "hotspot", method = "ngb") %>%
              addLegend(position = "bottomright",
                        colors  = leg_cols,
                        labels  = leg_labels,
                        opacity = input$hotspot_opacity %||% 0.95,
                        title   = leg_title,
                        layerId = "legend_main")
            # store for download
            map_raster_rv(terra::rast(rr))
            hide_progress()
            set_status("Hotspot map complete", done = TRUE)
          }
          incProgress(1)
        }, error = function(e) {
          exposure_msg_rv(paste("Map update failed:", conditionMessage(e)))
          leafletProxy("map") %>% clearImages() %>% clearControls()
          hide_progress()
          set_status("Map update failed")
        })
      }
    )
  })

  # ── Summary stat card above the map ──────────────────────────────────────
  summary_stats <- reactive({
    sel <- input$sel_vars %||% character(0)
    if (input$mode == "single") {
      sp_sel <- normalize_species_value(input$species)
      if (!length(sel) || !is_valid_species_selection(sp_sel)) return(NULL)
      cl <- get_species_year_cell_list(sp_sel, input$year, sel)
      cells <- unique(as.integer(unlist(cl, use.names = FALSE)))
      cells <- cells[!is.na(cells)]
      if (!length(cells)) return(NULL)
      list(kind = "single",
           n_cells = length(cells),
           n_species = 1L,
           species = species_display_label(sp_sel),
           n_vars = length(cl))
    } else {
      d <- tryCatch(filtered_df(), error = function(e) NULL)
      if (is.null(d) || !nrow(d)) return(NULL)
      top_order <- d %>% filter(!is.na(orderName)) %>%
        group_by(orderName) %>% summarize(n = n_distinct(spName), .groups = "drop") %>%
        arrange(desc(n)) %>% slice(1)
      top_family <- d %>% filter(!is.na(familyName)) %>%
        group_by(familyName) %>% summarize(n = n_distinct(spName), .groups = "drop") %>%
        arrange(desc(n)) %>% slice(1)
      list(kind = "hotspot",
           n_cells = dplyr::n_distinct(d$cell[!is.na(d$cell)]),
           n_species = dplyr::n_distinct(d$spName),
           top_order = if (nrow(top_order)) sprintf("%s (%d)", top_order$orderName, top_order$n) else "\u2014",
           top_family = if (nrow(top_family)) sprintf("%s (%d)", top_family$familyName, top_family$n) else "\u2014")
    }
  })

  output$summary_card <- renderUI({
    s <- summary_stats()
    fmt <- function(x) formatC(x, format = "d", big.mark = ",")
    stat <- function(num, lab, wide = FALSE, tip = NULL)
      div(class = if (wide) "stat wide" else "stat",
          span(class = "num", num),
          span(class = "lab", lab,
               if (!is.null(tip)) tags$span("?", class = "tip", `data-tip` = tip)))
    if (is.null(s)) {
      return(div(class = "summary-card",
                 div(class = "stat", span(class = "num", "\u2014"),
                     span(class = "lab", "No exposed cells in view"))))
    }
    if (identical(s$kind, "single")) {
      div(class = "summary-card",
          stat(fmt(s$n_cells), "Exposed cells", tip =
               "Number of distinct grid cells where this species is exposed in the selected year, counting a cell once even if several climate variables trip it."),
          stat(fmt(s$n_vars), "Variables exposed", tip =
               "How many of the selected climate variables expose this species in at least one cell this year."),
          stat(s$species, "Species", wide = TRUE, tip =
               "The species currently selected in single-species mode."))
    } else {
      sweep <- summary_sweep()
      sweep_ui <- if (!is.null(sweep) && nrow(sweep)) {
        parts <- sprintf("%d%%: %s", sweep$thr, formatC(sweep$n, format = "d", big.mark = ","))
        div(style = "font-size:11px;color:#666;margin:-2px 0 6px;text-align:center;",
            tags$span("Exposed species by threshold — ", style = "color:#888;"),
            paste(parts, collapse = "  \u00b7  "),
            tags$span("?", class = "tip",
              `data-tip` = "How many species would be on the map at other exposure thresholds for this year and variable selection (all taxa; order/family/threatened/DD filters not applied here). Shows the sensitivity of the map to the threshold slider."))
      } else NULL
      tagList(
        div(class = "summary-card",
            stat(fmt(s$n_cells), "Exposed cells", tip =
                 "Number of distinct grid cells drawn on the map — cells where at least one qualifying species is exposed this year (after the exposure-threshold and any taxon filters). Each cell counted once."),
            stat(fmt(s$n_species), "Exposed species", tip =
                 "Number of distinct species contributing to the map this year: those passing the exposure threshold (>= the slider % of their range exposed) and any active order/family/threatened/DD filters."),
            stat(s$top_order, "Top order", wide = TRUE, tip =
                 "The taxonomic order contributing the most distinct species to the map this year, with that species count in parentheses."),
            stat(s$top_family, "Top family", wide = TRUE, tip =
                 "The taxonomic family contributing the most distinct species to the map this year, with that species count in parentheses.")),
        sweep_ui
      )
    }
  })

  # Year shown as a title above the map (hotspot mode: the map year; single
  # mode: the selected year and species context is in the summary card).
  output$map_year_title <- renderUI({
    yr <- input$year
    if (is.null(yr) || !nzchar(as.character(yr))) return(NULL)
    lab <- if (isTRUE(input$change_mode) && identical(input$mode, "hotspot"))
      sprintf("%s \u2192 %s", input$baseline_year %||% "", yr) else as.character(yr)
    div(style = "text-align:center;font-size:22px;font-weight:700;color:#333;margin:2px 0 4px;",
        lab)
  })

  # Total exposed species per year at the current threshold, from sp_trend_df
  # (mean_prop = per-species-year exposed fraction). Threshold-aware; taxon
  # subfilters are NOT applied (sp_trend_df carries no taxonomy) — noted in UI.
  # Memoized per threshold so animation does not recompute each frame.
  .trend_ts_memo <- new.env(parent = emptyenv())
  trend_series <- reactive({
    if (is.null(sp_trend_df) || !nrow(sp_trend_df)) return(NULL)
    thr <- (as.numeric(input$excl_threshold %||% 10)) / 100
    key <- as.character(thr)
    if (!is.null(.trend_ts_memo[[key]])) return(.trend_ts_memo[[key]])
    d <- sp_trend_df
    if (is.finite(thr) && thr > 0.01) d <- d[!is.na(d$mean_prop) & d$mean_prop >= thr, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    ts <- d %>% group_by(year) %>% summarize(n_sp = dplyr::n_distinct(spName), .groups = "drop") %>%
      arrange(year)
    .trend_ts_memo[[key]] <- ts
    ts
  })

  output$trend_ts <- renderPlot({
    ts <- trend_series()
    if (is.null(ts) || !nrow(ts)) {
      op <- par(mar = c(0,0,0,0)); on.exit(par(op))
      plot.new(); text(0.5, 0.5, "No trend available", col = "#888"); return(invisible())
    }
    yr <- suppressWarnings(as.integer(input$year))
    op <- par(mar = c(2.6, 4.0, 0.4, 0.6), mgp = c(2.3, 0.5, 0), cex = 0.85, tcl = -0.3)
    on.exit(par(op))
    plot(ts$year, ts$n_sp, type = "l", lwd = 2, col = "#c0392b",
         xlab = "", ylab = "Exposed species", xaxs = "i", las = 1,
         panel.first = grid(col = "#eee", lty = 1))
    if (is.finite(yr)) {
      abline(v = yr, col = "#2166ac", lwd = 1.5, lty = 2)
      yv <- ts$n_sp[match(yr, ts$year)]
      if (!is.na(yv)) points(yr, yv, pch = 19, col = "#2166ac", cex = 1.1)
    }
  })

  # Threshold-sweep for the summary card: exposed-species count at a few
  # thresholds for the current year (all taxa). Cheap: uses filtered year slice.
  summary_sweep <- reactive({
    if (!identical(input$mode, "hotspot")) return(NULL)
    d <- tryCatch(year_df(), error = function(e) NULL)
    if (is.null(d) || !nrow(d) || !("propExposed" %in% names(d))) return(NULL)
    sel <- input$sel_vars %||% character(0)
    if (length(sel)) d <- d[d$var %in% sel, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    thr_grid <- c(1, 5, 10, 25)
    vapply(thr_grid, function(p) {
      dd <- d[is.na(d$propExposed) | d$propExposed >= p/100, , drop = FALSE]
      length(unique(dd$spName))
    }, integer(1)) -> counts
    data.frame(thr = thr_grid, n = counts)
  })

  # ── Year animation (hotspot mode): play/pause auto-advance ───────────────
  playing_rv <- reactiveVal(FALSE)
  observeEvent(input$play_years, { playing_rv(TRUE) })
  observeEvent(input$pause_years, { playing_rv(FALSE) })
  # Stop playback when leaving hotspot mode
  observeEvent(input$mode, { if (!identical(input$mode, "hotspot")) playing_rv(FALSE) }, ignoreInit = TRUE)
  observe({
    if (!isTRUE(playing_rv())) return()
    # Slider is in seconds/year; convert to ms. Reading it here (non-isolated)
    # makes a speed change take effect on the next scheduled frame.
    interval_ms <- as.numeric(input$play_speed %||% 0.75) * 1000
    if (!is.finite(interval_ms) || interval_ms < 100) interval_ms <- 750
    invalidateLater(interval_ms, session)
    isolate({
      cur <- suppressWarnings(as.integer(input$year))
      idx <- match(cur, years_avail)
      if (is.na(idx)) idx <- length(years_avail)
      nxt <- if (idx >= length(years_avail)) years_avail[1] else years_avail[idx + 1]
      updateSelectInput(session, "year", selected = nxt)
    })
  })

  click_data <- reactiveVal(NULL)
  click_info <- reactiveVal(list(cell = NA_integer_, lat = NA, lng = NA))
  hotspot_diag_data <- reactiveVal(NULL)
  polygon_cells_rv <- reactiveVal(NULL)  # cell IDs inside drawn polygon (year-independent)
  output$polygon_active <- renderText({ if (!is.null(polygon_cells_rv())) "yes" else "no" })
  outputOptions(output, "polygon_active", suspendWhenHidden = FALSE)
  # Table recomputes reactively whenever year/vars change while polygon is drawn
  polygon_table_data <- reactive({
    cell_ids <- polygon_cells_rv()
    req(!is.null(cell_ids))
    if (length(cell_ids) == 0) return(data.frame())
    fd <- filtered_df()
    d  <- fd %>% filter(cell %in% cell_ids)
    d  <- ensure_prop_exposed(d)
    if (nrow(d) == 0) return(data.frame())
    pd <- d %>%
      dplyr::select(spName, var,
                    any_of(c("propExposed")),
                    any_of(c("orderName", "familyName", "redlistCategory"))) %>%
      mutate(var = vapply(as.character(var), var_label, character(1))) %>%
      distinct() %>%
      arrange(desc(propExposed), var, spName) %>%
      rename(Species = spName, Variable = var)
    if ("propExposed"     %in% names(pd)) pd <- rename(pd, `Proportion exposed` = propExposed)
    if ("orderName"       %in% names(pd)) pd <- rename(pd, Order = orderName)
    if ("familyName"      %in% names(pd)) pd <- rename(pd, Family = familyName)
    if ("redlistCategory" %in% names(pd)) pd <- rename(pd, `IUCN status` = redlistCategory)
    pd
  })

  # Regional summary card for the selected area (current year, current filters).
  output$polygon_summary <- renderUI({
    cell_ids <- polygon_cells_rv()
    if (is.null(cell_ids) || !length(cell_ids)) return(NULL)
    fd <- tryCatch(filtered_df(), error = function(e) NULL)
    if (is.null(fd) || !nrow(fd)) return(div(class = "summary-card",
      div(class = "stat", span(class = "num", "\u2014"), span(class = "lab", "No exposed species in area"))))
    d <- fd %>% filter(cell %in% cell_ids)
    if (!nrow(d)) return(div(class = "summary-card",
      div(class = "stat", span(class = "num", "\u2014"), span(class = "lab", "No exposed species in area"))))
    n_sp   <- dplyr::n_distinct(d$spName)
    n_cell <- dplyr::n_distinct(d$cell)
    top_ord <- if ("orderName" %in% names(d)) {
      t <- d %>% distinct(spName, orderName) %>% count(orderName, sort = TRUE)
      if (nrow(t)) sprintf("%s (%d)", t$orderName[1], t$n[1]) else "\u2014"
    } else "\u2014"
    stat <- function(num, lab, wide = FALSE)
      div(class = if (wide) "stat wide" else "stat", span(class = "num", num), span(class = "lab", lab))
    div(class = "summary-card",
        stat(formatC(n_cell, format = "d", big.mark = ","), "Cells in area"),
        stat(formatC(n_sp,   format = "d", big.mark = ","), "Exposed species"),
        stat(top_ord, "Top order", wide = TRUE))
  })

  # Regional per-year trend: total exposed species within the selected cells per
  # year, from cell_trend_df (n_sp per cell/year). Reflects the default (unfiltered)
  # species count restricted to the area; taxon/threshold filters are not applied
  # here (cell_trend_df carries no per-species detail) — caption notes this.
  output$polygon_trend <- renderPlot({
    cell_ids <- polygon_cells_rv()
    if (is.null(cell_ids) || !length(cell_ids) || is.null(cell_trend_df)) {
      op <- par(mar = c(0,0,0,0)); on.exit(par(op)); plot.new(); text(0.5,0.5,"No area selected", col="#888"); return(invisible())
    }
    ct <- cell_trend_df[cell_trend_df$cell %in% cell_ids, , drop = FALSE]
    if (!nrow(ct)) { op <- par(mar=c(0,0,0,0)); on.exit(par(op)); plot.new(); text(0.5,0.5,"No exposure in area", col="#888"); return(invisible()) }
    ts <- ct %>% group_by(year) %>% summarize(n_sp = sum(n_sp), .groups = "drop") %>% arrange(year)
    yr <- suppressWarnings(as.integer(input$year))
    op <- par(mar = c(2.6, 4.0, 0.4, 0.6), mgp = c(2.3, 0.5, 0), cex = 0.85, tcl = -0.3); on.exit(par(op))
    plot(ts$year, ts$n_sp, type = "l", lwd = 2, col = "#c0392b", xlab = "", ylab = "Exposed species\u00d7cells",
         xaxs = "i", las = 1, panel.first = grid(col = "#eee", lty = 1))
    if (is.finite(yr)) {
      abline(v = yr, col = "#2166ac", lwd = 1.5, lty = 2)
      yv <- ts$n_sp[match(yr, ts$year)]; if (!is.na(yv)) points(yr, yv, pch = 19, col = "#2166ac", cex = 1.1)
    }
  })
  # Tracks the last raster(s) drawn on the map for download
  map_raster_rv <- reactiveVal(NULL)
  adm1_cache <- reactiveVal(NULL)
  adm2_cache <- reactiveVal(NULL)

  # Clear stale click/polygon state whenever the user switches mode
  observeEvent(input$mode, {
    click_info(list(cell = NA_integer_, lat = NA, lng = NA))
    click_data(NULL)
    hotspot_diag_data(NULL)
    polygon_cells_rv(NULL)
    polygon_empty_reason_rv(NULL)
    map_raster_rv(NULL)
    # Remove selected polygon overlay from the map
    leafletProxy("map") %>% clearGroup("drawn")
  }, ignoreInit = TRUE)

  observeEvent(input$pol_country, {
    iso <- toupper(trimws(input$pol_country %||% ""))
    polygon_cells_rv(NULL)
    polygon_empty_reason_rv(NULL)
    leafletProxy("map") %>% clearGroup("drawn")
    if (!nzchar(iso)) {
      adm1_cache(NULL)
      adm2_cache(NULL)
      updateSelectizeInput(session, "pol_state", choices = character(0), selected = "", server = TRUE)
      updateSelectizeInput(session, "pol_county", choices = character(0), selected = "", server = TRUE)
      return()
    }
    set_status(sprintf("Loading administrative units for %s...", iso))
    tryCatch({
      v1 <- get_gadm_cached(iso, 1L)
      d1 <- as.data.frame(v1)
      adm1_cache(d1)
      s1 <- sort(unique(as.character(d1$NAME_1)))
      s1 <- s1[!is.na(s1) & nzchar(s1)]
      updateSelectizeInput(session, "pol_state",
                           choices = c(list("(All states/provinces)" = ""), setNames(as.list(s1), s1)),
                           selected = "", server = TRUE)

      v2 <- tryCatch(get_gadm_cached(iso, 2L), error = function(e) NULL)
      if (!is.null(v2)) {
        d2 <- as.data.frame(v2)
        adm2_cache(d2)
      } else {
        adm2_cache(NULL)
      }
      updateSelectizeInput(session, "pol_county", choices = character(0), selected = "", server = TRUE)
      set_status(sprintf("Administrative units ready for %s", iso), done = TRUE)
    }, error = function(e) {
      adm1_cache(NULL)
      adm2_cache(NULL)
      updateSelectizeInput(session, "pol_state", choices = character(0), selected = "", server = TRUE)
      updateSelectizeInput(session, "pol_county", choices = character(0), selected = "", server = TRUE)
      exposure_msg_rv(paste("Could not load political units:", conditionMessage(e)))
      set_status("Administrative unit load failed")
    })
  }, ignoreInit = TRUE)

  observeEvent(input$pol_state, {
    d2 <- adm2_cache()
    st <- input$pol_state %||% ""
    if (is.null(d2) || !nrow(d2) || !nzchar(st)) {
      updateSelectizeInput(session, "pol_county",
                           choices = list("(All counties/districts)" = ""),
                           selected = "", server = TRUE)
      return()
    }
    c2 <- d2 %>%
      filter(NAME_1 == st) %>%
      pull(NAME_2) %>%
      as.character() %>%
      unique() %>%
      sort()
    c2 <- c2[!is.na(c2) & nzchar(c2)]
    updateSelectizeInput(session, "pol_county",
                         choices = c(list("(All counties/districts)" = ""), setNames(as.list(c2), c2)),
                         selected = "", server = TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$load_pol_unit, {
    req(input$mode == "hotspot")
    iso <- toupper(trimws(input$pol_country %||% ""))
    req(nzchar(iso))
    set_status("Loading selected boundary...")
    show_progress(10, "Fetching boundary...")
    tryCatch({
      st <- trimws(input$pol_state %||% "")
      ct <- trimws(input$pol_county %||% "")

      unit_v <- NULL
      label <- iso
      unit_key <- NULL
      if (nzchar(ct)) {
        show_progress(20, "Loading ADM2 boundary...")
        v2 <- get_gadm_cached(iso, 2L)
        d2 <- as.data.frame(v2)
        idx <- which(d2$NAME_1 == st & d2$NAME_2 == ct)
        if (!length(idx)) stop("Selected county/district not found")
        unit_v <- v2[idx, ]
        label <- paste0(ct, ", ", st, " (", iso, ")")
        unit_key <- paste0("ADM2|", iso, "|", st, "|", ct)
      } else if (nzchar(st)) {
        show_progress(20, "Loading state boundary...")
        v1 <- get_gadm_cached(iso, 1L)
        d1 <- as.data.frame(v1)
        idx <- which(d1$NAME_1 == st)
        if (!length(idx)) stop("Selected state/province not found")
        unit_v <- v1[idx, ]
        label <- paste0(st, " (", iso, ")")
        unit_key <- paste0("ADM1|", iso, "|", st)
      } else {
        show_progress(20, "Loading country boundary...")
        unit_v <- get_gadm_cached(iso, 0L)
        label <- iso
        unit_key <- paste0("ADM0|", iso)
      }

      show_progress(55, "Finding exposed cells in boundary...")
      if (!length(avail_cells_all)) {
        polygon_cells_rv(NULL)
        polygon_empty_reason_rv("No valid raster cells are available to intersect.")
        hide_progress()
        set_status("No available cells for selection")
        return()
      }
      cell_ids <- intersect_cached(unit_key, unit_v)
      show_progress(85, "Drawing boundary on map...")
      if (!length(cell_ids)) {
        polygon_cells_rv(integer(0))
        polygon_empty_reason_rv("No exposed raster cells intersect this selected boundary.")
        hide_progress()
        set_status("Boundary loaded, but no exposed cells intersect it", done = TRUE)
      } else {
        polygon_cells_rv(cell_ids)
        polygon_empty_reason_rv(NULL)
        hide_progress()
        set_status(sprintf("Boundary loaded: %s (%d cells)", label, length(cell_ids)), done = TRUE)
      }

      gj <- spatvector_to_leaflet_geojson(unit_v)
      leafletProxy("map") %>%
        clearGroup("drawn") %>%
        addGeoJSON(gj,
                   group = "drawn",
                   color = "#e67e00",
                   weight = 2,
                   fillColor = "#e67e00",
                   fillOpacity = 0.15)
      ext <- terra::ext(unit_v)
      leafletProxy("map") %>% fitBounds(ext[1], ext[3], ext[2], ext[4])

      click_data(NULL)
      hotspot_diag_data(NULL)
      click_info(list(cell = NA_integer_, lat = NA, lng = NA))
    }, error = function(e) {
      polygon_cells_rv(NULL)
      hide_progress()
      exposure_msg_rv(paste("Boundary selection failed:", conditionMessage(e)))
      set_status("Boundary selection failed")
    })
  })

  # Handle freehand polygon drawn by user on the map
  observeEvent(input$map_draw_new_feature, ignoreNULL = TRUE, {
    feat <- input$map_draw_new_feature
    coords_raw <- tryCatch(feat$geometry$coordinates[[1]], error = function(e) NULL)
    if (is.null(coords_raw) || length(coords_raw) < 3) return()

    coord_mat <- tryCatch({
      do.call(rbind, lapply(coords_raw, function(pt) c(as.numeric(pt[[1]]), as.numeric(pt[[2]]))))
    }, error = function(e) NULL)
    if (is.null(coord_mat) || nrow(coord_mat) < 3) return()

    # Close ring if needed
    if (!identical(coord_mat[1, ], coord_mat[nrow(coord_mat), ])) {
      coord_mat <- rbind(coord_mat, coord_mat[1, ])
    }

    poly_v <- tryCatch(
      terra::vect(list(list(coord_mat)), type = "polygons", crs = terra::crs(tpl)),
      error = function(e) NULL
    )
    if (is.null(poly_v)) return()

    avail_cells <- avail_cells_all
    if (!length(avail_cells)) { polygon_cells_rv(NULL); polygon_empty_reason_rv("No valid raster cells are available to intersect."); return() }

    pts <- avail_pts
    inside <- tryCatch(
      terra::is.related(pts, poly_v, "intersects"),
      error = function(e) rep(FALSE, length(avail_cells))
    )
    cell_ids <- avail_cells[inside]
    polygon_cells_rv(if (length(cell_ids)) cell_ids else integer(0))
    if (length(cell_ids) == 0) {
      polygon_empty_reason_rv("No exposed raster cells intersect the drawn polygon.")
    } else {
      polygon_empty_reason_rv(NULL)
    }

    click_data(NULL)
    hotspot_diag_data(NULL)
    click_info(list(cell = NA_integer_, lat = NA, lng = NA))
    set_status(sprintf("Drawn polygon: %d cells selected", length(cell_ids)), done = TRUE)
  })

  observeEvent(input$map_draw_deleted_features, {
    polygon_cells_rv(NULL)
    polygon_empty_reason_rv(NULL)
    set_status("Drawn polygon cleared", done = TRUE)
  })

  observeEvent(input$map_click, ignoreNULL = TRUE, {
    set_status("Loading clicked cell...")
    show_progress(5, "Reading click location...")
    withProgress(message = "Loading clicked cell...", value = 0, {
      tryCatch({
        click <- input$map_click
        if (is.null(click$lng) || is.null(click$lat) || !is.finite(click$lng) || !is.finite(click$lat)) {
          click_data(NULL)
          hotspot_diag_data(NULL)
          click_info(list(cell = NA_integer_, lat = NA, lng = NA))
          set_status("Clicked cell complete (invalid location)", done = TRUE)
          return()
        }
        cn <- tryCatch({
          as.integer(terra::cellFromXY(tpl, matrix(c(click$lng, click$lat), ncol = 2))[1])
        }, error = function(e) NA_integer_)
        if (is.na(cn)) {
          click_data(NULL)
          hotspot_diag_data(NULL)
          click_info(list(cell = NA_integer_, lat = NA, lng = NA))
          set_status("Clicked cell complete (outside raster)", done = TRUE)
          return()
        }

        click_info(list(cell = cn, lat = round(click$lat, 3), lng = round(click$lng, 3)))
        sel <- input$sel_vars %||% character(0)
        if (!length(sel)) {
          click_data(NULL)
          hotspot_diag_data(NULL)
          set_status("Clicked cell complete (no variables selected)", done = TRUE)
          return()
        }

        show_progress(30, "Fetching cell rows...")
        dsrc <- if (input$mode == "hotspot") filtered_df() else year_df()
        dsrc <- ensure_prop_exposed(dsrc)
        incProgress(0.4)
        show_progress(60, "Building species table...")
        d <- dsrc %>%
          filter(cell == cn, var %in% sel) %>%
          dplyr::select(spName, var,
                        any_of(c("propExposed")),
                        any_of(c("orderName", "familyName", "redlistCategory"))) %>%
          mutate(var = vapply(as.character(var), var_label, character(1))) %>%
          distinct()

        if (input$mode == "hotspot" && "propExposed" %in% names(d)) {
          d <- d %>% arrange(desc(propExposed), var, spName)
        } else {
          d <- d %>% arrange(var, spName)
        }

        d <- d %>%
          rename(Species = spName, Variable = var)
        if ("propExposed" %in% names(d)) d <- rename(d, `Proportion exposed` = propExposed)
        if ("orderName"       %in% names(d)) d <- rename(d, Order         = orderName)
        if ("familyName"      %in% names(d)) d <- rename(d, Family        = familyName)
        if ("redlistCategory" %in% names(d)) d <- rename(d, `IUCN status` = redlistCategory)
        click_data(if (!nrow(d)) NULL else d)

        if (input$mode == "hotspot") {
          diag_src <- dsrc %>%
            filter(cell == cn) %>%
            ensure_prop_exposed()

          diag_df <- diag_src %>%
            group_by(orderName, familyName) %>%
            summarize(
              Species_count     = n_distinct(spName),
              Mean_prop_exposed = round(mean(.data$propExposed, na.rm = TRUE), 4),
              .groups = "drop"
            ) %>%
            arrange(desc(Species_count), desc(Mean_prop_exposed))
          hotspot_diag_data(diag_df)
        } else {
          hotspot_diag_data(NULL)
        }
        incProgress(1)
        hide_progress()
        set_status("Clicked cell complete", done = TRUE)
      }, error = function(e) {
        click_data(NULL)
        hotspot_diag_data(NULL)
        exposure_msg_rv(paste("Cell click failed:", conditionMessage(e)))
        hide_progress()
        set_status("Clicked cell failed")
      })
    })
  })

  output$hotspot_diag <- renderDT({
    req(input$mode == "hotspot")
    req(hotspot_diag_data())
    datatable(hotspot_diag_data(), rownames = FALSE,
              options = list(pageLength = 8, dom = "tp"),
              class = "compact stripe hover")
  })

  output$click_title <- renderText({
    info <- click_info()
    if (is.na(info$cell)) return("Click a cell on the map to see exposed species")
    nd <- if (is.null(click_data())) 0 else nrow(click_data())
    sprintf("Cell %d  (%.3f°N, %.3f°E)  —  %d species × variable rows  |  year %s",
            info$cell, info$lat, info$lng, nd, input$year)
  })

  output$click_table <- renderDT({
    req(click_data())
    datatable(click_data(), rownames = FALSE, filter = "top",
              options = list(pageLength = 30, dom = "ftp", scrollX = TRUE,
                             columnDefs = list(
                               list(width = "220px", targets = 1),  # Variable
                               list(width = "72px",  targets = 2)   # Proportion exposed
                             )),
              class = "compact stripe hover")
  })

  output$cell_trend <- renderPlot({
    req(input$mode == "hotspot")
    info <- click_info()
    if (is.na(info$cell)) {
      plot.new(); title("Click a map cell to see trend"); return()
    }
    sel <- input$sel_vars
    if (is.null(sel) || !length(sel)) {
      plot.new(); title("Select at least one climate variable"); return()
    }

    # Use precomputed table: unique species per cell per year (all vars, no double-count)
    # Fill all years in data range with 0 for missing
    all_years_df <- data.frame(year = seq(min(years_avail), max(years_avail)))
    trend <- all_years_df %>%
      dplyr::left_join(
        cell_trend_df %>% filter(cell == info$cell) %>% rename(exposed_species = n_sp),
        by = "year"
      ) %>%
      mutate(exposed_species = ifelse(is.na(exposed_species), 0L, as.integer(exposed_species)))

    if (all(trend$exposed_species == 0)) {
      plot.new(); title("No trend data for current selection"); return()
    }

    par(mar = c(3.5, 4.5, 2.5, 1))
    barplot(trend$exposed_species,
            names.arg = trend$year,
            col    = "steelblue4",
            border = NA,
            xlab   = "Year",
            ylab   = "Exposed species",
            main   = paste0("Species exposed over time — cell ", info$cell),
            las    = 1)
    grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
  })

  output$species_range_trend <- renderPlot({
    req(input$mode == "single")
    sp <- normalize_species_value(input$species)
    req(nzchar(sp))
    sel <- input$sel_vars
    req(length(sel))

    # Use precomputed sp_trend_df: mean propExposed per species×var×year
    # Fill full year range with 0 for missing years
    all_years_df <- data.frame(year = seq(min(years_avail), max(years_avail)))
    trend_raw <- sp_trend_df %>%
      filter(spName == sp, var %in% sel) %>%
      group_by(year) %>%
      summarize(mean_prop_exposed = mean(mean_prop, na.rm = TRUE), .groups = "drop")
    trend <- all_years_df %>%
      dplyr::left_join(trend_raw, by = "year") %>%
      mutate(mean_prop_exposed = ifelse(is.na(mean_prop_exposed), 0, mean_prop_exposed))

    has_data <- trend %>% filter(mean_prop_exposed > 0)
    if (!nrow(has_data)) {
      plot.new()
      title(paste0("No exposure data for \"", sp, "\""))
      return()
    }

    yvals <- trend$mean_prop_exposed
    yvals[is.na(yvals)] <- 0
    ymax  <- max(yvals, na.rm = TRUE)
    par(mar = c(3.5, 4.5, 2.5, 1))
    barplot(yvals,
            names.arg = trend$year,
            col    = "steelblue4",
            border = NA,
            xlab   = "Year",
            ylab   = "Proportion of range exposed",
            main   = paste0(sp, "  \u2014  proportion of range exposed over time"),
            las    = 1,
            ylim   = c(0, max(ymax * 1.15, 0.01)))
    grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
  })

  filtered_export_df <- reactive({
    if (input$mode == "single") {
      # All years for this species + selected vars — enough to reproduce both
      # the summary table (group by year) and the trend plot (group by year x var).
      req(input$species)
      sel <- input$sel_vars %||% character(0)
      req(length(sel))
      load_species_history(normalize_species_value(input$species), sel) %>%
        arrange(year, var) %>%
        mutate(var = vapply(as.character(var), var_label, character(1)))
    } else {
      filtered_df() %>%
        mutate(var = vapply(as.character(var), var_label, character(1)))
    }
  })

  output$download_filtered <- downloadHandler(
    filename = function() {
      if (input$mode == "single")
        paste0("species_", gsub("[^A-Za-z0-9_]", "_", input$species %||% "unknown"), "_allyears.csv")
      else
        paste0("filtered_exposure_hotspot_", input$year, ".csv")
    },
    content = function(file) write.csv(filtered_export_df(), file, row.names = FALSE)
  )

  output$download_template <- downloadHandler(
    filename = function() "landTemplate.tif",
    content  = function(file) file.copy(template_path, file)
  )

  output$download_rscript <- downloadHandler(
    filename = function() "make_species_raster.R",
    content  = function(file) {
      sp  <- gsub("[^A-Za-z0-9_]", "_", input$species %||% "species_name")
      cat(file = file,
'# Reconstruct exposure rasters from the downloaded species CSV\n',
'# Requires: terra\n\n',
'library(terra)\n\n',
'# ---- 1. Load files ----\n',
paste0('d   <- read.csv("species_', sp, '_allyears.csv")\n'),
'tpl <- rast("landTemplate.tif")\n\n',
'# ---- 2. Choose year and variable ----\n',
'yr  <- max(d$year)          # or set to e.g. 2025\n',
'v   <- unique(d$var)[1]     # or any variable label in the CSV\n',
'd_sub <- d[d$year == yr & d$var == v, ]\n\n',
'# ---- 3. Build a raster (1 = exposed cell) ----\n',
'rv    <- rep(NA_real_, ncell(tpl))\n',
'rv[d_sub$cell] <- 1\n',
'r     <- setValues(tpl, rv)\n',
'names(r) <- paste(v, yr)\n\n',
'# ---- 4. Quick plot ----\n',
'plot(r, main = paste(v, yr))\n\n',
'# ---- 5. Save to GeoTIFF ----\n',
paste0('writeRaster(r, "', sp, '_", yr, "_", gsub(" ", "_", v), ".tif", overwrite = TRUE)\n'),
sep = "")
    }
  )

  # --- raster download (shared for both modes) ---
  .build_map_raster_for_download <- function() {
    r <- map_raster_rv()
    if (!is.null(r)) return(r)
    # fallback: rebuild from current state
    if (input$mode == "single" && nzchar(input$species %||% "")) {
      sel <- input$sel_vars %||% character(0)
      cell_list <- get_species_year_cell_list(input$species, input$year, sel)
      if (!length(cell_list)) return(NULL)
      layers <- list()
      for (v in intersect(sel, names(cell_list))) {
        cv <- unique(as.integer(cell_list[[v]]))
        cv <- cv[!is.na(cv) & cv >= 1 & cv <= ncell(tpl)]
        if (!length(cv)) next
        rv <- rep(NA_real_, ncell(tpl)); rv[cv] <- 1
        layers[[make.names(v)]] <- setValues(tpl, rv)
      }
      if (!length(layers)) return(NULL)
      tryCatch(do.call(c, layers), error = function(e) layers[[1]])
    } else if (input$mode == "hotspot") {
      d <- filtered_df()
      if (!nrow(d)) return(NULL)
      rr <- build_hotspot_raster(d)
      if (is.null(rr)) return(NULL)
      terra::rast(rr)
    } else NULL
  }

  .map_raster_filename <- function() {
    if (input$mode == "single")
      paste0("raster_", gsub("[^A-Za-z0-9_]", "_", input$species %||% "species"), "_", input$year, ".tif")
    else
      paste0("raster_hotspot_", input$year, ".tif")
  }
  output$download_map_raster <- downloadHandler(
    filename = .map_raster_filename,
    content = function(file) {
      r <- .build_map_raster_for_download()
      if (is.null(r)) { file.create(file); return() }
      suppressWarnings(terra::writeRaster(r, file, overwrite = TRUE))
    }
  )
  # Download button placed directly below the map (same raster as shown).
  output$download_map_below <- downloadHandler(
    filename = .map_raster_filename,
    content = function(file) {
      r <- .build_map_raster_for_download()
      if (is.null(r)) { file.create(file); return() }
      suppressWarnings(terra::writeRaster(r, file, overwrite = TRUE))
    }
  )
  # Summary stats (the card numbers + current control state) as a one-row CSV.
  output$download_summary <- downloadHandler(
    filename = function() sprintf("exposure_summary_%s_%s.csv",
                                  input$mode %||% "view", input$year %||% "yr"),
    content = function(file) {
      s <- tryCatch(summary_stats(), error = function(e) NULL)
      row <- list(
        mode                 = input$mode %||% "",
        year                 = input$year %||% "",
        variables            = paste(input$sel_vars %||% character(0), collapse = ";"),
        exposure_threshold_pct = as.numeric(input$excl_threshold %||% 10),
        only_threatened      = isTRUE(input$flt_threatened),
        only_data_deficient  = isTRUE(input$flt_data_deficient),
        order_filter         = paste(input$flt_order  %||% character(0), collapse = ";"),
        family_filter        = paste(input$flt_family %||% character(0), collapse = ";"),
        fixed_color_scale    = isTRUE(input$fix_color_scale),
        change_mode          = isTRUE(input$change_mode),
        baseline_year        = if (isTRUE(input$change_mode)) input$baseline_year %||% "" else ""
      )
      if (!is.null(s) && identical(s$kind, "hotspot")) {
        row$exposed_cells <- s$n_cells; row$exposed_species <- s$n_species
        row$top_order <- s$top_order;  row$top_family <- s$top_family
      } else if (!is.null(s) && identical(s$kind, "single")) {
        row$species <- s$species; row$exposed_cells <- s$n_cells; row$variables_exposed <- s$n_vars
      }
      utils::write.csv(as.data.frame(row, stringsAsFactors = FALSE), file, row.names = FALSE)
    }
  )

  output$download_map_raster_hs <- downloadHandler(
    filename = function() {
      if (input$mode == "single")
        paste0("raster_", gsub("[^A-Za-z0-9_]", "_", input$species %||% "species"), "_", input$year, ".tif")
      else
        paste0("raster_hotspot_", input$year, ".tif")
    },
    content = function(file) {
      r <- .build_map_raster_for_download()
      if (is.null(r)) { file.create(file); return() }
      suppressWarnings(terra::writeRaster(r, file, overwrite = TRUE))
    }
  )

  output$download_click <- downloadHandler(
    filename = function() paste0("clicked_cell_", input$year, ".csv"),
    content = function(file) {
      d <- click_data()
      if (is.null(d)) d <- data.frame()
      write.csv(d, file, row.names = FALSE)
    }
  )

  output$polygon_table <- renderDT({
    req(input$mode == "hotspot")
    cell_ids <- polygon_cells_rv()
    d <- tryCatch(polygon_table_data(), error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0) {
      msg <- polygon_empty_reason_rv()
      if (is.null(msg) || !nzchar(msg)) {
        if (is.null(cell_ids) || length(cell_ids) == 0) {
          msg <- "No exposed raster cells intersect the selected area."
        } else {
          msg <- "Cells intersect the selected area, but no species pass the current year/variable/filter settings."
        }
      }
      return(datatable(
        data.frame(Message = msg),
        rownames = FALSE, options = list(dom = "t", ordering = FALSE), class = "compact"
      ))
    }
    datatable(d, rownames = FALSE, filter = "top",
              options = list(pageLength = 30, dom = "ftp", scrollX = TRUE,
                             columnDefs = list(
                               list(width = "220px", targets = 1),
                               list(width = "72px",  targets = 2)
                             )),
              class = "compact stripe hover")
  })

  output$download_polygon <- downloadHandler(
    filename = function() paste0("polygon_species_", input$year, ".csv"),
    content = function(file) {
      d <- tryCatch(polygon_table_data(), error = function(e) NULL)
      if (is.null(d)) d <- data.frame()
      write.csv(d, file, row.names = FALSE)
    }
  )

  observeEvent(input$clear_polygon, {
    polygon_cells_rv(NULL)
    polygon_empty_reason_rv(NULL)
    leafletProxy("map") %>% clearGroup("drawn")
    set_status("Selection cleared", done = TRUE)
  })

  observeEvent(input$reset_all, {
    playing_rv(FALSE)
    updateSliderInput(session, "play_speed", value = 0.75)
    updateRadioButtons(session, "mode", selected = "hotspot")
    updateSelectInput(session, "year", selected = default_year)
    updateCheckboxInput(session, "change_mode", value = FALSE)
    updateSelectInput(session, "baseline_year", selected = min(years_avail))
    updateCheckboxGroupInput(session, "sel_vars", selected = vars_avail)
    updateCheckboxInput(session, "flt_exposed_only", value = FALSE)
    updateSelectizeInput(session, "species", selected = species_avail[1], server = TRUE)
    updateCheckboxInput(session, "flt_threatened", value = FALSE)
    updateCheckboxInput(session, "flt_data_deficient", value = FALSE)
    updateCheckboxGroupInput(session, "flt_groups",
                             selected = c("Amphibians", "Birds", "Mammals", "Reptiles"))
    updateSliderInput(session, "excl_threshold", value = 10)
    updateCheckboxInput(session, "fix_color_scale", value = FALSE)
    updateSelectizeInput(session, "flt_order", selected = "", server = TRUE)
    updateSelectizeInput(session, "flt_family", selected = "", server = TRUE)
    updateSelectizeInput(session, "pol_country", selected = "USA", server = TRUE)
    updateSelectizeInput(session, "pol_state", selected = "", server = TRUE)
    updateSelectizeInput(session, "pol_county", selected = "", server = TRUE)
    polygon_cells_rv(NULL)
    polygon_empty_reason_rv(NULL)
    click_data(NULL)
    hotspot_diag_data(NULL)
    click_info(list(cell = NA_integer_, lat = NA, lng = NA))
    map_raster_rv(NULL)
    exposure_msg_rv(NULL)
    hide_progress()
    leafletProxy("map") %>% clearGroup("drawn")
    set_status("All controls reset", done = TRUE)
  })

  observeEvent(input$shp_upload, {
    files <- input$shp_upload
    req(!is.null(files))
    shp_row <- files[grepl("\\.shp$", files$name, ignore.case = TRUE), ]
    if (nrow(shp_row) == 0) {
      showNotification("No .shp file found — please select all shapefile components.", type = "error")
      return()
    }
    # Copy all uploaded files into a temp dir preserving original extensions
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    for (i in seq_len(nrow(files))) {
      file.copy(files$datapath[i], file.path(tmp_dir, files$name[i]))
    }
    uploaded_paths <- file.path(tmp_dir, files$name)
    md5 <- tools::md5sum(uploaded_paths)
    shp_hash <- paste(paste(names(md5), unname(md5), sep = "="), collapse = "|")
    shp_path <- file.path(tmp_dir, shp_row$name[1])
    set_status("Reading uploaded shapefile...")
    show_progress(15, "Reading shapefile...")
    tryCatch({
      if (exists(shp_hash, envir = shp_cache, inherits = FALSE)) {
        show_progress(60, "Using cached shapefile intersection...")
        cached <- get(shp_hash, envir = shp_cache, inherits = FALSE)
        cell_ids <- cached$cell_ids
        gj <- cached$geojson
        ext <- cached$ext
      } else {
        shp_v <- terra::vect(shp_path)
        if (!identical(terra::crs(shp_v), terra::crs(tpl))) {
          show_progress(40, "Shapefile CRS differs — reprojecting to raster CRS...")
          set_status("Shapefile projection differs from raster; reprojecting now...")
          showNotification("Uploaded shapefile projection differs from raster. Reprojecting to match raster CRS...", type = "message", duration = 4)
          shp_v <- terra::project(shp_v, terra::crs(tpl))
        } else {
          show_progress(40, "Projection matches raster CRS...")
        }

        show_progress(60, "Finding exposed cells in shape...")
        if (!length(avail_cells_all)) {
          cell_ids <- integer(0)
        } else {
          inside <- tryCatch(
            terra::is.related(avail_pts, shp_v, "intersects"),
            error = function(e) rep(FALSE, length(avail_cells_all))
          )
          cell_ids <- avail_cells_all[inside]
        }
        gj <- spatvector_to_leaflet_geojson(shp_v)
        ext <- terra::ext(shp_v)
        assign(shp_hash, list(cell_ids = cell_ids, geojson = gj, ext = ext), envir = shp_cache)
      }

      polygon_cells_rv(if (length(cell_ids)) cell_ids else integer(0))
      if (length(cell_ids) == 0) {
        polygon_empty_reason_rv("No exposed raster cells intersect this uploaded shape.")
      } else {
        polygon_empty_reason_rv(NULL)
      }
      leafletProxy("map") %>%
        clearGroup("drawn") %>%
        addGeoJSON(gj, group = "drawn",
                   color = "#e67e00", weight = 2,
                   fillColor = "#e67e00", fillOpacity = 0.15)
      leafletProxy("map") %>% fitBounds(ext[1], ext[3], ext[2], ext[4])

      click_data(NULL); hotspot_diag_data(NULL)
      click_info(list(cell = NA_integer_, lat = NA, lng = NA))
      hide_progress()
      set_status(sprintf("Shapefile loaded: %d cells selected", length(cell_ids)), done = TRUE)
    }, error = function(e) {
      hide_progress()
      showNotification(paste("Shapefile read failed:", conditionMessage(e)), type = "error")
      set_status("Shapefile read failed")
    })
  })
}

shinyApp(ui, server)

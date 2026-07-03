#!/usr/bin/env Rscript
# Regenerate exposureApp data products applying the CORRECTED >=1% range-exposure
# filter. Mirrors the app-data block of src/r/10Tables.r, but builds the
# species x year keep-set from the UNROUNDED propExposed (>= 0.01) instead of the
# rounded value (> 0.01). See 10Tables.r comment at app_keep_sp_year.
#
# Inputs  (read-only):
#   output/v8/Int_V8/AllCellExposureSpXVar.qs      raw cell x year x var x sp + rangeSize
#   exposureApp/inst/extdata/allExpForShiny_species.qs   taxonomy for all 32,345 species
# Outputs (overwrite, git-tracked -> revert with `git checkout` if needed):
#   exposureApp/inst/extdata/allExpForShiny_by_year/allExpForShiny_<yr>.qs
#   exposureApp/inst/extdata/allExpForShiny.qs
#   exposureApp/inst/extdata/allExpForShiny_species.qs
#   exposureApp/inst/extdata/allExpForShiny_trends_cache_v1.qs
#   exposureApp/inst/extdata/allExpForShiny_species_year_cells_v1.qs
#   exposureApp/inst/extdata/allExpForShiny_manifest_v1.qs
# (hotspot cache rebuilt separately via build_hotspot_cache.R)

suppressPackageStartupMessages({ library(qs2); library(dplyr); library(stringr) })

PROJ  <- "/Users/cory.merow/Dropbox/Projects/2025_Exposure"
intDir <- file.path(PROJ, "output/v8/Int_V8")
ext    <- file.path(PROJ, "exposureApp/inst/extdata")
ydir   <- file.path(ext, "allExpForShiny_by_year")

message("Reading source AllCellExposureSpXVar.qs ...")
allExp <- qs2::qs_read(file.path(intDir, "AllCellExposureSpXVar.qs"))

# --- true (unrounded) range-exposed fraction, per species x year ---------------
pe <- allExp %>%
  group_by(spName, year) %>%
  summarise(propExposed = dplyr::n_distinct(cell) / dplyr::first(rangeSize),
            .groups = "drop")

# --- CORRECTED keep-set: unrounded fraction >= 1% ------------------------------
keep <- pe %>% filter(!is.na(year), propExposed >= 0.01) %>% select(spName, year)
message(sprintf("keep set: %d species x year combos (>=1%% of range)", nrow(keep)))

# taxonomy (already title-cased, covers all species)
spMeta <- qs2::qs_read(file.path(ext, "allExpForShiny_species.qs"))

# rounded display value, per species x year
pe_disp <- pe %>% mutate(propExposed = round(propExposed, 2))

message("Building app_rows ...")
app_rows <- allExp %>%
  select(cell, year, var, spName) %>%
  semi_join(keep,    by = c("spName", "year")) %>%
  left_join(pe_disp, by = c("spName", "year")) %>%
  left_join(spMeta,  by = "spName") %>%
  filter(!is.na(year), !is.na(var), !is.na(cell)) %>%
  mutate(
    spName          = as.character(spName),
    var             = as.character(var),
    group           = as.character(group),
    orderName       = str_to_title(orderName),
    familyName      = str_to_title(familyName),
    redlistCategory = as.character(redlistCategory),
    propExposed     = as.numeric(propExposed)
  ) %>%
  arrange(year, spName, var, cell)
message(sprintf("app_rows: %d rows, %d species", nrow(app_rows), dplyr::n_distinct(app_rows$spName)))

# --- derived products ----------------------------------------------------------
app_species_meta <- app_rows %>%
  distinct(spName, group, orderName, familyName, redlistCategory)

cell_trend_df <- app_rows %>% distinct(spName, cell, year) %>%
  count(cell, year, name = "n_sp")
sp_trend_df <- app_rows %>% group_by(spName, var, year) %>%
  summarize(mean_prop = mean(propExposed, na.rm = TRUE), .groups = "drop")
species_year_cells <- app_rows %>% distinct(year, spName, var, cell) %>%
  arrange(year, spName, var, cell) %>%
  group_by(year, spName, var) %>%
  summarize(cells = list(as.integer(cell)), .groups = "drop")

# --- write year shards ---------------------------------------------------------
if (!dir.exists(ydir)) dir.create(ydir, recursive = TRUE)
unlink(list.files(ydir, pattern = "^allExpForShiny_.*\\.qs$", full.names = TRUE))
app_years <- sort(unique(app_rows$year))
app_year_files <- setNames(character(length(app_years)), as.character(app_years))
for (yr in app_years) {
  yf <- file.path(ydir, paste0("allExpForShiny_", yr, ".qs"))
  qs2::qs_save(app_rows %>% filter(year == yr), file = yf)
  app_year_files[as.character(yr)] <- file.path("allExpForShiny_by_year", basename(yf))
}

# --- write monolith + caches + meta + manifest ---------------------------------
qs2::qs_save(app_species_meta, file.path(ext, "allExpForShiny_species.qs"))
qs2::qs_save(app_rows,         file.path(ext, "allExpForShiny.qs"))
qs2::qs_save(list(version = 1L, cell_trend_df = cell_trend_df, sp_trend_df = sp_trend_df),
             file.path(ext, "allExpForShiny_trends_cache_v1.qs"))
qs2::qs_save(species_year_cells, file.path(ext, "allExpForShiny_species_year_cells_v1.qs"))
qs2::qs_save(
  list(version = 1L, years_avail = app_years, year_files = as.list(app_year_files),
       species_file = "allExpForShiny_species.qs",
       trends_cache_file = "allExpForShiny_trends_cache_v1.qs",
       avail_cells_all = sort(unique(app_rows$cell))),
  file.path(ext, "allExpForShiny_manifest_v1.qs")
)
message("Done regenerating app data products.")

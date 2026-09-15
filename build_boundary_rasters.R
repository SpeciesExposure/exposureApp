# build_boundary_rasters.R --- rasterize GADM 4.1 ADM0/1/2 onto the app's 0.25-degree
# template grid so political-unit summaries need no downloads at run time.
#
# Output (inst/extdata):
#   gadm_ADM0.tif / gadm_ADM1.tif / gadm_ADM2.tif : integer unit id per cell (NA = none)
#   gadm_ADM0.csv / gadm_ADM1.csv / gadm_ADM2.csv : id -> codes, names, parent, n_cells
#
# Cell assignment: a cell belongs to the unit containing its centre. Cells whose
# centre is not inside any unit (coastal cells with an ocean centre, small islands)
# are then given to a unit that touches them, so that every app cell near land
# has an owner. Units that own no cell after the centre pass (small islands, thin
# districts) get first claim on touched cells; remaining cells go to any touching
# unit, ties resolved by rasterize().
# Units that still have no cells are left out of the tables (they can never be selected).
#
# Run once from the package root:  Rscript build_boundary_rasters.R
suppressPackageStartupMessages(library(terra))

gadm_path <- "/Users/cory.merow/Documents/SDMs/Exposure_2023/gadm_410-levels.gpkg"
out_dir   <- file.path("inst", "extdata")
tpl <- rast(file.path(out_dir, "landTemplate.tif"))

rasterize_level <- function(layer, key_cols) {
  t0 <- Sys.time()
  v <- vect(gadm_path, layer = layer)
  v$id <- seq_len(nrow(v))
  message(sprintf("%s: %d features read in %.0f s", layer, nrow(v), as.numeric(Sys.time() - t0, units = "secs")))

  r <- rasterize(v, tpl, field = "id", touches = FALSE)
  have <- unique(na.omit(values(r)[, 1]))
  missing <- setdiff(v$id, have)
  if (length(missing)) {
    r <- cover(r, rasterize(v[v$id %in% missing, ], tpl, field = "id", touches = TRUE))
  }
  r <- cover(r, rasterize(v, tpl, field = "id", touches = TRUE))   # cover() fills only cells still NA
  vals <- values(r)[, 1]
  n_cells <- tabulate(vals[!is.na(vals)], nbins = nrow(v))

  tab <- as.data.frame(v)[, c("id", key_cols)]
  tab$n_cells <- n_cells
  # GADM has a few units with blank names and one ADM1 row with no country code:
  # drop rows without codes, and label nameless units by their GADM id so they
  # remain selectable and distinguishable.
  code_cols <- grep("^GID_", key_cols, value = TRUE)
  tab <- tab[complete.cases(tab[, code_cols, drop = FALSE]), ]
  for (nm in grep("^NAME_|^COUNTRY$", key_cols, value = TRUE)) {
    gid <- if (nm == "COUNTRY") "GID_0" else sub("NAME_", "GID_", nm)
    bad <- is.na(tab[[nm]]) | !nzchar(trimws(tab[[nm]])) | tab[[nm]] == "NA"
    tab[[nm]][bad] <- paste0("(unnamed ", tab[[gid]][bad], ")")
  }
  message(sprintf("%s: %d units with cells, %d without, %d cells assigned, %.0f s total",
                  layer, sum(n_cells > 0), sum(n_cells == 0), sum(n_cells), as.numeric(Sys.time() - t0, units = "secs")))

  stopifnot(nrow(v) < 65535)
  writeRaster(r, file.path(out_dir, sprintf("gadm_%s.tif", sub("_", "", layer))),
              overwrite = TRUE, datatype = "INT2U", gdal = c("COMPRESS=DEFLATE"))
  write.csv(tab[tab$n_cells > 0, ], file.path(out_dir, sprintf("gadm_%s.csv", sub("_", "", layer))), row.names = FALSE)
  invisible(tab)
}

rasterize_level("ADM_0", c("GID_0", "COUNTRY"))
rasterize_level("ADM_1", c("GID_0", "COUNTRY", "GID_1", "NAME_1"))
rasterize_level("ADM_2", c("GID_0", "COUNTRY", "GID_1", "NAME_1", "GID_2", "NAME_2"))
message("done")

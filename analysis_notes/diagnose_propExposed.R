## diagnose_propExposed.R
## ---------------------------------------------------------------------------
## Two questions about the hotspot map's species counts:
##   (1) Does the map apply a ">1% of range exposed" threshold when counting?
##   (2) Why does `propExposed` run 0-4 instead of 0-1?
##
## Run from the repo root with the app data directory available.
## Uses only base R + qs2 + dplyr (the app's own deps).
## ---------------------------------------------------------------------------

suppressMessages({library(qs2); library(dplyr)})

intDir <- Sys.getenv("INT_DIR",
  unset = "exposureApp/inst/extdata")           # adjust if needed
byDir  <- file.path(intDir, "allExpForShiny_by_year")

## ---- (1) How the hotspot count is defined -------------------------------
## In app.R the per-cell species count is a PLAIN distinct-species count:
##   cell_trend_df:  distinct(spName, cell, year) |> count(cell, year)
##   build_hotspot_raster: group_by(cell) |> summarize(n = n_distinct(spName))
## There is NO propExposed / range-fraction gate anywhere in the count.
## A species is counted in a cell whenever a row exists for it (a selected
## climate variable crossed that species' historical threshold that year).

## ---- (2) What propExposed actually is ------------------------------------
d  <- qs_read(file.path(byDir, "allExpForShiny_2025.qs"))
rc <- qs_read(file.path(intDir, "range_cells_mb.qs"))   # per-species range cells

## propExposed is defined in src/r/10Tables.r line 45:
##   group_by(spName, year) |> mutate(propExposed = n() / rangeSize)
## allExp has one row per spName x var x cell x year, so n() sums exposed
## (cell x variable) events across ALL variables in that year. It is a
## SPECIES x YEAR scalar (identical across every var and cell of that species).

## (a) show it is constant across vars, and that it SUMS across vars.
sp <- "Lepilemur_edwardsi"; rsz <- length(rc[[sp]])
one <- d |> filter(spName == sp) |>
  group_by(var) |>
  summarize(propExposed = first(propExposed),
            exposed_cells_this_var = n_distinct(cell), .groups = "drop")
cat("rangeSize (range_cells_mb):", rsz, "cells\n")
print(as.data.frame(one))
cat(sprintf("2025: sum of exposed cells across vars = %d ; /rangeSize = %.3f = propExposed\n",
    sum(one$exposed_cells_this_var), sum(one$exposed_cells_this_var)/rsz))
## -> propExposed identical across vars; equals (rows that year, all vars)/rangeSize.
##    It exceeds 1.0 because a cell is counted once PER VARIABLE and 7 vars can
##    sum to more cell-events than the species has range cells. NOT cumulative
##    across years (there is no cumsum anywhere in the pipeline).

## (b) confirm it is a single-year quantity: 2024 vs 2025 across BOTH vars.
for (y in c(2024, 2025)) {
  s <- qs_read(file.path(byDir, sprintf("allExpForShiny_%d.qs", y))) |> filter(spName == sp)
  bv <- s |> group_by(var) |> summarize(cells = n_distinct(cell), .groups = "drop")
  cat(sprintf("%d: propExposed=%.2f  = (%s)/%d = %.2f\n", y, s$propExposed[1],
      paste(bv$cells, collapse = "+"), rsz, sum(bv$cells)/rsz))
}
## 2024: 0.95 = (15+4)/20 ; 2025: 1.05 = (13+8)/20 -- each computed from that
## year's cells only. (An earlier draft wrongly called this cumulative; it is not.)

## ---- (3) Inclusion floor: YES, a >1% range-fraction filter is applied ------
## src/r/10Tables.r lines 90-96 keep only species-years with propExposed > 0.01
## before writing the shards + cell_trend_df (semi_join on app_keep_sp_year).
## So a species appears on the hotspot map in a year only if >1% of its range
## is exposed that year. The floor is visible in the data:
pe <- d |> distinct(spName, propExposed) |> pull(propExposed)
pe <- pe[!is.na(pe)]
cat(sprintf("species-level propExposed: min=%.3f  #<=0.01=%d  #==0.02=%d  max=%.3f\n",
    min(pe), sum(pe <= 0.01), sum(pe == 0.02), max(pe)))
## Expect min=0.02, #<=0.01 = 0: that 0.02 is filter(propExposed > 0.01)
## combined with round(propExposed, 2) (line 70) -- the smallest 2-dp value
## strictly greater than 0.01 is 0.02. NOT mere rounding granularity.
## NOTE: the "no magnitude filter" comment (line 51) refers to not filtering by
## exposure MAGNITUDE, not to the range-fraction filter, which IS applied.

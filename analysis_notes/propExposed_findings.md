# `propExposed` and hotspot species counts — findings

Reproduce with [`diagnose_propExposed.R`](diagnose_propExposed.R)
(`INT_DIR=exposureApp/inst/extdata Rscript ...`).

## 1. Yes — a >1%-of-range filter is applied, but in the pipeline, not the app

There are two layers, and it is important to keep them separate:

**The app's counting code has no range-fraction gate.** The per-cell species
count is a plain distinct-species count:

- `cell_trend_df` (fast/precomputed path): `distinct(spName, cell, year) |> count(cell, year, name = "n_sp")`
- `build_hotspot_raster` (exact path): `group_by(cell) |> summarize(n = n_distinct(spName))`

`propExposed` is carried and displayed (mean per cell in the click table) but
is **never referenced by the app's counting logic**. So *within the app* there
is no threshold.

**But the data the app receives is already filtered at >1%.** The pipeline
script `src/r/10Tables.r` (lines 90–96) keeps only species-years above the
threshold before writing the shards and `cell_trend_df`:

```r
app_keep_sp_year <- bb %>%
  distinct(spName, year, propExposed) %>%
  filter(!is.na(year), propExposed > 0.01) %>%   # >1% of range exposed that year
  dplyr::select(spName, year)
app_filtered <- bb %>% semi_join(app_keep_sp_year, by = c('spName', 'year'))
```

Both output paths — the by-year shards (`allExpForShiny_<yr>.qs`) and
`cell_trend_df` (the fast-path hotspot count) — are built from `app_filtered`,
so **any species-year with `propExposed <= 0.01` is dropped before the app ever
sees it.** The net effect is exactly what the question asked: a species is
represented on the hotspot map in a given year only if its `propExposed` for
that year exceeds 1%.

Two caveats on interpreting that threshold:

1. It is applied to `propExposed` as defined in the same script. **Before** the
   fix in §4, `propExposed = n()/rangeSize` summed cell-events across variables
   (0–4), so ">1%" was against an inflated, not-quite-range-fraction quantity.
   **After** the fix (`n_distinct(cell)/rangeSize`), the same `> 0.01` filter
   becomes a clean "more than 1% of range cells exposed that year" gate.
2. The filter is at the **species-year** level, not per cell: once a species
   clears >1% for a year, *all* of its exposed cells that year are kept and each
   contributes to the per-cell count.

`range_cells_mb` (per-species range size) is used directly only in
single-species mode; the range-size denominator for `propExposed` comes from
the `rangeSize` column inside the pipeline.

## 2. Why `propExposed` runs 0–4 instead of 0–1 — it sums across variables

**Source located:** `src/r/10Tables.r`, line 45:

```r
allExp <- allExp %>%
  group_by(spName, year) %>%
  mutate(propExposed = n() / rangeSize) %>%
  ungroup()
```

This is the script that writes both `allExpForShiny.qs` and the
`allExpForShiny_by_year/*.qs` shards the app reads (same file, lines ~82–123).
So the app's `propExposed` column is defined exactly here.

`allExp` has **one row per `spName × var × cell × year`** (built in
`6_MakeDF.r` from `unnest(exposedYearCells)`; the per-year exposed-cell lists
come from `2_Functions.r`). Grouping by `spName, year` and taking `n()`
therefore counts **every exposed (cell × variable) event that year, summed
across all 7 climate variables**, then divides by the range size. A single cell
that trips 3 different variables contributes 3 to `n()`, and the 7 variables
can collectively count more exposed cell-events than the species has range
cells — so the ratio exceeds 1.0 (max observed 4.0).

**Correction to an earlier draft of this note:** this is **NOT** cumulative
across years — it is a single-year quantity that sums across variables. My
earlier "cumulative over years" claim was wrong; it was based on tracking one
variable (`temp__12_up`: 15→13 cells, 2024→2025) while missing that a second
variable (`temp__3__max_up`: 4→8) pushed the across-variable total up 19→21.
Verified directly: for *Lepilemur_edwardsi* (range = 20), 2025 `propExposed`
= 1.05 = 21/20, where 21 = 13 (`temp__12_up`) + 8 (`temp__3__max_up`) exposed
cells **in 2025 alone**; 2024 = 0.95 = 19/20 = (15 + 4)/20. No year-to-year
accumulation is involved, and there is no `cumsum` anywhere in the pipeline.

**What "proportion of range exposed that year" should be.** The intended
per-year fraction is *(number of distinct range cells exposed by any variable
in that year) / rangeSize*, which is bounded 0–1. The current formula differs
in two ways: (i) it counts a cell once **per variable** rather than once if
exposed by any variable, and (ii) `rangeSize` here is the `range_cells_mb`
count, which must be the correct denominator. To fix, change line 45 to count
**distinct exposed cells** per species-year instead of raw rows:

```r
allExp <- allExp %>%
  group_by(spName, year) %>%
  mutate(propExposed = dplyr::n_distinct(cell) / rangeSize) %>%
  ungroup()
```

(and re-run `10Tables.r` to regenerate the shards). `n_distinct(cell)` can
still equal `rangeSize` at most, so the result is bounded 0–1. Decide first
whether "exposed" should mean *exposed by any variable* (use `n_distinct(cell)`
as above) or *retain per-variable resolution* (then compute the fraction per
`spName, year, var` with `group_by(spName, year, var)` and it stays ≤1 per
variable).

## 3. There IS an inclusion floor at 1% — and it explains the 0.02 minimum

**Correction to an earlier draft of this note, which wrongly stated there was
no floor.** `10Tables.r` lines 90–96 apply `filter(propExposed > 0.01)` at the
species-year level before writing the shards (see §1). The comment near the
`propExposed` definition ("no magnitude filter…", line 51) refers only to *not*
filtering by exposure **magnitude** (the severity/`mag` of the crossing); it
does **not** mean there is no range-fraction filter — the `> 0.01`
range-fraction filter is applied a few lines later.

This floor is directly visible in the data: in the 2025 shard the minimum
distinct species-level `propExposed` is exactly **0.02**, there are **no values
≤ 0.01**, and 5002 species sit exactly at 0.02. That 0.02 minimum is the
footprint of `filter(propExposed > 0.01)` combined with `round(propExposed, 2)`
(line 70) — the smallest 2-dp value strictly greater than 0.01 is 0.02 — **not**
mere rounding granularity as an earlier draft claimed.

## 4. Fix applied and its blast radius

**Fix (applied):** `src/r/10Tables.r` line 49 now uses
`propExposed = dplyr::n_distinct(cell) / rangeSize` grouped by `spName, year`,
so each range cell is counted at most once per year regardless of how many
variables expose it. Verified over all shards: the new definition has
**max = 1.00, zero species-years > 1** (the old `n()/rangeSize` gave max 4.0
with 1997 species-years over 1). **The shards must be regenerated by re-running
`10Tables.r`** for the app to pick this up.

**Why this error appeared only in the app table, not in the figures.** The bug
is specific to grouping by `spName, year` and counting **raw rows** (`n()`),
which sums exposed cell-events across all variables. Every figure/manuscript
script instead computes exposure **per variable** — grouping includes `var`, or
the fraction is built from a single variable's `X<year>` column — so a cell is
never double-counted across variables and the ratio stays ≤1. Confirmed
empirically: the per-variable grouping gives max = 1.00, zero over.

Other places that touch a `propExposed = n()/rangeSize` construction, and
whether they carry the same over-count:

- **Carries the identical over-count:** `docs/paper_code_data_summary.Rmd`
  line ~2251 — a code/data-summary document that reproduces the app-table
  build with the *same* `group_by(spName, year) + n()/rangeSize`. Its
  `propExposed` (and the reported `propExposed > 5%` / `> 1%` row and species
  counts) inherit the over-count. If those numbers appear in the
  paper/supplement, apply the same `n_distinct(cell)` fix and re-knit.
- **Does NOT carry the same inflation:** `src/r/6MapJointExposure.r` line 162
  uses `group_by(spName, year, var)` with `n()/first(rangeSize)` — because the
  grouping includes `var` it never sums across variables, so it is not the
  0–4 inflation. (Within a single variable, `n()` still counts rows rather than
  distinct cells, so it could exceed the range only if a cell recurs within one
  variable-year — worth a glance, but a different and much smaller concern.)

Scripts that were **already correct** (per-variable, bounded): `12WhereInRange.r`
(`group_by(spName, var)`, `n_distinct(cell)`), and all `MS2023` /
`MS_2025_Bioscience` figure scripts (per-variable `X<year>/rangeSize`).

**Not automatically covered by the pipeline fix:** any exposure numbers already
baked into a saved figure, table, or the manuscript text were produced by their
own (per-variable) scripts and are unaffected by this bug — but the app's
displayed "Proportion exposed" column and the summary card were affected until
the shards are rebuilt.

# AI handoff: species exposure to climate extremes, 1941–2025

## Read this first: what this file is

This file is an **AI handoff**. It is a briefing written so that an AI assistant
(such as Claude, ChatGPT or Gemini), or a person, can start working with this
project's data without any prior knowledge of it.

**How to use it.** Give this file to an AI assistant, by attaching it or pasting
its text, and then ask your question in plain language. For example: "Using the
files described here, write R code to list the species exposed in Brazil in
2024." The assistant will then know what each data file contains, what the
numbers mean, and which mistakes to avoid. You can also simply read it yourself;
it is ordinary text.

**What it is not.** It is not the paper and it does not replace the paper's
Methods. Where this file and the published paper disagree, the paper is the
authority for what was reported, and the code is the authority for what was
computed.

**How far to trust it.** Every file description and every code example below
was run against the shipped data on 2026-09-29 (app version 0.13.3). Statements
about how the data were generated were checked against the analysis scripts.
Items that could not be verified are marked as such in section 9.

---

## 1. The project in brief

The project asks, for every year from 1941 to 2025, which terrestrial vertebrate
species experienced climate conditions more extreme than anything recorded in
their range during a historical baseline period, and where.

- **Species:** 32,345 terrestrial vertebrates: 10,316 birds, 9,202 reptiles,
  7,498 amphibians, 5,329 mammals.
- **Climate:** ERA5 monthly reanalysis, 2 m air temperature and total
  precipitation.
- **Grid:** 0.25 degrees (about 28 km at the equator), land cells only.
- **Result:** for each species, year and climate variable, the set of grid cells
  in the species' range where that year's value went beyond the species'
  historical limit.

Two products come out of it:

1. A paper, "The 2025 Biodiversity Exposure Report" (Merow et al., in revision at
   *BioScience* as of September 2026; the title may change).
2. An R Shiny app, the `exposureApp` R package, which ships the derived data
   described here. Repository: <https://github.com/SpeciesExposure/exposureApp>.
   Project site: <https://speciesexposure.github.io/#home>.

---

## 2. Definitions you must get right

### 2.1 Year

A "year" is a **water year**: year Y runs from October of Y−1 to September of Y.
So "2025" means October 2024 to September 2025. Partial years 1940 and 2026 are
not included.

### 2.2 Climate variables

Eight variables. Four climate measures, each taken at the tail that represents
an extreme:

| Code | Meaning | Extreme |
|---|---|---|
| `temp__12_up` | Mean temperature over the 12 months of the water year | unusually hot year |
| `temp__12_lo` | same | unusually cold year |
| `temp__3__max_up` | Warmest 3-month period of the year (mean of 3 consecutive months) | unusually hot warm season |
| `temp__3__min_lo` | Coldest 3-month period of the year | unusually cold cold season |
| `precip__12_up` | Total precipitation over the water year | unusually wet year |
| `precip__12_lo` | same | unusually dry year |
| `precip__3__max_up` | Wettest 3-month period of the year | unusually wet wet season |
| `precip__3__min_lo` | Driest 3-month period of the year | unusually dry dry season |

The 3-month periods are consecutive months inside one water year. They do not
span September to October, so there are 10 candidate periods per year.

### 2.3 Exposure

For one species and one variable:

1. **In each grid cell of the species' range**, take the 99th percentile of the
   variable over the baseline years 1941–2022 (the 1st percentile for the
   `_lo` variables).
2. **Across the cells of the range**, take the 99th percentile of those cell
   values (the 1st percentile for `_lo`). This single number is the species'
   historical limit for that variable.
3. A cell is **exposed** in a given year if its value that year is strictly
   greater than the limit (`_up`) or strictly less than it (`_lo`).

Percentiles use R's default method (`quantile(..., type = 7)`). The baseline is
fixed at 1941–2022. It does not move with the year being evaluated.

Exposure measures a hazard: conditions outside what the species is known to
have experienced. It is not a measurement of mortality, population decline or
extinction risk.

### 2.4 Fraction of range exposed: there are two, and they differ

This is the most common source of mismatched numbers.

| | Single-variable fraction | Pooled fraction |
|---|---|---|
| What it is | Cells exposed to one variable, divided by range cells | Cells exposed to any of the 8 variables, each cell counted once, divided by range cells |
| Used by | **The paper.** A species counts as exposed in a year if this fraction is at least 0.25 for at least one variable | **The app.** This is the `propExposed` column and the app's threshold slider |
| Size | | Always at least as large as the single-variable fraction |
| In the shipped data | Not stored. Compute it as shown in section 4.4 | Stored, rounded to 2 decimals |

For 2025 the paper's rule gives **2,662** species. The pooled fraction at the
same 0.25 cut-off gives 2,789. Both are correct answers to different questions.

### 2.5 The 1% floor

The shipped data contain only species-years in which the pooled fraction is at
least 0.01. Species-years below that are absent. This was done to keep the
files small. It means:

- The smallest value of `propExposed` is 0.01.
- You cannot compute counts for thresholds below 1% from these files.
- All 32,345 species appear at least once across the 85 years.

### 2.6 Range size

`rangeSize` is the number of 0.25-degree land cells in the species' range. It
is the denominator of both fractions. It is in the species table from app
version 0.13.3 onward.

---

## 3. What ships with the app

All files are in `inst/extdata/` of the package. `.qs` files are read with
`qs2::qs_read()`. Rasters are read with `terra::rast()`.

### 3.1 The grid

`landTemplate.tif`: 588 rows by 1,440 columns, 846,720 cells, 0.25 degrees,
longitude −180 to 180, latitude −61 to 86, EPSG:4326. Land cells have the value
0 and everything else is NA. There are 244,432 land cells.

**Cell ids** everywhere in the data are the row-major cell numbers of this
raster, starting at 1 in the north-west corner, as used by `terra`. Convert
with `terra::xyFromCell(tpl, cell)` and `terra::cellFromXY(tpl, xy)`.

### 3.2 Exposure data

| File | Size | Content |
|---|---|---|
| `allExpForShiny_by_year/allExpForShiny_<year>.qs` | 98 MB total, 85 files | One row per species, variable and exposed cell for that year. Columns: `cell`, `year`, `var`, `spName`, `propExposed`, `group`, `orderName`, `familyName`, `redlistCategory`. The 2025 file has 1,041,735 rows and 17,861 species |
| `allExpForShiny_species.qs` | 0.4 MB | One row per species, 32,345 rows. Columns: `spName`, `group`, `orderName`, `familyName`, `redlistCategory`, `rangeSize` |
| `allExpForShiny_species_year_cells_v1.qs` | 80 MB | One row per species, year and variable, 2,707,589 rows. Columns: `year`, `spName`, `var`, and `cells`, a list column holding the integer cell ids exposed. Same information as the yearly files, in compact form, for all years at once |
| `allExpForShiny_trends_cache_v1.qs` | 8 MB | A list with two tables. `cell_trend_df`: `cell`, `year`, `n_sp`, the number of species exposed in that cell that year. `sp_trend_df`: `spName`, `var`, `year`, `mean_prop` |
| `allExpForShiny_manifest_v1.qs` | 56 KB | A list: `years_avail`, `year_files`, `species_file`, `trends_cache_file`, and `avail_cells_all`, the 215,104 cells that appear anywhere in the data |
| `allExpForShiny_hotspot_cache_v1.qs` | 0.6 MB | Precomputed species counts per cell for 2021–2025, used only to draw the app's first screen quickly. Keys `"2025"` (1% floor) and `"2025@10"` (pooled fraction at least 0.10) |
| `range_cells_mb.qs` | 9 MB | For 15,645 birds and mammals, the cells where the species was exposed in any year 1941–2025. See the warning in section 8 |
| `allExpForShiny.qs` | 97 MB | The older single-file form: all yearly files stacked, 57,074,639 rows. The app does not read it when the manifest is present |

Notes on columns:

- `spName` is `Genus_species` with an underscore.
- `group` is one of `Amphibians`, `Birds`, `Mammals`, `Reptiles`.
- `redlistCategory` is the IUCN Red List category with the first space replaced
  by an underscore, for example `Least_Concern`, `Critically_Endangered`,
  `Extinct_in the Wild`. It is NA for 895 species with no assessment.
- `propExposed` is the pooled fraction of section 2.4, one value per species
  and year, repeated on every row of that species-year. It does not change
  with `var`.
- `mean_prop` in `sp_trend_df` is that same pooled value. Although the table
  has a `var` column, `mean_prop` is identical across variables within a
  species-year. It is not a per-variable fraction.

### 3.3 Regions

All on the same grid as the template.

| File | Content |
|---|---|
| `ecoregions2017_ECO_ID.tif`, `ecoregions2017_meta.csv` | Ecoregion id per cell, and a table of 795 ecoregions with `ECO_ID`, `ECO_NAME`, `BIOME_NAME`, `REALM` (RESOLVE Ecoregions 2017). This is the same cell-to-ecoregion assignment used for the paper's ecoregion results |
| `gadm_ADM0.tif`, `gadm_ADM0.csv` | Country id per cell; 256 countries. Columns `id`, `GID_0` (ISO3 code), `COUNTRY`, `n_cells` |
| `gadm_ADM1.tif`, `gadm_ADM1.csv` | State or province; 3,272 units. Adds `GID_1`, `NAME_1` |
| `gadm_ADM2.tif`, `gadm_ADM2.csv` | County or district; 29,325 units. Adds `GID_2`, `NAME_2` |

Boundaries are GADM 4.1. A cell belongs to the unit that contains its centre.
Coastal cells whose centre is in the sea are given to a unit that touches
them. Units too small to own any cell are not listed: about 38% of GADM
counties and 11% of states. 98% of data cells have a county and effectively
all have a country and a state.

---

## 4. Working with the data

All examples are R and were run against the shipped files.

### 4.1 Setup

```r
library(qs2); library(dplyr); library(terra)
ext  <- system.file("extdata", package = "exposureApp")   # or the path to inst/extdata
tpl  <- rast(file.path(ext, "landTemplate.tif"))
meta <- qs_read(file.path(ext, "allExpForShiny_species.qs"))
d25  <- qs_read(file.path(ext, "allExpForShiny_by_year", "allExpForShiny_2025.qs"))
```

### 4.2 Cell id to coordinates and back

```r
xyFromCell(tpl, 493319)                  # lon 29.625, lat 0.375
cellFromXY(tpl, cbind(29.625, 0.375))    # 493319
```

### 4.3 The app's hotspot map: species per cell

Species exposed per cell in 2025, keeping species whose pooled fraction is at
least 0.10, which is the app's default.

```r
hs <- d25 %>% filter(propExposed >= 0.10) %>%
  group_by(cell) %>% summarize(n_sp = n_distinct(spName), .groups = "drop")
nrow(hs)        # 5,821 cells
max(hs$n_sp)    # 111 species in the richest cell

r <- rast(tpl); values(r) <- NA_integer_; r[hs$cell] <- hs$n_sp
names(r) <- "n_species_2025_thr10"
writeRaster(r, "hotspot_2025_thr10.tif", overwrite = TRUE)
```

4,652 species contribute to this map.

### 4.4 The paper's count: species exposed in a year

At least 25% of the range exposed to at least one single variable.

```r
exposed_2025 <- d25 %>%
  count(spName, var, name = "exposed_cells") %>%
  inner_join(meta, by = "spName") %>%
  mutate(fraction = exposed_cells / rangeSize) %>%
  group_by(spName) %>% slice_max(fraction, n = 1, with_ties = FALSE) %>% ungroup() %>%
  filter(fraction >= 0.25)
nrow(exposed_2025)            # 2,662
table(exposed_2025$group)     # Amphibians 1060, Birds 255, Mammals 315, Reptiles 1032
```

The `var` column of the result is the variable with the largest fraction, which
the paper calls the dominant variable.

### 4.5 The same count for every year

```r
syc <- qs_read(file.path(ext, "allExpForShiny_species_year_cells_v1.qs"))
per_year <- syc %>%
  mutate(n = lengths(cells)) %>%
  inner_join(meta %>% select(spName, rangeSize), by = "spName") %>%
  filter(n / rangeSize >= 0.25) %>%
  group_by(year) %>% summarize(n_species = n_distinct(spName), .groups = "drop")
```

This gives 499 for 2021, 2,266 for 2023, 6,077 for 2024 and 2,662 for 2025.

### 4.6 One species through time

```r
sp <- "Panthera_leo"
syc %>% filter(spName == sp) %>%
  group_by(year) %>%
  summarize(cells = length(unique(unlist(cells))), .groups = "drop") %>%
  mutate(pooled_fraction = cells / meta$rangeSize[meta$spName == sp])
```

For the lion in 2025 this gives 90 exposed cells of 3,304, a pooled fraction of
0.027.

### 4.7 Species exposed in a country or an ecoregion

```r
# country
g0 <- rast(file.path(ext, "gadm_ADM0.tif"))
t0 <- read.csv(file.path(ext, "gadm_ADM0.csv"), na.strings = character(0))
cells <- which(values(g0)[, 1] == t0$id[t0$GID_0 == "MDG"])           # 933 cells
d25 %>% filter(cell %in% cells, propExposed >= 0.10) %>% distinct(spName) %>% nrow()   # 557

# ecoregion
e  <- rast(file.path(ext, "ecoregions2017_ECO_ID.tif"))
em <- read.csv(file.path(ext, "ecoregions2017_meta.csv"))
ec <- which(values(e)[, 1] == em$ECO_ID[em$ECO_NAME == "Madagascar humid forests"])    # 151 cells
d25 %>% filter(cell %in% ec, propExposed >= 0.10) %>% distinct(spName) %>% nrow()      # 243
```

Read `gadm_*.csv` with `na.strings = character(0)`. Otherwise R reads the
country code `NA` for Namibia, and any unit literally named "NA", as missing.

A species is listed for a region if any of its exposed cells fall in the
region. The fraction attached to it is still the fraction of its whole range,
not of the part inside the region.

---

## 5. Reproducing the paper's numbers from the shipped data

| Quantity | Paper | From shipped data | Note |
|---|---|---|---|
| Species exposed in 2025 | 2,662 | 2,662 | Section 4.4. Identical species set |
| By group, 2025 | 1,060 / 1,032 / 315 / 255 | same | Amphibians, reptiles, mammals, birds |
| Species exposed in 2024 | 6,078 | 6,077 | One species, *Thamnophilus sticturus*, is at 0.2496 in the shipped data and 0.25 in the analysis table |
| Species exposed in 2021 | 499 | 499 | |
| Every year 1990–2025 | | matches in 35 of 36 years | 2024 is the only difference |

The paper reports counts from 1990 onward. The shipped data go back to 1941.

Maps in the paper that show species per cell count every species with at least
one exposed cell there, without the 25% rule, so they show more species than
the headline counts.

---

## 6. How the data were generated

The analysis code is in a separate repository (`2025_Exposure`, scripts under
`src/r/`). The steps, in order:

| Step | Script | What it does |
|---|---|---|
| 1 | `1.8a_ClimatePrep.r` | Reads ERA5 monthly temperature and precipitation, puts them on the 0.25-degree grid, builds the land template from the ERA5 land–sea mask (land fraction above 0.5), and computes the water-year variables of section 2.2. Output: one table per variable with a row per cell and a column per year, 1941–2025 |
| 2 | `1.7MakeAoHRanges.r` | Records how species ranges were put on the grid. Mammals and birds: Area of Habitat rasters (Lumbierres et al. 2022), aggregated by a factor of 16 and assigned to land cells. Migratory birds use breeding and resident ranges only. Amphibians and reptiles: IUCN Red List range polygons, rasterized so that any cell a polygon touches is included. Output: one file per species holding its cell ids |
| 3 | `2_Metadata.r` | Builds the species attributes table: IUCN category, taxonomy, range size, ecoregions |
| 4 | `2CalcExposure.r` with functions in `1.5_Functions.r` | For every species and variable, computes the historical limit of section 2.3 and the cells exceeding it in each year 1941–2025. This is the longest step and takes many hours on a multi-core machine |
| 5 | `3_MakeDF.r` | Collects results. Writes `AllCellExposureSpXVar.qs`, every exposed cell-year for the 8 variables, and `singleMaxExposure_v1.qs`, one row per species and year where the largest single-variable fraction is at least 0.25, for 1990–2025. The paper's counts come from the second file |
| 6 | `10Tables.r` | Builds the app data: computes the pooled fraction, applies the 1% floor on the unrounded value, rounds `propExposed` to 2 decimals, and writes the yearly files, species table, trends cache, cell-list table and manifest |
| 7 | `prep_range_cells_mb.R` | Writes `range_cells_mb.qs` |
| 8 | In the app repository: `build_hotspot_cache.R`, `build_boundary_rasters.R`, `add_range_size.R` | Precomputed landing view, region rasters from GADM 4.1, and a repair script that adds `rangeSize` to the species table if it is missing |

Sensitivity runs with other baseline end years (1997, 2009, 2015, 2024) and
with a calendar year in place of the water year exist as separate run variants.
The shipped data are the main run: baseline 1941–2022, water year.

Input data sources: ERA5 monthly means from the Copernicus Climate Data Store;
IUCN Red List spatial data; Area of Habitat maps from Lumbierres et al. (2022);
RESOLVE Ecoregions 2017; GADM 4.1.

---

## 7. The app

Run it from an installed package with `exposureApp::exposureApp()`, or from a
checkout:

```r
Sys.setenv(INT_DIR = "<path>/inst/extdata"); shiny::runApp("<path>/inst/app")
```

- **Hotspot mode** maps the number of exposed species per cell for one year.
- **Single-species mode** maps one species' exposed cells, coloured by variable.
- **Filters:** climate variables, taxonomic group, order, family, threatened
  only, Data Deficient only, and the exposure threshold, which is the pooled
  fraction of section 2.4, default 10%, minimum 1%.
- **Change against a baseline year** maps the difference in species counts
  between two years. It uses all species and variables at the 1% floor and
  ignores the filters.
- **Region summaries** by country, state, county, ecoregion, a polygon drawn on
  the map, or an uploaded shapefile.
- **Downloads:** map as GeoTIFF, summary statistics, species tables as CSV.
- **Shareable links:** the address bar records the view. Parameters: `mode`,
  `year`, `vars`, `excl` (threshold in percent), `groups`, `order`, `family`,
  `threat`, `dd`, `change`, `baseline`, `fix`, `species`, `country`, `eco`.

The app's default view, pooled fraction at least 10%, shows 4,652 species for
2025. That is a different rule from the paper's 2,662. See section 2.4.

---

## 8. Known caveats and pitfalls

1. **Two fractions.** `propExposed` is pooled across variables. The paper's
   rule uses the largest single-variable fraction. Section 2.4.
2. **`propExposed` is rounded to 2 decimals.** Do not back out range sizes or
   apply fine thresholds from it. Use cell counts and `rangeSize`.
3. **1% floor.** Species-years with a pooled fraction below 0.01 are absent.
4. **Baseline years are in the data.** Years 1941–2022 were used to set the
   limits, so exposure in those years is exceedance of a 99th percentile
   within the sample that defined it. Expect a low, non-zero level in every
   baseline year by construction. Years after 2022 are outside the baseline.
5. **`range_cells_mb.qs` is not the species' range.** It holds the cells where
   the species was exposed in at least one year. The app draws it as a grey
   underlay labelled as such. The true range cells are not shipped; only their
   count, `rangeSize`, is.
6. **Large ranges dominate maps at low thresholds.** A species with a very
   large range that barely passes 1% still puts hundreds of cells on the map.
   Cell counts on the map track range size far more than they track the
   fraction exposed.
7. **Very large ranges.** 70 species have more than 50,000 range cells and 8
   have more than 100,000, up to 214,679, which is 88% of all land cells. All
   are birds or mammals, several of them migratory birds. Their fractions are
   correspondingly small. The cause has not been investigated.
8. **Regions are assigned by cell.** Region outlines in the app follow cell
   edges, and units smaller than a cell are not available.
9. **Static ranges.** Ranges are fixed for the whole period.
10. **Precipitation units.** The variables are sums or means of ERA5 monthly
    values. The physical unit was not established from the code. Exposure is
    relative to each species' own limit, so it does not depend on the unit.
11. **IUCN category strings** contain an underscore in place of the first
    space only. Compare after replacing underscores with spaces.
12. **Ecoregions.** 1.7% of data cells have no ecoregion in the lookup and are
    in no ecoregion summary.

---

## 9. What was and was not verified for this document

**Verified by running code on the shipped files:** every file's size, row
count and columns; the grid; every example and number in sections 4 and 5;
the species, group and IUCN totals; the region tables.

**Verified by reading the analysis scripts:** the water-year rule; the
variable definitions; the two-stage 99th percentile limits and percentile
method; the strict exceedance test; the 25% rule; the 1% floor; the land mask;
how ranges were rasterized.

**Not verified:** the ERA5 precipitation unit; the native resolution of the
Area of Habitat rasters before aggregation; which of two equivalent scripts
produced the shipped yearly files; the cause of the very large ranges in
caveat 7; the one-species difference for 2024.

---

## 10. Glossary

| Term | Meaning |
|---|---|
| Cell | One 0.25-degree grid square, identified by its number in `landTemplate.tif` |
| Baseline | 1941–2022, the years used to set each species' limits |
| Historical limit | The species' threshold for one variable, from section 2.3 |
| Exposed cell | A range cell whose value in a year is beyond the limit |
| Single-variable fraction | Exposed cells for one variable divided by range cells |
| Pooled fraction, `propExposed` | Cells exposed to any variable divided by range cells |
| Exposed species (paper) | Single-variable fraction at least 0.25 for some variable. Earlier reports in this series called this "highly exposed" |
| Dominant variable | The variable with the largest single-variable fraction |
| Hotspot | A cell with many exposed species |
| Water year | October of the previous year to September of the named year |

---

## 11. Citation and contact

Merow C, Wilson A, Maitner BS, Pillet M, Nikolopoulos E, Serra-Diaz JM, Meyer
AS, Trisos C, Pigot A, Urban MC. The 2025 Biodiversity Exposure Report. In
revision, *BioScience*. Check the project site for the final reference.

App and data: <https://github.com/SpeciesExposure/exposureApp>.
Corresponding author: Cory Merow.

Document written 2026-09-29 for app version 0.13.3.

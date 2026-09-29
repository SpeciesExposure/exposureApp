# exposureApp

`exposureApp` is an R package that launches an interactive Shiny app for exploring the exposure of terrestrial vertebrate species to climate extremes, for every year from 1941 to 2025. The derived exposure data ship with the package, so the app runs locally with no downloads.

Project site and web dashboard: <https://speciesexposure.github.io/#home>

## What the app shows

- **Hotspot maps.** The number of exposed species in each 0.25-degree land cell for a chosen year.
- **Single-species maps.** The exposed cells of one species, coloured by climate variable, with its history through time.
- **Change against a baseline year.** The difference in species counts between two years.
- **Region summaries.** By country, state or province, county, ecoregion, a polygon drawn on the map, or an uploaded shapefile.
- **Filters.** Climate variable, taxonomic group, order, family, IUCN threat status, and the minimum fraction of a species' range that must be exposed.
- **Downloads.** Maps as GeoTIFF and species tables as CSV.

## What "exposed" means

A species is exposed in a grid cell in a given year when that year's temperature or precipitation in the cell is more extreme than the limit of what the species experienced across its range during the historical baseline, 1941 to 2022. Exposure is a measure of hazard. It does not by itself show that a population was harmed.

Eight climate variables are used, covering annual and seasonal temperature and precipitation. A "year" is a water year: year Y runs from October of Y-1 to September of Y.

## Requirements

- R 4.2 or later
- System libraries needed by spatial R packages such as `terra`

The app uses these R packages:

```r
install.packages(c(
  "shiny", "qs2", "dplyr", "terra", "leaflet",
  "leaflet.extras", "jsonlite", "DT", "raster"
))
```

## Install

Install from the `main` branch. The download is about 300 MB because it includes the data.

```r
install.packages(
  "https://github.com/SpeciesExposure/exposureApp/archive/refs/heads/main.tar.gz",
  repos = NULL,
  type = "source"
)
```

## Launch

```r
exposureApp::exposureApp()
```

To run from a checkout of this repository without installing:

```r
Sys.setenv(INT_DIR = "<path to checkout>/inst/extdata")
shiny::runApp("<path to checkout>/inst/app")
```

Launch options:

```r
exposureApp::exposureApp(
  intDir = "/path/to/data",   # defaults to the data shipped with the package
  launch.browser = TRUE,
  host = "127.0.0.1",
  port = 3838
)
```

## Working with the data directly

The data files are in `inst/extdata/`. They can be read in R without running the app.

[AI_HANDOFF.md](AI_HANDOFF.md) documents every shipped file, the definitions behind the numbers, tested R examples, how the data were generated, and known caveats. It is written so that it can be given to an AI assistant, and it can equally be read by a person.

A short example, the number of exposed species per cell in 2025:

```r
library(qs2); library(dplyr); library(terra)

ext  <- system.file("extdata", package = "exposureApp")
rows <- qs_read(file.path(ext, "allExpForShiny_by_year", "allExpForShiny_2025.qs"))
tpl  <- rast(file.path(ext, "landTemplate.tif"))

n <- rows %>% filter(propExposed >= 0.10) %>% distinct(cell, spName) %>% count(cell)
r <- tpl; values(r) <- NA; r[n$cell] <- n$n
plot(r)
```

## Data sources

- Climate: ERA5 monthly reanalysis, Copernicus Climate Data Store
- Species ranges: IUCN Red List spatial data, BirdLife International, and Area of Habitat maps from Lumbierres et al. (2022)
- Ecoregions: RESOLVE Ecoregions 2017
- Political boundaries: GADM 4.1

## License

MIT. See [LICENSE](LICENSE).

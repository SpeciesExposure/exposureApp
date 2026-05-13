# exposureApp

`exposureApp` is an R package that launches an interactive Shiny app for exploring species climate exposure, including hotspot and single-species views.

## Requirements

- R (>= 4.2)
- System libraries needed by spatial R packages (for example `terra`)

The app runtime uses these R packages:

- `shiny`
- `qs2`
- `dplyr`
- `terra`
- `leaflet`
- `leaflet.extras`
- `jsonlite`
- `DT`
- `raster`

## Install

Install from the `main` branch tarball:

```r
install.packages(
	"https://github.com/SpeciesExposure/exposureApp/archive/refs/heads/main.tar.gz",
	repos = NULL,
	type = "source"
)
```

If needed, install missing runtime packages:

```r
install.packages(c(
	"shiny", "qs2", "dplyr", "terra", "leaflet",
	"leaflet.extras", "jsonlite", "DT", "raster"
))
```

## Launch

Start the app with packaged `extdata`:

```r
exposureApp::exposureApp()
```

Use custom data directory and launch options:

```r
exposureApp::exposureApp(
	intDir = "/path/to/int_data",
	launch.browser = TRUE,
	host = "127.0.0.1",
	port = 3838
)
```

`intDir` should contain required files such as `allExpForShiny.qs`, `landTemplate.tif`, and species metadata used by the app.

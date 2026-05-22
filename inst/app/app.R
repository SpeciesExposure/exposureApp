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
stopifnot(file.exists(allcells_path), file.exists(template_path))

log_msg("Loading exposure data: ", basename(allcells_path), " ...")
all_df <- read_qs(allcells_path)
log_msg("Loading template raster ...")
tpl    <- rast(template_path)

index_path <- file.path(intDir, "allExpForShiny_year_index_v1.qs")

# ── 2. Slim to only columns the app needs (keep all_df untouched) ────────────
# Keep group and prejoined species attributes when available
app_df <- all_df[, intersect(c("spName", "cell", "var", "year", "group",
                               "orderName", "familyName", "redlistCategory", "propExposed"), names(all_df))]

# Convert high-cardinality strings to factors: less memory, faster grouping
app_df$spName <- as.factor(app_df$spName)
app_df$var    <- as.factor(app_df$var)
log_msg(sprintf("app_df: %d rows, %.0f MB (full all_df was %.0f MB)",
                nrow(app_df),
                as.numeric(object.size(app_df)) / 1e6,
                as.numeric(object.size(all_df)) / 1e6))

to_title_case <- function(x) {
  x <- as.character(x)
  x <- trimws(tolower(x))
  x <- gsub("(^|[[:space:]-])([[:alpha:]])", "\\1\\U\\2", x, perl = TRUE)
  x[nchar(x) == 0] <- NA_character_
  x
}

# ── 3. Use species attributes already present in allExpForShiny ───────────────
if (!("orderName" %in% names(app_df))) app_df$orderName <- NA_character_
if (!("familyName" %in% names(app_df))) app_df$familyName <- NA_character_
if (!("redlistCategory" %in% names(app_df))) app_df$redlistCategory <- NA_character_
if (!("propExposed" %in% names(app_df))) app_df$propExposed <- NA_real_

# Enforce display style: apply regex only on unique values (fast for millions of rows)
.uniq_fmt <- function(x) {
  u <- unique(as.character(x))
  m <- setNames(to_title_case(u), u)
  m[as.character(x)]
}
app_df$orderName  <- .uniq_fmt(app_df$orderName)
app_df$familyName <- .uniq_fmt(app_df$familyName)

# left_join can coerce factors; restore factor columns used by controls/filters
app_df$spName <- as.factor(app_df$spName)
app_df$var    <- as.factor(app_df$var)

# ── 4. Pre-index by year for fast filtering ───────────────────────────────────
# Use row-index split so app_df is stored once (faster startup, lower memory)
# Persist this index so most launches skip rebuild work.
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

# ── 4b. Precompute trend lookup tables — cached to disk for fast startup ─────
trends_cache_path <- file.path(intDir, "allExpForShiny_trends_cache_v1.qs")
cell_trend_df <- NULL
sp_trend_df   <- NULL
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
}
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

rm(all_df); gc()   # free the original (large) object now that we are done with it

log_msg(sprintf("Ready. %d rows across %d years, %d species.",
                nrow(app_df), length(year_row_idx), length(levels(app_df$spName))))

years_avail    <- sort(as.integer(names(year_row_idx)))
species_avail  <- sort(unique(as.character(app_df$spName)))
vars_raw       <- unique(as.character(app_df$var))
vars_avail     <- vars_raw[order(ifelse(grepl("^temp", vars_raw), 0L, 1L), vars_raw)]
orders_avail   <- sort(na.omit(unique(app_df$orderName)))
families_avail <- sort(na.omit(unique(app_df$familyName)))
# IUCN filter currently disabled in app UI/server filtering; keep for easy restore:
# iucn_avail     <- sort(na.omit(unique(app_df$redlistCategory)))
default_year   <- if (2025 %in% years_avail) 2025 else max(years_avail)

species_lookup <- app_df %>%
  dplyr::select(spName, orderName, familyName) %>%
  distinct()

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

# ── helpers ───────────────────────────────────────────────────────────────────
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

safe_read_json <- function(url) {
  tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
}

BOUNDARY_CACHE_DIR <- file.path(intDir, "boundary_cache")
if (!dir.exists(BOUNDARY_CACHE_DIR)) dir.create(BOUNDARY_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

geoboundaries_country_index <- local({
  x <- safe_read_json("https://www.geoboundaries.org/api/current/gbOpen/ALL/ADM0/")
  if (is.null(x) || !is.data.frame(x)) return(data.frame())
  out <- x %>%
    dplyr::transmute(
      iso3 = as.character(boundaryISO),
      country = as.character(boundaryName)
    ) %>%
    filter(!is.na(iso3), nchar(iso3) == 3, !is.na(country), nzchar(country)) %>%
    distinct()
  out[order(out$country, out$iso3), , drop = FALSE]
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
                                       searchField = "label"))
      ),
      conditionalPanel("input.mode == 'hotspot'",
        h5(tags$span("Filter species (hotspot mode)", tags$span("?", class="tip", `data-tip`="Restrict which species count toward the hotspot map. Leave all blank to include all species."))),
        checkboxInput("flt_threatened", "Only threatened (CR/EN/VU)", value = FALSE),
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
        # IUCN filter temporarily disabled; keep code for easy restore.
        # if (length(iucn_avail))
        #   fluidRow(
        #     column(4, tags$label("IUCN status", style = "padding-top:7px;font-weight:600")),
        #     column(8, selectizeInput("flt_iucn", NULL,
        #                    choices = c("", iucn_avail), selected = "", multiple = TRUE,
        #                    options = list(placeholder = "All")))
        #   )
        # else helpText("IUCN filter unavailable — no metadata loaded."),
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
                  buttonLabel = "Browse…", placeholder = "No file selected")
      )
    ),
    mainPanel(width = 9,
      leafletOutput("map", height = "540px"),
      conditionalPanel("input.mode == 'hotspot'",
        div(style = "display:flex;justify-content:flex-end;margin-top:4px;margin-bottom:2px;",
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

  ensure_prop_exposed <- function(d) {
    if (!"propExposed" %in% names(d)) d$propExposed <- NA_real_
    d$propExposed <- suppressWarnings(as.numeric(d$propExposed))
    d
  }

  is_threatened <- function(x) {
    xx <- tolower(trimws(as.character(x)))
    xx %in% c("critically endangered", "endangered", "vulnerable", "cr", "en", "vu")
  }

  apply_hotspot_filters <- function(d, ord, fam, iuc, threatened = FALSE) {
    ord <- ord[nzchar(ord)]
    fam <- fam[nzchar(fam)]
    iuc <- iuc[nzchar(iuc)]
    if (length(ord)) d <- d %>% filter(orderName %in% ord)
    if (length(fam)) d <- d %>% filter(familyName %in% fam)
    if (length(iuc)) d <- d %>% filter(redlistCategory %in% iuc)
    if (isTRUE(threatened)) d <- d %>% filter(is_threatened(redlistCategory))
    d
  }

  updateSelectizeInput(session, "species",
                       choices = species_avail,
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

  avail_cells_all <- sort(unique(app_df$cell[!is.na(app_df$cell) & app_df$cell >= 1 & app_df$cell <= ncell(tpl)]))
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
    cur <- isolate(input$species)
    updateSelectizeInput(session, "species",
                         choices = spp,
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
    idx <- year_row_idx[[as.character(input$year)]]
    if (is.null(idx) || !length(idx)) return(app_df[0, , drop = FALSE])
    app_df[idx, , drop = FALSE]
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
                                 # input$flt_iucn %||% character(0),
                                 character(0),
                                 input$flt_threatened %||% FALSE)
    }
    d
  }) |> bindCache(input$year, input$sel_vars, input$mode,
                  input$flt_order, input$flt_family, input$flt_threatened)

  species_summary_df <- reactive({
    req(input$mode == "single", input$species)
    sel <- input$sel_vars
    req(length(sel))
    d <- app_df %>%
      filter(spName == input$species, var %in% sel)
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
    leaflet() %>% addProviderTiles(providers$Esri.WorldStreetMap) %>%
      setView(lng = 0, lat = 15, zoom = 2) %>%
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
          d   <- filtered_df()
          sel <- input$sel_vars %||% character(0)
          incProgress(0.2)

          show_progress(20, "Checking data...")
          if (!length(sel) || nrow(d) == 0) {
            exposure_msg_rv("No exposed cells for the selected year/filter/variable combination.")
            proxy %>% clearImages() %>% clearControls()
            hide_progress()
            set_status("Map update complete (no exposed cells)", done = TRUE)
            return()
          }

          if (input$mode == "single") {
            req(input$species)
            d_sp <- d %>% filter(spName == input$species)
            if (!nrow(d_sp)) {
              exposure_msg_rv(paste0("No exposed cells for '", input$species,
                                     "' in ", input$year,
                                     " with the selected variable(s)."))
              proxy %>% clearImages() %>% clearControls() %>% setView(lng = 0, lat = 15, zoom = 2)
              set_status("Species render complete (no exposed cells)", done = TRUE)
              return()
            }
            layers <- build_species_layers(d_sp, sel)
            if (!length(layers)) {
              exposure_msg_rv(paste0("No exposed cells for '", input$species,
                                     "' in ", input$year,
                                     " with the selected variable(s)."))
              proxy %>% clearImages() %>% clearControls() %>% setView(lng = 0, lat = 15, zoom = 2)
              map_raster_rv(NULL)
              set_status("Species render complete (no exposed cells)", done = TRUE)
              return()
            }
            exposure_msg_rv(NULL)
            proxy %>% clearImages() %>% clearControls()
            drawn <- character(0)
            # store cropped terra rasters for download
            terra_layers <- list()
            for (nm in names(layers)) {
              proxy %>% addRasterImage(layers[[nm]], colors = var_col(nm),
                                       opacity = 0.7, layerId = paste0("lyr_", nm),
                                       method = "ngb")
              drawn <- c(drawn, nm)
              # keep a terra version (full-extent, single layer)
              cv2 <- unique(d_sp$cell[as.character(d_sp$var) == nm])
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
            all_cells <- unique(d_sp$cell[!is.na(d_sp$cell) & d_sp$cell >= 1 & d_sp$cell <= ncell(tpl)])
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
          } else {
            exposure_msg_rv(NULL)
            rr <- build_hotspot_raster(d)
            if (is.null(rr)) {
              proxy %>% clearImages() %>% clearControls()
              set_status("Hotspot map complete (no exposed cells)", done = TRUE)
              return()
            }
            vals <- raster::values(rr)
            vals <- vals[is.finite(vals)]
            if (!length(vals)) {
              proxy %>% clearImages() %>% clearControls()
              set_status("Hotspot map complete (no exposed cells)", done = TRUE)
              return()
            }
            minN <- suppressWarnings(min(vals, na.rm = TRUE))
            maxN <- suppressWarnings(max(vals, na.rm = TRUE))
            if (!is.finite(minN) || !is.finite(maxN) || maxN == 0) {
              proxy %>% clearImages() %>% clearControls()
              set_status("Hotspot map complete (no exposed cells)", done = TRUE)
              return()
            }
            if (minN == maxN) {
              minN <- minN - 0.5
              maxN <- maxN + 0.5
            }
            n_cols <- max(2, min(100, as.integer(maxN - minN + 1)))
            cm_cols <- colorRampPalette(c("steelblue4", "steelblue1", "gold", "red1", "red4"))(n_cols)
            pal <- colorNumeric(cm_cols, domain = c(minN, maxN), na.color = "transparent")
            proxy %>% clearImages() %>% clearControls()
            proxy %>%
              addRasterImage(rr, colors = pal, opacity = input$hotspot_opacity %||% 0.95, layerId = "hotspot", method = "ngb") %>%
              addLegend(position = "bottomright", pal = pal, values = vals,
                        opacity = input$hotspot_opacity %||% 0.95,
                        title = "# Species<br>exposed", layerId = "legend_main")
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
    sp <- input$species %||% ""
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
      app_df %>%
        filter(spName == input$species, var %in% sel) %>%
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
      d_sp <- filtered_df() %>% filter(spName == input$species)
      if (!nrow(d_sp)) return(NULL)
      sel <- input$sel_vars %||% character(0)
      layers <- list()
      for (v in intersect(sel, as.character(unique(d_sp$var)))) {
        cv <- unique(d_sp$cell[as.character(d_sp$var) == v])
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

  output$download_map_raster <- downloadHandler(
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
    updateRadioButtons(session, "mode", selected = "hotspot")
    updateSelectInput(session, "year", selected = default_year)
    updateCheckboxGroupInput(session, "sel_vars", selected = vars_avail)
    updateCheckboxInput(session, "flt_exposed_only", value = FALSE)
    updateSelectizeInput(session, "species", selected = species_avail[1], server = TRUE)
    updateCheckboxInput(session, "flt_threatened", value = FALSE)
    updateSelectizeInput(session, "flt_order", selected = "", server = TRUE)
    updateSelectizeInput(session, "flt_family", selected = "", server = TRUE)
    # updateSelectizeInput(session, "flt_iucn", selected = "", server = TRUE)
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

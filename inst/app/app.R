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
iucn_avail     <- sort(na.omit(unique(app_df$redlistCategory)))
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
var_choices <- setNames(vars_avail, vapply(vars_avail, var_label, character(1)))

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

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML(
    "body{font-size:13px}
     #click_panel{margin-top:12px;padding:10px;border:1px solid #ddd;border-radius:4px;background:#fafafa}
     #click_panel h4{margin-top:4px}
     #polygon_panel{margin-top:12px;padding:10px;border:1px solid #e67e00;border-radius:4px;background:#fffaf5}
     #polygon_panel h4{margin-top:4px}"
  ))),

  titlePanel("Species exposure to extreme climate"),
  sidebarLayout(
    sidebarPanel(width = 3,
      radioButtons("mode", "Map mode",
                   choices  = c("Single species" = "single",
                                "Hotspot: # species/cell" = "hotspot"),
                   selected = "hotspot"),
      hr(),
      selectInput("year", "Year", choices = years_avail, selected = default_year),
      checkboxGroupInput("sel_vars", "Climate variables", choices = var_choices, selected = vars_avail),
      fluidRow(
        column(6, actionButton("vars_all", "Select all", width = "100%")),
        column(6, actionButton("vars_none", "Deselect all", width = "100%"))
      ),
      hr(),
      conditionalPanel("input.mode == 'single'",
        h5("Species"),
        checkboxInput("flt_exposed_only", "Only show species with exposure in selected year/variable(s)", value = FALSE),
        selectizeInput("species", NULL,
                       choices  = NULL,
                       selected = character(0),
                       options  = list(maxItems = 1, placeholder = "Type species name...",
                                       searchField = "label"))
      ),
      conditionalPanel("input.mode == 'hotspot'",
        h5("Filter species (hotspot mode)"),
        checkboxInput("flt_threatened", "Only threatened (CR/EN/VU)", value = FALSE),
        if (length(orders_avail))
          selectizeInput("flt_order", "Order",
                         choices = c("", orders_avail), selected = "", multiple = TRUE,
                         options = list(placeholder = "All orders"))
        else helpText("Order filter unavailable — no metadata loaded."),
        if (length(families_avail))
          selectizeInput("flt_family", "Family",
                         choices = c("", families_avail), selected = "", multiple = TRUE,
                         options = list(placeholder = "All families"))
        else helpText("Family filter unavailable — no metadata loaded."),
        if (length(iucn_avail))
          selectizeInput("flt_iucn", "IUCN status",
                         choices = c("", iucn_avail), selected = "", multiple = TRUE,
                         options = list(placeholder = "All categories"))
        else helpText("IUCN filter unavailable — no metadata loaded.")
      )
    ),
    mainPanel(width = 9,
      leafletOutput("map", height = "540px"),
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
        div(id = "polygon_panel",
            h4("Species in drawn polygon"),
            p(style = "font-size:12px;color:#555",
              "Use the polygon tool (top-left map toolbar) to select a region. ",
              "Shows all unique species exposed for the current year & variables within the polygon."),
            actionButton("clear_polygon", "Clear polygon", class = "btn-sm"),
            br(), br(),
            DTOutput("polygon_table"),
            br(),
            downloadButton("download_polygon",      "Download polygon table"),
            downloadButton("download_map_raster_hs", "Download map raster")
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  exposure_msg_rv <- reactiveVal(NULL)
  status_msg_rv <- reactiveVal("Ready")

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
                                 input$flt_iucn %||% character(0),
                                 input$flt_threatened %||% FALSE)
    }
    d
  }) |> bindCache(input$year, input$sel_vars, input$mode,
                  input$flt_order, input$flt_family, input$flt_iucn, input$flt_threatened)

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
      addDrawToolbar(
        targetGroup          = "drawn",
        polylineOptions      = FALSE,
        rectangleOptions     = FALSE,
        circleOptions        = FALSE,
        markerOptions        = FALSE,
        circleMarkerOptions  = FALSE,
        polygonOptions = drawPolygonOptions(
          showArea     = FALSE,
          shapeOptions = drawShapeOptions(fillOpacity = 0.15, color = "#e67e00", weight = 2)
        ),
        editOptions = editToolbarOptions(remove = TRUE, edit = FALSE)
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
              addRasterImage(rr, colors = pal, opacity = 0.75, layerId = "hotspot", method = "ngb") %>%
              addLegend(position = "bottomright", pal = pal, values = vals,
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
    req(!is.null(cell_ids), length(cell_ids) > 0)
    fd <- filtered_df()
    d  <- fd %>% filter(cell %in% cell_ids)
    d  <- ensure_prop_exposed(d)
    req(nrow(d) > 0)
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

  # Clear stale click/polygon state whenever the user switches mode
  observeEvent(input$mode, {
    click_info(list(cell = NA_integer_, lat = NA, lng = NA))
    click_data(NULL)
    hotspot_diag_data(NULL)
    polygon_cells_rv(NULL)
    map_raster_rv(NULL)
    # Remove drawn polygons from the map
    leafletProxy("map") %>% clearGroup("drawn")
  }, ignoreInit = TRUE)

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
    d <- tryCatch(polygon_table_data(), error = function(e) NULL)
    req(!is.null(d))
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

  # Helper: reset the draw toolbar so it no longer intercepts map clicks
  reset_draw_toolbar <- function() {
    leafletProxy("map") %>%
      clearGroup("drawn") %>%
      removeDrawToolbar(clearFeatures = TRUE) %>%
      addDrawToolbar(
        targetGroup         = "drawn",
        polylineOptions     = FALSE,
        rectangleOptions    = FALSE,
        circleOptions       = FALSE,
        markerOptions       = FALSE,
        circleMarkerOptions = FALSE,
        polygonOptions = drawPolygonOptions(
          showArea     = FALSE,
          shapeOptions = drawShapeOptions(fillOpacity = 0.15, color = "#e67e00", weight = 2)
        ),
        editOptions = editToolbarOptions(remove = TRUE, edit = FALSE)
      )
  }

  observeEvent(input$clear_polygon, {
    polygon_cells_rv(NULL)
    reset_draw_toolbar()
  })

  # Also handle delete via the toolbar's own delete button
  observeEvent(input$map_draw_deleted_features, {
    polygon_cells_rv(NULL)
    reset_draw_toolbar()
  })

  observeEvent(input$map_draw_new_feature, ignoreNULL = TRUE, {
    req(input$mode == "hotspot")
    feat <- input$map_draw_new_feature
    req(!is.null(feat))
    # hide click panels while polygon is active
    click_data(NULL)
    hotspot_diag_data(NULL)
    click_info(list(cell = NA_integer_, lat = NA, lng = NA))
    tryCatch({
      geojson_str <- jsonlite::toJSON(feat, auto_unbox = TRUE)
      tmp_geojson <- tempfile(fileext = ".geojson")
      writeLines(as.character(geojson_str), tmp_geojson)
      poly_vect <- terra::vect(tmp_geojson)

      # Compute polygon cell IDs against ALL cells ever in app_df (year-independent)
      # so the table can reactively refilter when year/vars change.
      avail_cells <- unique(app_df$cell)
      avail_cells <- avail_cells[!is.na(avail_cells) & avail_cells >= 1 &
                                    avail_cells <= ncell(tpl)]
      if (!length(avail_cells)) { polygon_cells_rv(NULL); return() }

      xy  <- terra::xyFromCell(tpl, avail_cells)
      pts <- terra::vect(xy, type = "points", crs = terra::crs(tpl))
      inside <- tryCatch(
        terra::is.related(pts, poly_vect, "intersects"),
        error = function(e) rep(FALSE, length(avail_cells))
      )
      cell_ids <- avail_cells[inside]
      if (!length(cell_ids)) { polygon_cells_rv(NULL); return() }

      polygon_cells_rv(cell_ids)
      set_status(sprintf("Polygon: %d cells selected", length(cell_ids)), done = TRUE)
    }, error = function(e) {
      polygon_cells_rv(NULL)
      exposure_msg_rv(paste("Polygon query failed:", conditionMessage(e)))
    })
  })
}

shinyApp(ui, server)

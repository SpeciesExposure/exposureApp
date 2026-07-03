#' Launch the Species Climate Exposure Shiny App
#'
#' @param intDir Path to the directory containing the app data files
#'   (sharded \code{allExpForShiny_by_year/}, \code{allExpForShiny_manifest_v1.qs},
#'   \code{landTemplate.tif}, etc.).
#'   Defaults to the package's built-in \code{extdata} directory.
#' @param launch.browser Logical; open the app in the default browser. Default \code{TRUE}.
#' @param host Host IP address to listen on. Default \code{"127.0.0.1"}.
#' @param port Port number. Default uses \code{getOption("shiny.port")}.
#' @return Launches a Shiny application (does not return a value).
#' @export
exposureApp <- function(intDir = NULL,
                        launch.browser = TRUE,
                        host = "127.0.0.1",
                        port = getOption("shiny.port")) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required.")
  }

  if (is.null(intDir)) {
    intDir <- system.file("extdata", package = "exposureApp")
  }

  if (!nzchar(intDir) || !dir.exists(intDir)) {
    stop("Could not find app data directory. Provide intDir explicitly.")
  }

  Sys.setenv(INT_DIR = intDir)

  app_dir <- system.file("app", package = "exposureApp")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Could not find packaged app directory.")
  }

  shiny::runApp(appDir = app_dir,
                launch.browser = launch.browser,
                host = host,
                port = port)
}

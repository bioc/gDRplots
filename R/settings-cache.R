.settings_cache <- new.env(parent = emptyenv())

#' @keywords internal
.get_setting <- function(key) {
  if (is.null(.settings_cache[[key]])) {
    .settings_cache[[key]] <- gDRutils::get_settings_from_json(
      key, system.file(package = "gDRplots", "settings.json")
    )
  }
  .settings_cache[[key]]
}

#' Get \code{ggplot2::theme_bw()}, built once per session
#'
#' @return \code{ggplot2::theme_bw()}
#' @keywords internal
.theme_bw <- function() {
  if (is.null(.settings_cache[["theme_bw"]])) {
    .settings_cache[["theme_bw"]] <- ggplot2::theme_bw()
  }
  .settings_cache[["theme_bw"]]
}

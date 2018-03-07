## User guide launcher
#' @title View ASICS User's Guide
#' @name ASICSUsersGuide
#' @export
#'
#' @importFrom utils browseURL
#'
#' @description Open the ASICS User's Guide (with default browser)
#'
#' @param view if \code{TRUE}, the user's guide will be opened with default
#' browser
#' @details The function \code{vignette("ASICS")} will find the short
#' ASICS vignette that describes the main functions and how to obtain the ASICS
#' User's Guide.
#' The User's Guide is not itself a true vignette because it is not
#' automatically generated during the package build process.
#' If the operating system is not Windows, then the HTML viewer used is the one
#' given by \code{Sys.getenv("R_BROWSER")}. The HTML viewer can be changed using
#' \code{Sys.setenv(R_BROWSER = )}.
#'
#' @return Character string giving the file location. If \code{view = TRUE},
#' the HTML viewer is started and the User’s Guide is opened, as a side effect.
#'
#' @examples
#' # To get the location
#' ASICSUsersGuide(view = FALSE)
#'
#' # To open in a HTML viewer
#' \dontrun{ASICSUsersGuide()}

ASICSUsersGuide <- function(view=TRUE) {
  f <- system.file("doc", "ASICSUsersGuide.html", package = "ASICS")
  if(view) {
    if(.Platform$OS.type == "windows")
      shell.exec(f)
    else
      browseURL(paste0("file://", f))
  }
  return(f)
}

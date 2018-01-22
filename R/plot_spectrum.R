#' Plot original and recomposed spectra
#'
#' Plot original spectrum and recomposed one obtained from ASICS quantification.
#' It is possible to add one pure metabolite spectrum.
#'
#' @param ASICS.results an object of class \linkS4class{ASICSResults}.
#' @param idx index of the spectrum to plot. Default to 1.
#' @param xlim,ylim boundaries for x and y, respectively.
#' @param pure.library pure library used for quantification. Default to
#' \code{NULL} (the library included in the package is used).
#' @param add.metab name of one metabolite to add to the plot. Default to
#' \code{NULL} (no pure spectrum added to the plot).
#' @return a ggplot plot original and recomposed spectra of one sample
#' in one figure. In addition, one pure metabolite spectrum (as
#' provided in the reference library) and the deformed one can be superimposed
#' to the plot.
#' @importFrom ggplot2 ggplot aes_string geom_line theme_bw labs
#' @importFrom ggplot2 scale_x_reverse ylim scale_color_manual
#' @export
plot_spectrum <- function(ASICS.results, idx = 1, xlim = c(0, 10), ylim = NULL,
                          pure.library = NULL, add.metab = NULL) {

  # Original spectrum and reconstituted one
  spectra <- data.frame(grid = rep(ASICS.results@ppm.grid, 2),
                        mixture = c(ASICS.results@spectra[, idx],
                                    ASICS.results@recomposed.spectra[, idx]),
                        which_mix  = c(rep("Original spectrum",
                                           length(ASICS.results@ppm.grid)),
                                       rep("Recomposed spectrum",
                                           length(ASICS.results@ppm.grid))))

  if (is.null(ylim)) ylim <- c(0, max(spectra$mixture))

  if (is.null(pure.library)) {
    pure_lib <- ASICS::pure_library
  } else {
    pure_lib <- pure.library
  }

  # Add a pure spectrum if the user wants
  if (!is.null(add.metab) && add.metab %in% pure_lib@sample.name) {

    # original pure spectra
    metab_to_add_pure <-
      pure_lib@spectra[, pure_lib@sample.name == add.metab]
    metab_to_add_pure <- metab_to_add_pure *
      max(spectra$mixture) / max(metab_to_add_pure)

    spectra <-
      rbind(spectra,
            data.frame(grid = ASICS.results@ppm.grid,
                       mixture = metab_to_add_pure,
                       which_mix  = rep(add.metab,
                                        length(ASICS.results@ppm.grid))))

    if (add.metab %in% ASICS.results@deformed.library$metabolite_name[
      ASICS.results@deformed.library$sample == ASICS.results@sample.name[idx]]) {

      metab_to_add_def_without0 <- ASICS.results@deformed.library[
        ASICS.results@deformed.library$sample == ASICS.results@sample.name[idx] &
          ASICS.results@deformed.library$metabolite_name == add.metab , ]

      metab_to_add_def <- merge(data.frame(ppm_grid = as.character(ASICS.results@ppm.grid)),
                                metab_to_add_def_without0[, c("intensity", "ppm_grid")],
                                all.x = TRUE)
      metab_to_add_def[is.na(metab_to_add_def)] <- 0

      spectra <-
        rbind(spectra,
              data.frame(grid = ASICS.results@ppm.grid,
                         mixture = metab_to_add_def$intensity,
                         which_mix  = rep(paste(add.metab, "deformed"),
                                          length(ASICS.results@ppm.grid))))
    }
  }

  spectra$which_mix <- relevel(spectra$which_mix, "Original spectrum")
  # Plot with ggplot2
  p1 <- ggplot(spectra) +
    geom_line(aes_string(x = "grid", y = "mixture",
                         colour = "which_mix"), na.rm = TRUE) +
    theme_bw() +
    labs(x = "Chemical shift (ppm)", y = "Intensity") +
    scale_x_reverse(limits = rev(xlim)) +
    ylim(ylim) +
    scale_color_manual(name = "", values = c("black", "blue", "red", "orange"))

  return(p1)
}

#' @importFrom stats relevel
.plotSpectrum <- function(ASICS.results, idx = 1, xlim = c(0, 10), ylim = NULL,
                          pure.library = NULL, add.metab = NULL) {

  # original spectrum and reconstituted one
  spectra <- data.frame(grid = rep(ASICS.results@ppm.grid, 2),
                        mixture = c(ASICS.results@spectra[, idx],
                                    ASICS.results@reconstructed.spectra[, idx]),
                        which_mix  = c(rep("Original spectrum",
                                           length(ASICS.results@ppm.grid)),
                                       rep("Reconstructed spectrum",
                                           length(ASICS.results@ppm.grid))))

  if (is.null(ylim)) ylim <- c(0, max(spectra$mixture))

  if (is.null(pure.library)) {
    pure_lib <- ASICS::pure_library
  } else {
    pure_lib <- pure.library
    pure_lib@spectra <-
      pure_lib@spectra[pure_lib@ppm.grid >= min(ASICS.results@ppm.grid) &
                         pure_lib@ppm.grid <= max(ASICS.results@ppm.grid), ]
    pure_lib@ppm.grid <-
      pure_lib@ppm.grid[pure_lib@ppm.grid >= min(ASICS.results@ppm.grid) &
                          pure_lib@ppm.grid <= max(ASICS.results@ppm.grid)]
  }

  # add a pure spectrum if the user wants
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
      ASICS.results@deformed.library$sample ==
      ASICS.results@sample.name[idx]]) {

      metab_to_add_def_without0 <- ASICS.results@deformed.library[
        ASICS.results@deformed.library$sample ==
          ASICS.results@sample.name[idx] &
          ASICS.results@deformed.library$metabolite_name == add.metab , ]

      metab_to_add_def <-
        merge(data.frame(ppm_grid = as.character(ASICS.results@ppm.grid)),
              metab_to_add_def_without0[, c("intensity", "ppm_grid")],
              all.x = TRUE)
      metab_to_add_def[is.na(metab_to_add_def)] <- 0

      spectra <-
        rbind(spectra,
              data.frame(grid = ASICS.results@ppm.grid,
                         mixture = metab_to_add_def$intensity,
                         which_mix  = rep(paste("Pre-processed", add.metab),
                                          length(ASICS.results@ppm.grid))))
    }
  }

  spectra$which_mix <- relevel(spectra$which_mix, "Original spectrum")
  # plot with ggplot2
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

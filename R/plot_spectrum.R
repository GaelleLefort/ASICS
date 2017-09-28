## Plot ASICS results
plot_spectrum <- function(res_ASICS, xmin = 0, xmax = 10, ymin = 0, ymax = NULL,
                         add_metab = NULL){

  # Original spectrum and reconstituted one
  spectra <- data.frame(grid = rep(res_ASICS@ppm_grid, 2),
                        mixture = c(res_ASICS@original_mixture,
                                    res_ASICS@reconstituted_mixture),
                        which_mix  = c(rep("Original spectrum",
                                           length(res_ASICS@ppm_grid)),
                                       rep("Estimated spectrum",
                                           length(res_ASICS@ppm_grid))))

  if(is.null(ymax)) ymax <- max(spectra$mixture)

  # Add a pure spectrum if the user wants
  if(!is.null(add_metab)){
    metab_to_add <- ASICS::pure_library$spectra[, ASICS::pure_library$name == add_metab]
    metab_to_add <- metab_to_add * max(spectra$mixture) / max(metab_to_add)
    spectra <- rbind(spectra, data.frame(grid = rep(res_ASICS@ppm_grid, 2),
                                         mixture = metab_to_add,
                                         which_mix  = rep(add_metab,
                                                          length(res_ASICS@ppm_grid))))
  }

  spectra$which_mix <- relevel(spectra$which_mix, "Original spectrum")
  # Plot with ggplot2
  p1 <- ggplot2::ggplot(spectra) +
    ggplot2::geom_line(ggplot2::aes_string(x = "grid", y = "mixture",
                                           colour = "which_mix")) +
    ggplot2::theme_bw() +
    ggplot2::xlab("ppm") +
    ggplot2::ylab("Standardized intensity") +
    ggplot2::scale_x_reverse(limits = c(xmax, xmin)) +
    ggplot2::ylim(c(ymin, ymax)) +
    ggplot2::scale_color_manual(name = "", values = c("black", "blue", "red"))

  print(p1)
}

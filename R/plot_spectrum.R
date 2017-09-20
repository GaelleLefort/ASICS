#' Plot ASICS results
#'
#' Description
#' @param res_ASICS result of ASICS function
#' @param xmin,xmax x minimum and maximum limits
#' @param ymin,ymax y minimum and maximum limits
#' @param add_metab name of one metabolite to add to the plot
#' @return Produces a plot of original and recontituted spectra.
#' @export
#' @examples
#' #Compute quantification on one spectrum
#' result <- ASICS(path = system.file("extdata", "example_spectra",
#'                                    "AG_faq_Beck01", package = "ASICS"),
#'                 exclusion.areas = matrix(c(4.5,5.1,5.5,6.5),
#'                                          ncol = 2, byrow = TRUE),
#'                 max.shift = 0.02, which.spectra = "last",
#'                 library.metabolites = NULL, threshold.noise = 0.02)
#'
#' #Produces a plot of original, recontituted and Lactate spectra
#' plot_spectrum(result, xmin = 1, xmax = 1.5, ymax = 10, add_metab = "Lactate")

plot_spectrum <- function(res_ASICS, xmin = 0, xmax = 10, ymin = 0, ymax = NULL,
                         add_metab = NULL){

  # Original spectrum and reconstituted one
  spectra <- data.frame(grid = rep(res_ASICS$Grid, 2),
                        mixture = c(res_ASICS$Mixture,
                                    res_ASICS$Estimated_mixture),
                        which_mix  = c(rep("Original spectrum",
                                           length(res_ASICS$Grid)),
                                       rep("Estimated spectrum",
                                           length(res_ASICS$Grid))))

  if(is.null(ymax)) ymax <- max(spectra$mixture)

  # Add a pure spectrum if the user wants
  if(!is.null(add_metab)){
    spectra <- rbind(spectra, data.frame(grid = rep(res_ASICS$Grid, 2),
                                         mixture = pure_library$spectra[,
                                              pure_library$name == add_metab],
                                         which_mix  = rep(add_metab,
                                             length(res_ASICS$Grid))))
  }

  spectra$which_mix <- relevel(spectra$which_mix, "Original spectrum")
  # Plot with ggplot2
  p1 <- ggplot2::ggplot(spectra) +
    ggplot2::geom_line(ggplot2::aes(x = grid, y = mixture, colour = which_mix)) +
    ggplot2::theme_bw() +
    ggplot2::xlab("ppm") +
    ggplot2::ylab("Standardized intensity") +
    ggplot2::scale_x_reverse(limits = c(xmax, xmin)) +
    ggplot2::ylim(c(ymin, ymax)) +
    ggplot2::scale_color_manual(name = "", values = c("black", "blue", "red"))

  return(p1)
}

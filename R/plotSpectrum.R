#' Plot a spectrum
#'
#' Description
#' @param resASICS result of ASICS function
#' @param xmin x minimum limit
#' @param xmax x maximum limit
#' @param ymin y minimum limit
#' @param ymax y maximum limit
#' @param DecalageLib chemical shift between library of pure spectra and complex mixture
#' @param addMetab name of one metabolite that we want to add the pure spectra
#' @keywords NMR quantification metabolites plot
#' @export
#' @examples
#' plotSpectrum(result, xmin = 1, xmax = 1.5, ymax = 10, addMetab = "Lactate",
#'  DecalageLib = 0)
plotSpectrum <- function(resASICS, xmin = 0, xmax = 10, ymin = 0, ymax = NULL,
                         DecalageLib = 0, addMetab = NULL){
  spectra <- data.frame(grid = resASICS$Grid, mixture = resASICS$Mixture,
                        est_mixture  = resASICS$Estimated_mixture)

  if(is.null(ymax)) ymax <- max(c(spectra$mixture, spectra$est_mixture))

  p1 <- ggplot(spectra) + geom_line(aes(x = grid, y = mixture, color = "Original spectrum")) +
    geom_line(aes(x = grid, y = est_mixture, color = "Estimated mixture")) +
    theme_bw() + xlab("ppm") + ylab("Standardized intensity") +
    scale_x_reverse(limits = c(xmax, xmin)) + ylim(c(ymin, ymax)) +
    scale_color_manual(name = "", values=c("Original spectrum"="black", "Estimated mixture"="blue"),
                       guide = guide_legend(reverse=TRUE))


  if(!is.null(addMetab)){
    p1 <- p1 + geom_line(aes(x = Library$Grid + DecalageLib,
                             y = Library$Metab[, Library$Name == addMetab],
                             color = addMetab)) +
      scale_color_manual(name = "", values=c("black", "blue", "red"),
                         guide = guide_legend(reverse=TRUE))
  }

  return(p1)
}

#' Tile plot
#'
#' Tile plot of spectra to see if an alignment is needed or the result of an
#' alignment.
#'
#' @param spectra_obj An object of class \linkS4class{Spectra} obtained with the
#' function \link{createSpectra}.
#' @param xlim Boundaries for x.
#'
#' @return
#' A \code{\link[ggplot2]{ggplot}} plot of original and reconstructed
#' spectra of one sample in the same figure for \linkS4class{ASICSResults}
#' object. In addition, one pure metabolite spectrum (as provided in the
#' reference library) and the deformed one can be superimposed to the plot.
#'
#' @examples
#' # Import data and create object
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- importSpectraBruker(current_path)
#' spectra_obj <- createSpectra(spectra_data)
#' plotAlignment(spectra_obj, xlim = c(3,4))
#'
#'
#' @seealso \code{\link{alignSpectra}}
#' @importFrom stats reshape
#' @export
plotAlignment <- function(spectra_obj, xlim = c(0, 10)) {
  # reshape data
  data_to_plot <- reshape(as.data.frame(as.matrix(spectra_obj@spectra)),
                          varying = seq_len(ncol(spectra_obj@spectra)),
                          times = spectra_obj@sample.name,
                          ids = as.numeric(spectra_obj@ppm.grid),
                          v.names = "intensity", timevar = "sample_name",
                          idvar = "ppm_grid", direction = "long",
                          new.row.names = seq_len(ncol(spectra_obj@spectra) *
                                                    nrow(spectra_obj@spectra)))
  data_to_plot$intensity <- ifelse(data_to_plot$intensity > 0,
                                   data_to_plot$intensity, 0)
  data_to_plot$log_intensity <- log(data_to_plot$intensity)

  # plot
  graph <- ggplot(data_to_plot, aes_string(x = "ppm_grid", y = "sample_name")) +
    geom_tile(aes_string(fill = "log_intensity"), na.rm = TRUE) + theme_bw() +
    scale_x_reverse(limits = rev(xlim)) +
    labs(x = "Chemical shift (ppm)", y = "Sample names") +
    scale_fill_distiller(type = "div", palette = "Spectral",
                         name = "log intensity")

  return(graph)
}


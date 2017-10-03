#' Pure spectra library
#'
#' The 1D 1H NMR spectra of 175 reference compounds were collected to build the
#' spectral library. These compounds have been prepared and recorded using a
#' Bruker Avance III HD spectrometer in the MetaToul - AXIOM Site at Toulouse
#' (France).
#'
#' @name pure_library
#'
#' @docType data
#'
#' @format A list of 4 elements:
#' \describe{
#'   \item{name}{names of the metabolites}
#'   \item{grid}{common grid for all spectra}
#'   \item{spectra}{a data frame with each pure metabolite spectrum in column}
#'   \item{nb_protons}{number of protons of each metabolite}
#' }
#'
#' @references Tardivel P., Canlet C., Lefort G., Tremblay-Franco M., Debrauwer 
#' L., Concordet D., Servien R. (2017). ASICS: an automatic method for 
#' identification and quantification of metabolites in complex 1D 1H NMR 
#' spectra. \emph{Metabolomics}, \strong{13}(10): 109.
#' \url{https://doi.org/10.1007/s11306-017-1244-5}

NULL


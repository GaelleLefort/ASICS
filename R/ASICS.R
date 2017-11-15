#' Automatic Statistical Identification in Complex Spectra
#'
#' Quantification of 1D 1H NMR spectra with ASICS method using a library of 175
#' pure metabolite spectra. The method is presented in Tardivel et al. (2017).
#'
#' @param spectra_obj an object of class \linkS4class{Spectra} obtained with the
#' function \link{preprocessing_spectra}
#' @param exclusion.areas definition domain of spectra to exclude (ppm)
#' @param max.shift maximum chemical shift allowed (in ppm)
#' @param pure.library an object of class \linkS4class{PureLibrary} containing
#' the references (pure metabolite spectra). If \code{NULL}, the library
#' included in the package is used
#' @param threshold.noise threshold for signal noise
#' @param seed random seed to control randomness in the algorithm (used in the
#' estimation of significativity of a given metabolite concentration)
#'
#' @return An object of type \code{\link{ASICSResults}}.
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
#'
#' @seealso \code{\link{resASICS-class}} \code{\link{ASICS_multiFiles}}
#' \code{\link{pure_library}}
#'
#' @references Tardivel P., Canlet C., Lefort G., Tremblay-Franco M., Debrauwer
#' L., Concordet D., Servien R. (2017). ASICS: an automatic method for
#' identification and quantification of metabolites in complex 1D 1H NMR
#' spectra. \emph{Metabolomics}, \strong{13}(10): 109.
#' \url{https://doi.org/10.1007/s11306-017-1244-5}
#'

ASICS <- function(spectra_obj,
                  exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                  max.shift = 0.02, pure.library = NULL,
                  threshold.noise = 0.02, seed = 1234) {

  if(!is.matrix(exclusion.areas) | ncol(exclusion.areas) != 2){
    stop("'exclusion.areas' needs to be a matrix with 2 columns.")
  }

  if(max.shift < 0){
    stop("'max.shift' must be non negative.")
  }

  if(threshold.noise < 0){
    stop("'threshold.noise' must be non negative.")
  }

  if(class(pure.library) != "PureLibrary" & !is.null(pure.library)){
    stop(paste("'pure.library' needs to be either NULL or an object of class",
               "'PureLibrary'."))
  }

  list_spec <- lapply(seq_along(spectra_obj), function(x) spectra_obj[x])

  res_estimation_list <- bplapply(list_spec,
                          ASICS.internal, exclusion.areas, max.shift,
                          pure.library, threshold.noise, seed,
                          BPPARAM = MulticoreParam(progressbar = TRUE,
                                                   tasks = length(spectra_obj)))

  res_estimation <- do.call(c, res_estimation_list)

  return(res_estimation)

}


#' @importFrom methods new
ASICS.internal <- function(spectrum_obj,
                           exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                           max.shift = 0.02, pure.library = NULL,
                           threshold.noise = 0.02, seed = 1234){

  ##Seed
  set.seed(seed)


  #-----------------------------------------------------------------------------
  #### Remove areas from spectrum and library ####
  cleaned_obj <- remove_areas(spectrum_obj, exclusion.areas = exclusion.areas,
                              pure.library = pure.library)

  cleaned_spectrum <- cleaned_obj$cleaned_spectrum
  cleaned_library <- cleaned_obj$cleaned_library

  #number of points on library grid corresponding to maximum shift
  nb_points_shift <- floor(max.shift / (cleaned_library@ppm.grid[2] -
                                          cleaned_library@ppm.grid[1]))


  #-----------------------------------------------------------------------------
  #### Cleaning step: remove metabolites that cannot belong to the mixture ####
  cleaned_library <- clean_library(cleaned_spectrum, cleaned_library,
                                   threshold.noise, nb_points_shift)


  #-----------------------------------------------------------------------------
  #### Find the best translation between each pure spectra and mixture ####
  #and sort metabolites by regression residuals

  #Compute weights
  s1 <- 0.172 #standard deviation of multiplicative noise
  s2 <- 0.15 #standard deviation of additive noise
  noises <- abs(cleaned_spectrum@spectra) * s1 ^ 2 + s2 ^ 2
  mixture_weights <- as.numeric(1 / noises)

  res_translation <- translate_library(cleaned_spectrum, cleaned_library,
                                       mixture_weights, nb_points_shift)

  sorted_library <- res_translation$sorted_library
  shift <- res_translation$shift


  #-----------------------------------------------------------------------------
  #### Localized deformations of pure spectra ####
  deformed_library <- deform_library(cleaned_spectrum, sorted_library,
                                     mixture_weights, nb_points_shift,
                                     max.shift, shift)


  #-----------------------------------------------------------------------------
  #### Threshold and concentration optimisation for each metabolites ####
  system.time(res_opti <- concentration_opti(cleaned_spectrum, deformed_library,
                                 noises, mixture_weights))
  final_library <- res_opti$final_library


  #-----------------------------------------------------------------------------
  #### Results ####

  #Reconstituted mixture with estimated coefficents and
  est_mixture <- final_library@spectra %*% res_opti$B_final
  pure_lib_final_conc <- final_library
  pure_lib_final_conc@spectra <- final_library@spectra %*%
    diag(res_opti$B_final)

  #Compute relative concentration of identified metabolites
  relative_concentration <- res_opti$B_final / final_library@nb.protons
  #Sort library according to relative concentration
  sorted_idx <- sort(relative_concentration, decreasing = TRUE,
                     index.return = TRUE)$ix
  pure_lib_final_sorted <- pure_lib_final_conc[sorted_idx]

  present_metab <- data.frame(Metabolites = pure_lib_final_sorted@sample.name,
                              relative_concentration[sorted_idx])
  colnames(present_metab)[2] <- spectrum_obj@sample.name


  #List to return
  res_object <- new(Class = "ASICSResults",
                    sample.name = spectrum_obj@sample.name,
                    ppm.grid = spectrum_obj@ppm.grid,
                    spectra = spectrum_obj@spectra,
                    recomposed.spectra = est_mixture,
                    quantification = present_metab,
                    deformed.library = list(pure_lib_final_sorted))

  return(res_object)
}

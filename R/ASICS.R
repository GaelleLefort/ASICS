#' Automatic Statistical Identification in Complex Spectra
#'
#' Quantification of 1D 1H NMR spectrum with ASICS method using a library of 175
#' pure metabolite spectra. The method is presented in Tardivel et al. (2017).
#'
#' @param path folder path of the Bruker files
#' @param exclusion.areas definition domain of spectra to exclude (ppm)
#' @param max.shift maximum chemical shift allowed (in ppm)
#' @param which.spectra if more than one spectra by sample, spectra used to
#' perform the quantification (either \code{"first"}, \code{"last"} or its
#' number). Default to \code{"last"}
#' @param library.metabolites path of the library containing the references
#' (pure metabolite spectra). If \code{NULL}, the library included in the
#' package is used
#' @param threshold.noise threshold for signal noise
#' @param seed random seed to control randomness in the algorithm (used in the
#' estimation of significativity of a given metabolite concentration)
#' @param nb.iter.signif number of iterations for the estimation of
#' significativity of a given metabolite concentration. Default to 400
#' 
#' @return An object of type \code{\link{resASICS-class}}
#' 
#' @importFrom methods new
#' @importFrom stats relevel
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
#' @examples
#' \dontshow{
#' lib_file <- system.file("extdata", "library_for_examples.rda",
#'                         package = "ASICS")
#' cur_path <- system.file("extdata", "example_spectra", "AG_faq_Beck01",
#'                         package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' result <- ASICS(path = cur_path, exclusion.areas = to_exclude,
#'                 nb.iter.signif = 10, library.metabolites = lib_file)
#' }
#' \dontrun{
#' cur_path <- system.file("extdata", "example_spectra", "AG_faq_Beck01",
#'                         package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' result <- ASICS(path = cur_path, exclusion.areas = to_exclude)
#' }


ASICS <- function(path, exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                  max.shift = 0.02, which.spectra = "last",
                  library.metabolites = NULL, threshold.noise = 0.02,
                  seed = 1234, nb.iter.signif = 400){

  ## checking the validity of the parameters
  if(!dir.exists(path)){
    stop("Path of the Bruker file doesn't exist!")
  }

  if(!is.matrix(exclusion.areas) | ncol(exclusion.areas) != 2){
    stop("'exclusion.areas' needs to be a matrix with 2 columns.")
  }

  if(max.shift < 0){
    stop("'max.shift' must be non negative.")
  }

  if(threshold.noise < 0){
    stop("'threshold.noise' must be non negative.")
  }

  ##Seed and variables declaration
  set.seed(seed)


  #-----------------------------------------------------------------------------
  #### Import complexe mixture and pure spectra library ####
  import <- load_data(path = path, exclusion.areas = exclusion.areas,
                      max.shift = max.shift, which.spectra = which.spectra,
                      library.metabolites = library.metabolites)

  mixture <- import$mixture
  pure_library <- import$pure_library

  #number of points on library grid corresponding to maximum shift
  nb_points_shift <- floor(max.shift / (pure_library$grid[2] -
                                          pure_library$grid[1]))


  #-----------------------------------------------------------------------------
  #### Cleaning step: remove metabolites that cannot belong to the mixture ####
  pure_lib_clean <- clean_library(mixture, pure_library, threshold.noise,
                                  nb_points_shift)


  #-----------------------------------------------------------------------------
  #### Find the best translation between each pure spectra and mixture ####
  #and sort metabolites by regression residuals
  res_translation <- translate_library(mixture, pure_lib_clean, nb_points_shift)

  pure_lib_sorted <- res_translation$pure_lib_sorted
  shift <- res_translation$shift


  #-----------------------------------------------------------------------------
  #### Localized deformations of pure spectra ####
  pure_lib_deformed <- deform_library(mixture, pure_lib_sorted, nb_points_shift,
                                      max.shift, shift)


  #-----------------------------------------------------------------------------
  #### Threshold and concentration optimisation for each metabolites ####
  res_opti <- concentration_opti(mixture, pure_lib_deformed, nb.iter.signif)
  pure_lib_final <- res_opti$pure_lib_final


  #-----------------------------------------------------------------------------
  #### Results ####

  #Reconstituted mixture with estimated coefficents
  est_mixture <- as.numeric(pure_lib_final$spectra %*% res_opti$B_final)

  #Compute relative concentration of identified metabolites
  relative_concentration <- res_opti$B_final / pure_lib_final$nb_protons
  #Sort library according to relative concentration
  sorted_idx <- sort(relative_concentration, decreasing = TRUE,
                     index.return = TRUE)$ix
  pure_lib_final_sorted <- subset_library(pure_lib_final, sorted_idx)

  present_metab <- data.frame(Name = pure_lib_final_sorted$name,
                              Relative_Concentration =
                                relative_concentration[sorted_idx])


  #List to return
  res_object <- new(Class = "resASICS",
                    original_mixture = mixture,
                    reconstituted_mixture = est_mixture,
                    ppm_grid = pure_lib_final$grid,
                    present_metabolites = present_metab)

  return(res_object)
}

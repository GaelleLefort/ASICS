#' Automatic Statistical Identification in Complex Spectra
#'
#' Description
#' @param path folder path of the Bruker files
#' @param exclusion.areas exclusion areas of the quantification
#' @param max.shift maximum chemical shift allowed (in ppm)
#' @param which.spectra if more than one spectra by sample, spectra to choose
#' (either "first", "last" or its number)
#' @param library.metabolites path of the library of standard if not the default
#' one
#' @param threshold.noise threshold for signal noise
#' @return A object of type resASICS
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' result <- ASICS(path = system.file("extdata", "example_spectra",
#'                                    "AG_faq_Beck01", package = "ASICS"),
#'                 exclusion.areas = matrix(c(4.5,5.1,5.5,6.5),
#'                                          ncol = 2, byrow = TRUE))
#' }

ASICS <- function(path, exclusion.areas = matrix(c(4.5, 5.1), ncol = 2,
                                                 nrow = 1),
                  max.shift = 0.02, which.spectra = "last",
                  library.metabolites = NULL, threshold.noise = 0.02){

  ## checking the validity of the parameters
  if(!dir.exists(path)){
    stop("Path of the Bruker file does'nt exist !")
  }

  if(!is.matrix(exclusion.areas) | ncol(exclusion.areas) != 2){
    stop("exclusion.areas need to be a matrix of 2 columns.")
  }

  if(max.shift < 0){
    stop("max.shift must be positive.")
  }

  if(threshold.noise < 0){
    stop("threshold.noise must be positive.")
  }

  ##Seed and variables declaration
  set.seed(12345)


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
  res_opti <- concentration_opti(mixture, pure_lib_deformed)
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


  #Compute identification threshold of non identified metabolites
  identification_threshold <-
    round((res_opti$threshold[!res_opti$identified_metab] *
             pure_lib_final_sorted$nb_protons[1]) /
            (pure_lib_deformed$nb_protons[!res_opti$identified_metab]),
          digits = 4)

  non_identified_metab <-
    data.frame(Name = pure_lib_deformed$name[!res_opti$identified_metab],
               Threshold = identification_threshold)


  #List to return
  res_object <- new(Class = "resASICS",
                    original_mixture = mixture,
                    reconstituted_mixture = est_mixture,
                    grid = pure_lib_final$grid,
                    present_metabolites = present_metab,
                    non_identified_metabolites = non_identified_metab)

  return(res_object)
}

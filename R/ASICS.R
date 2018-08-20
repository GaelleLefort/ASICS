#' Automatic Statistical Identification in Complex Spectra
#'
#' Quantification of 1D 1H NMR spectra with ASICS method using a library of
#' pure metabolite spectra. The method is presented in Tardivel et al. (2017).
#'
#' @param spectra_obj An object of class \linkS4class{Spectra} obtained with the
#' function \link{createSpectra}.
#' @param exclusion.areas Definition domain of spectra that has to be excluded
#' for the quantification (ppm). By default, the water region is excluded
#' (4.5-5.1 ppm).
#' @param max.shift Maximum chemical shift allowed (in ppm). Default to 0.02.
#' @param pure.library An object of class \linkS4class{PureLibrary} containing
#' the reference spectra (pure metabolite spectra). If \code{NULL}, the library
#' included in the package (that contains 191 reference spectra) is used.
#' @param threshold.noise Threshold for signal noise. Default to 0.02.
#' @param seed Random seed to control randomness in the algorithm (used in the
#' estimation of the significativity of a given metabolite concentration).
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#'
#' @return An object of type \linkS4class{ASICSResults} containing the
#' quantification results.
#'
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SnowParam
#' @importFrom BiocParallel snowWorkers
#' @importFrom stats reshape
#' @export
#'
#' @seealso \linkS4class{ASICSResults} \code{\link{pure_library}}
#' \code{\link{createSpectra}}
#'
#' @references Tardivel P., Canlet C., Lefort G., Tremblay-Franco M., Debrauwer
#' L., Concordet D., Servien R. (2017). ASICS: an automatic method for
#' identification and quantification of metabolites in complex 1D 1H NMR
#' spectra. \emph{Metabolomics}, \strong{13}(10): 109.
#' \url{https://doi.org/10.1007/s11306-017-1244-5}
#'
#' @examples
#' # Import data and create object
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_obj <- createSpectra(spectra_data)
#'
#' # Estimation of relative quantification of Lactate and L-Alanine
#' to_exclude <- matrix(c(4.5, 10), ncol = 2)
#' pure_lib <- pure_library[getSampleName(pure_library) %in%
#'                          c("Lactate", "L-Alanine")]
#' resASICS <- ASICS(spectra_obj[1], exclusion.areas = to_exclude,
#'                   pure.library = pure_lib)
ASICS <- function(spectra_obj,
                  exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                  max.shift = 0.02, pure.library = NULL,
                  threshold.noise = 0.02, seed = 1234,
                  ncores = 1) {

  if(!is.null(exclusion.areas) &&
     (!is.matrix(exclusion.areas) | ncol(exclusion.areas) != 2)){
    stop("'exclusion.areas' must be a matrix with 2 columns.")
  }

  if(max.shift < 0){
    stop("'max.shift' must be non negative.")
  }

  if(threshold.noise < 0){
    stop("'threshold.noise' must be non negative.")
  }

  if(class(pure.library) != "PureLibrary" & !is.null(pure.library)){
    stop(paste("'pure.library' must be either NULL or an object of class",
               "'PureLibrary'."))
  }

  # number of cores
  ncores <- min(ncores, length(spectra_obj))
  if (.Platform$OS.type == "windows") {
    para_param <- SnowParam(workers = ncores,
                            progressbar = TRUE,
                            tasks = length(spectra_obj))
  } else {
    para_param <- MulticoreParam(workers = ncores,
                                 progressbar = TRUE,
                                 tasks = length(spectra_obj))
  }

  list_spec <- lapply(seq_along(spectra_obj), function(x) spectra_obj[x])

  res_estimation_list <- bplapply(list_spec,
                                  .ASICSInternal, exclusion.areas, max.shift,
                                  pure.library, threshold.noise, seed,
                                  BPPARAM = para_param)

  res_estimation <- do.call(c, res_estimation_list)

  return(res_estimation)

}


#' @importFrom methods new
.ASICSInternal <- function(spectrum_obj,
                           exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                           max.shift = 0.02, pure.library = NULL,
                           threshold.noise = 0.02, seed = 1234){

  # seed
  set.seed(seed)

  #-----------------------------------------------------------------------------
  #### Remove areas from spectrum and library ####
  cleaned_obj <- .removeAreas(spectrum_obj, exclusion.areas = exclusion.areas,
                              pure.library = pure.library)

  cleaned_spectrum <- cleaned_obj$cleaned_spectrum
  cleaned_library <- cleaned_obj$cleaned_library

  # number of points on library grid corresponding to maximum shift
  nb_points_shift <- floor(max.shift / (cleaned_library@ppm.grid[2] -
                                          cleaned_library@ppm.grid[1]))

  #-----------------------------------------------------------------------------
  #### Cleaning step: remove metabolites that cannot belong to the mixture ####
  cleaned_library <- .cleanLibrary(cleaned_spectrum, cleaned_library,
                                   threshold.noise, nb_points_shift)

  #-----------------------------------------------------------------------------
  #### Find the best translation between each pure spectra and mixture ####
  #and sort metabolites by regression residuals

  # compute weights
  s1 <- 0.172 #standard deviation of multiplicative noise
  s2 <- 0.15 #standard deviation of additive noise
  noises <- abs(cleaned_spectrum@spectra) * s1 ^ 2 + s2 ^ 2
  mixture_weights <- as.numeric(1 / noises)

  res_translation <- .translateLibrary(cleaned_spectrum, cleaned_library,
                                       mixture_weights, nb_points_shift)

  sorted_library <- res_translation$sorted_library
  shift <- res_translation$shift

  #-----------------------------------------------------------------------------
  #### Localized deformations of pure spectra ####
  deformed_library <- .deformLibrary(cleaned_spectrum, sorted_library,
                                     mixture_weights, nb_points_shift,
                                     max.shift, shift)

  #-----------------------------------------------------------------------------
  #### Threshold and concentration optimisation for each metabolites ####
  res_opti <- .concentrationOpti(cleaned_spectrum, deformed_library,
                                 noises, mixture_weights)
  final_library <- res_opti$final_library

  #-----------------------------------------------------------------------------
  #### Results ####

  # reconstituted mixture with estimated coefficents
  est_mixture <- final_library@spectra %*% res_opti$B_final
  pure_lib_final_conc <- final_library
  pure_lib_final_conc@spectra <- final_library@spectra %*%
    diag(res_opti$B_final)

  # compute relative concentration of identified metabolites
  relative_concentration <- res_opti$B_final / final_library@nb.protons
  # sort library according to relative concentration
  sorted_idx <- sort(relative_concentration, decreasing = TRUE,
                     index.return = TRUE)$ix
  pure_lib_final_sorted <- pure_lib_final_conc[sorted_idx]

  present_metab <- data.frame(relative_concentration[sorted_idx])
  rownames(present_metab) <- pure_lib_final_sorted@sample.name
  colnames(present_metab) <- spectrum_obj@sample.name

  # change format of pure library
  temp_df <- as.data.frame(pure_lib_final_sorted@spectra)
  rownames(temp_df) <- pure_lib_final_sorted@ppm.grid
  colnames(temp_df) <- pure_lib_final_sorted@sample.name
  pure_lib_format <-
    reshape(temp_df, idvar = "ppm_grid", ids = row.names(temp_df),
            times = names(temp_df), timevar = "metabolite_name",
            varying = list(names(temp_df)), direction = "long",
            v.names = "intensity")
  rownames(pure_lib_format) <- NULL
  pure_lib_format <- pure_lib_format[pure_lib_format$intensity != 0, ]
  pure_lib_format <- cbind(sample = rep(spectrum_obj@sample.name,
                                        nrow(pure_lib_format)),
                           pure_lib_format)

  # Object to return
  res_object <- new(Class = "ASICSResults",
                    sample.name = cleaned_spectrum@sample.name,
                    ppm.grid = cleaned_spectrum@ppm.grid,
                    spectra = cleaned_spectrum@spectra,
                    reconstructed.spectra = est_mixture,
                    quantification = present_metab,
                    deformed.library = pure_lib_format)

  return(res_object)
}

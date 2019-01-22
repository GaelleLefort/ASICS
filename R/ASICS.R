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
#' @param combine Logical. If \code{TRUE}, information from all spectra are
#' taken into account to align individual library.
#' @param seed Random seed to control randomness in the algorithm (used in the
#' estimation of the significativity of a given metabolite concentration).
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#' @param verbose A boolean value to allow print out process information.
#'
#' @return An object of type \linkS4class{ASICSResults} containing the
#' quantification results.
#'
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SerialParam
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
#' current_path <- system.file("extdata", package = "ASICS")
#' spectra_data <- importSpectra(name.dir = current_path,
#'                      name.file = "spectra_example.txt", type.import = "txt")
#' spectra_obj <- createSpectra(spectra_data)
#'
#' # Estimation of relative quantification of Lactate and L-Alanine
#' to_exclude <- matrix(c(4.5, 10), ncol = 2)
#' pure_lib <- pure_library[getSampleName(pure_library) %in%
#'                          c("Lactate", "L-Alanine")]
#' resASICS <- ASICS(spectra_obj[1], exclusion.areas = to_exclude,
#'                   pure.library = pure_lib, combine = FALSE)
ASICS <- function(spectra_obj,
                  exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                  max.shift = 0.02, pure.library = NULL,
                  threshold.noise = 0.02, combine = TRUE, seed = 1234,
                  ncores = 1, verbose = TRUE) {

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

  if(!is(pure.library, "PureLibrary") & !is.null(pure.library)){
    stop(paste("'pure.library' must be either NULL or an object of class",
               "'PureLibrary'."))
  }

  res_estimation <- .ASICSInternal(spectra_obj, exclusion.areas, max.shift,
                                   pure.library, threshold.noise,  seed, ncores,
                                   combine, verbose)

  return(res_estimation)
}



#' @importFrom methods new
.ASICSInternal <- function(spectra_obj_raw,
                           exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                           max.shift = 0.02, pure.library = NULL,
                           threshold.noise = 0.02, seed = 1234, ncores = 1,
                           combine = TRUE, verbose = TRUE){

  # seed and parallel environment
  set.seed(seed)

  ncores <- min(ncores, length(spectra_obj_raw))
  if (.Platform$OS.type == "windows" | ncores == 1) {
    para_param <- SerialParam(progressbar = verbose)
  } else {
    para_param <- MulticoreParam(workers = ncores,
                                 progressbar = verbose,
                                 tasks = length(spectra_obj_raw),
                                 manager.hostname = "localhost")
  }

  # default library or not
  if(is.null(pure.library)){
    pure.library <- ASICS::pure_library
  }

  # spectra object as a list where each element are 1 spectrum
  spectra_list <- lapply(seq_along(spectra_obj_raw),
                         function(x) spectra_obj_raw[x])


  #-----------------------------------------------------------------------------
  #### Remove areas from spectrum and library ####
  if (verbose) cat("Remove areas from spectrum and library \n")
  spectra_obj <- bplapply(spectra_list, .removeAreas, exclusion.areas,
                          pure.library, BPPARAM = para_param)

  # number of points on library grid corresponding to maximum shift
  if (length(spectra_list) == 1 | !combine) {
    nb_points_shift <-
      floor(max.shift / (spectra_obj[[1]][["cleaned_library"]]@ppm.grid[2] -
                           spectra_obj[[1]][["cleaned_library"]]@ppm.grid[1]))
  } else {
    max.shift <- seq_len(5) * max.shift / 5
    nb_points_shift <-
      floor(max.shift / (spectra_obj[[1]][["cleaned_library"]]@ppm.grid[2] -
                           spectra_obj[[1]][["cleaned_library"]]@ppm.grid[1]))
  }

  #-----------------------------------------------------------------------------
  #### Cleaning step: remove metabolites that cannot belong to the mixture ####
  if (verbose) cat("Remove metabolites that cannot belong to the mixture \n")
  spectra_obj <- bplapply(spectra_obj, .cleanLibrary, threshold.noise,
                          nb_points_shift[length(nb_points_shift)],
                          BPPARAM = para_param)

  #-----------------------------------------------------------------------------
  #### Find the best translation between each pure spectra and mixture ####
  #and sort metabolites by regression residuals

  # compute weights
  s1 <- 0.172 #standard deviation of multiplicative noise
  s2 <- 0.15 #standard deviation of additive noise
  if (verbose) cat("Compute weights \n")
  spectra_obj <-
    bplapply(spectra_obj,
             function(x){x[["mixture_weights"]] <-
               as.numeric(1 / (abs(x[["cleaned_spectrum"]]@spectra) *
                                 s1 ^ 2 + s2 ^ 2)); return(x)},
             BPPARAM = para_param)

  if (verbose) cat("Translate library \n")
  if (length(spectra_list) == 1 | !combine) {
    spectra_obj <- bplapply(spectra_obj, .translateLibrary,
                            nb_points_shift[length(nb_points_shift)],
                            max.shift[length(max.shift)],
                              BPPARAM = para_param)
  } else {
    # spectra binning
    spec_bin <- binning(data.frame(as.matrix(getSpectra(spectra_obj_raw))),
                        exclusion.areas = exclusion.areas,
                        ncores = ncores, verbose = FALSE)
    spec_bin <- spec_bin[rowSums(spec_bin) != 0, ]
    spectra_obj <- .translateLibrary_combineVersion(spectra_obj, max.shift,
                                                 nb_points_shift, spec_bin,
                                                 pure.library, para_param,
                                                 verbose)
  }

  #-----------------------------------------------------------------------------
  #### Localized deformations of pure spectra ####
  if (verbose) cat("Deform library peaks \n")
  spectra_obj <- bplapply(spectra_obj, .deformLibrary, BPPARAM = para_param)

  #-----------------------------------------------------------------------------
  #### Threshold and concentration optimisation for each metabolites ####
  if (verbose) cat("Compute quantifications \n")
  spectra_obj <- bplapply(spectra_obj, .concentrationOpti, BPPARAM = para_param)

  #-----------------------------------------------------------------------------
  #### Results ####
  if (verbose) cat("Format results... \n")
  sample_name <-
    unlist(vapply(spectra_obj,
                  function(x) return(x[["cleaned_spectrum"]]@sample.name),
                  "character"))
  spectra <-
    do.call("cbind",
            lapply(spectra_obj,
                   function(x) return(x[["cleaned_spectrum"]]@spectra)))
  rec_spectra <-
    do.call("cbind", lapply(spectra_obj,
                            function(x) return(x[["est_mixture"]])))

  rel_conc <- lapply(spectra_obj, function(x) {
    x[["relative_concentration"]]$row_names <-
      x[["cleaned_library"]]@sample.name ;
    return(x[["relative_concentration"]])})
  metab_conc <- join_all(rel_conc, by = "row_names", type = "full")
  rownames(metab_conc) <- metab_conc$row_names
  metab_conc$row_names <- NULL
  metab_conc[is.na(metab_conc)] <- 0
  pure_lib_format <- do.call("rbind",
                             lapply(spectra_obj,
                                    function(x) return(x[["format_library"]])))
  # Object to return
  res_object <- new(Class = "ASICSResults",
                    sample.name = sample_name,
                    ppm.grid = spectra_obj[[1]][["cleaned_spectrum"]]@ppm.grid,
                    spectra = spectra,
                    reconstructed.spectra = rec_spectra,
                    quantification = metab_conc,
                    deformed.library = pure_lib_format)

  return(res_object)
}

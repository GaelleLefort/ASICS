#' Automatic Statistical Identification in Complex Spectra for many files
#'
#' Compute ASICS on multiple spectra of a folder.
#'
#' @param name.dir folder path of the Bruker files
#' @param exclusion.areas areas to exclude (in ppm)
#' @param max.shift maximum chemical shift allowed (in ppm)
#' @param which.spectra if more than one spectra by sample, spectra to choose
#' (either "first", "last" or its number)
#' @param library.metabolites path of the library containing the standard of pure metabolites, if not the default
#' one
#' @param threshold.noise threshold for signal noise
#' @param ncores number of cores to use
#' @return A list containing ASICS results for each spectrum.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' @seealso \link{ASICS}, \link{resASICS-class}
#' @examples
#' \dontrun{
#' res_multi <- ASICS_multiFiles(name.dir = system.file("extdata",
#'                                                      "example_spectra",
#'                                                      package = "ASICS"),
#'                               exclusion.areas = matrix(c(4.5,5.1,5.5,6.5),
#'                                                     ncol = 2, byrow = TRUE),
#'                               ncores = 2)
#' }

ASICS_multiFiles <- function(name.dir,
                             exclusion.areas = matrix(c(4.5, 5.1), ncol = 2,
                                                            nrow = 1),
                             max.shift = 0.02, which.spectra = "last",
                             library.metabolites = NULL, threshold.noise = 0.02,
                             ncores = 1){

  ## checking the validity of the parameters
  if(!dir.exists(name.dir)){
    stop("Path of the Bruker file does'nt exist !")
  }

  dossiers <- dir(name.dir)

  # Parallel environment
  cl <- parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)

  # Progress bar
  pb <- txtProgressBar(max = length(dossiers), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # ASICS for each spectrum
  `%dopar%` <- foreach::`%dopar%`
  i <- NULL
  res_estimation <- foreach::foreach(i = 1:length(dossiers),
                                     .options.snow = opts,
                                     .packages = "ASICS") %dopar% {
     path <- file.path(name.dir, dossiers[i])
     res <- NULL
     res <- try(ASICS(path = path, exclusion.areas = exclusion.areas,
                      max.shift = max.shift, which.spectra = which.spectra,
                      library.metabolites = library.metabolites,
                      threshold.noise = threshold.noise))
     res
  }

  # Close parallel environment
  close(pb)
  parallel::stopCluster(cl)

  names(res_estimation) <- dossiers
  return(res_estimation)
}


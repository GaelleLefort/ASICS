#' Automatic Statistical Identification in Complex Spectra for many files
#'
#' Compute ASICS on multiple spectra of a folder.
#'
#' @param name.dir folder path of the Bruker files
#' @param exclusion.areas definition domain of spectra to exclude (ppm)
#' @param max.shift maximum chemical shift allowed (in ppm)
#' @param which.spectra if more than one spectra by sample, spectra used to 
#' perform the quantification (either \code{"first"}, \code{"last"} or its 
#' number). Default to \code{"last"}
#' @param library.metabolites path of the library containing the references 
#' (pure metabolite spectra). If \code{NULL}, the library included in the 
#' package is used
#' @param threshold.noise threshold for signal noise
#' @param ncores number of cores to use. Default to 1 (no parallel processing)
#' 
#' @return A list containing \code{ASICS} results for each spectrum.
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
#' @export
#' 
#' @seealso \code{\link{ASICS}}, \code{\link{resASICS-class}} 
#' \code{\link{pure_library}}
#' 
#' @examples
#' \dontrun{
#' cur_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' res_multi <- ASICS_multiFiles(name.dir = cur_path, 
#'                               exclusion.areas = to_exclude,
#'                               ncores = 2)
#' }

ASICS_multiFiles <- function(name.dir,
                             exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                             max.shift = 0.02, which.spectra = "last",
                             library.metabolites = NULL, threshold.noise = 0.02,
                             ncores = 1){

  ## checking the validity of the parameters
  if(!dir.exists(name.dir)){
    stop("Path of the Bruker file doesn't exist!")
  }

  cur_dir <- dir(name.dir)

  # Parallel environment
  cl <- parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)

  # Progress bar
  pb <- txtProgressBar(max = length(cur_dir), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # ASICS for each spectrum
  `%dopar%` <- foreach::`%dopar%`
  i <- NULL
  res_estimation <- foreach::foreach(i = 1:length(cur_dir),
                                     .options.snow = opts,
                                     .packages = "ASICS") %dopar% {
     path <- file.path(name.dir, cur_dir[i])
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

  names(res_estimation) <- cur_dir
  return(res_estimation)
}


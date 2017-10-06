#' Automatic Statistical Identification in Complex Spectra for many files
#'
#' Compute ASICS on multiple spectra of a folder.
#'
#' @param name.dir folder path of the Bruker files
#' @param ncores number of cores to use. Default to 1 (no parallel processing)
#' @param ... further arguments to be passed to the function
#' \code{\link{ASICS}} for specifying the parameters of the algorithm
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
#' \dontshow{
#' lib_file <- system.file("extdata", "library_for_examples.rda",
#'                         package = "ASICS")
#' cur_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' result <- ASICS_multiFiles(name.dir = cur_path,
#'                            exclusion.areas = to_exclude,
#'                            nb.iter.signif = 10, which.spectra = 2,
#'                            library.metabolites = lib_file, ncores = 1)
#' }
#' \dontrun{
#' cur_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' res_multi <- ASICS_multiFiles(name.dir = cur_path,
#'                               exclusion.areas = to_exclude,
#'                               ncores = 1)
#' }

ASICS_multiFiles <- function(name.dir, ncores = 1, ...){
  ## ASICS parameters
  param.args <- list(...)

  ## checking the validity of the parameters
  if(!dir.exists(name.dir)){
    stop("Path of the Bruker file doesn't exist!")
  }

  cur_dir <- dir(name.dir)

  # Progress bar
  pb <- txtProgressBar(max = length(cur_dir), style = 3)

  # Parallel environment
  if(ncores > 1) {
    cl <- snow::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`

    # Progress bar
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    `%dopar%` <- foreach::`%do%`

    # Progress bar
    opts <- NULL
  }

  # ASICS for each spectrum
  i <- NULL
  res_estimation <- foreach::foreach(i = 1:length(cur_dir),
                                     .options.snow = opts,
                                     .packages = "ASICS") %dopar% {
     path <- file.path(name.dir, cur_dir[i])
     res <- NULL
     res <- suppressWarnings(try(do.call("ASICS", c(path, param.args)),
                                 silent = TRUE))

     if(ncores == 1){
       setTxtProgressBar(pb, i)
     }

     res
  }

  # Close parallel environment
  if(ncores != 1) {
    close(pb)
    snow::stopCluster(cl)
  }

  names(res_estimation) <- cur_dir
  return(res_estimation)
}


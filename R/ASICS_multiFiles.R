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
#' @importFrom BiocParallel bplapply MulticoreParam
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
#'                               exclusion.areas = to_exclude)
#' }

ASICS_multiFiles <- function(name.dir, ...){
  ## ASICS parameters
  param.args <- list(...)

  ## checking the validity of the parameters
  if(!dir.exists(name.dir)){
    stop("Path of the Bruker file doesn't exist!")
  }

  cur_dir <- as.list(file.path(name.dir, dir(name.dir)))

  res_estimation <- bplapply(cur_dir,
                             function(x) suppressWarnings(
                               try(do.call("ASICS", c(x, param.args)),
                                                              silent = TRUE)),
                             BPPARAM = MulticoreParam(progressbar = TRUE,
                                                      tasks = length(cur_dir)))

  names(res_estimation) <- dir(name.dir)
  return(res_estimation)
}


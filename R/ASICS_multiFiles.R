#' Automatic Statistical Identification in Complex Spectra for many files
#'
#' Description
#' @param name folder path of the Bruker files
#' @param ZoneAEnlever exclusion areas to remove before the quantification
#' @param DecalageMax maximum chemical shift allowed
#' @param DecalageLib chemical shift between library of pure spectra and complex mixture
#' @param sample if more than one spectra by sample, spectra to choose (either "first",
#' "last" or its number
#' @param seed seed
#' @param libraryMetab path of the library of standard if not the default one
#' @param seuilBruit threshold noise
#' @param ncores number of cores to use
#' @keywords NMR quantification metabolites
#' @export
#' @examples
#' ASICS_multiFiles(nameDir = "./spectres_exemple/",
#' ZoneAEnlever = matrix(c(4.5, 5.1), ncol = 2,
#' nrow = 1),
#' DecalageMax = 0.021, DecalageLib = -0.002,
#' ncores = 4, sample = "last", seed = 12345,
#' libraryMetab = NULL)
ASICS_multiFiles <- function(nameDir, ZoneAEnlever = matrix(c(4.5, 5.1), ncol = 2,
                                                            nrow = 1),
                             DecalageMax = 0.021, DecalageLib = 0,
                             ncores = 1, sample = "last", seed = 12345,
                             libraryMetab = NULL, seuilBruit = 0.02){

  dossiers <- dir(nameDir)

  #Parallel
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)

  ##Progress bar
  pb <- txtProgressBar(max = length(dossiers), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  resEstimation <- foreach(i=1:length(dossiers), .options.snow = opts,
                           .packages = "ASICS") %dopar% {
     name <- file.path(nameDir, dossiers[i])
     res <- NULL
     res <- try(ASICS(name = name, ZoneAEnlever = ZoneAEnlever, DecalageMax = DecalageMax, DecalageLib = DecalageLib,
                      sample = sample, seed = seed, libraryMetab = libraryMetab, seuilBruit = seuilBruit))
     res
                           }
  close(pb)
  stopCluster(cl)

  names(resEstimation) <- dossiers
  return(resEstimation)
}

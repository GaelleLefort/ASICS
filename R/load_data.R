#' Import from Bruker files
#'
#' Import spectra from Bruker files contained in a unique folder. This folder
#' contains subfolders for each sample.
#'
#' @param name.dir Path of the folder containing one subfolder by sample. Each
#' subfolder contains Bruker files.
#' @param which.spectra If more than one spectra by sample, spectra to import
#' (either \code{"first"}, \code{"last"} or its number). Default to
#' \code{"last"}
#' @param ppm.grid Numeric vector of a unique grid (definition domain) for all
#' spectra (in ppm). Default to \code{NULL} (grid of default pure library is
#' used).
#' @param sample.names Character vector of sample names. Default to \code{NULL}
#' (folder names are used).
#' @param parallel if \code{TRUE}, apply function in parallel. Default to
#' \code{TRUE}.
#'
#' @return A data frame with spectra in columns and chemical shifts (in p.p.m.)
#' in rows.
#'
#' @importFrom BiocParallel bplapply bptry MulticoreParam bpok
#' @importFrom utils head tail
#' @export
#'
#spectra_dat <- import_spectra_bruker(system.file("extdata", "example_spectra", package = "ASICS"))
import_spectra_bruker <- function(name.dir, which.spectra = "last",
                                  ppm.grid = NULL, sample.names = NULL,
                                  parallel = TRUE){

  if(!dir.exists(name.dir)){
    stop("Path of the Bruker files doesn't exist!")
  }

  if (length(which.spectra) != 1 | (!(which.spectra %in% c("first", "last")) &
      !is.numeric(which.spectra))) {
    stop(paste("'which.spectra' must be of length 1 and either 'first', 'last'",
               " or a number."))
  }

  if (!is.null(ppm.grid) & !is.numeric(ppm.grid)) {
    stop("'ppm.grid' must be NULL or a numeric vector.")
  }

  # path of all samples
  cur_dir <- as.list(file.path(name.dir, dir(name.dir)))

  if (!is.null(sample.names) & length(sample.names) != length(cur_dir)) {
    stop(paste0("'sample.names' must be NULL or have the same length than the ",
         "number of folder in ", name.dir, "."))
  }

  # path of the chosen spectra
  if (length(which.spectra) == 1 && which.spectra == "first") {
    cur_dir_spec <- as.list(file.path(cur_dir,
                                      sapply(as.list(sapply(cur_dir, dir)),
                                             head, 1)))
  } else if (length(which.spectra) == 1 && which.spectra == "last") {
    cur_dir_spec <- as.list(file.path(cur_dir,
                                      sapply(as.list(sapply(cur_dir, dir)),
                                             tail, 1)))
  } else {
    cur_dir_spec <- as.list(file.path(cur_dir, which.spectra))
  }

  # default ppm grid and sample names
  if (is.null(ppm.grid)) ppm.grid <- ASICS::pure_library@ppm.grid
  if (is.null(sample.names)) sample.names <- dir(name.dir)

  if(any(duplicated(sample.names))){
    stop("Sample names need to be unique.")
  }

  # number of cores
  if (parallel) {
    ncores <- multicoreWorkers()
  } else {
    ncores <- 1
  }

  #import spectra
  imported_spectra <-
    bptry(bplapply(cur_dir_spec, read_NMR_bruker, ppm.grid,
                   BPPARAM = MulticoreParam(workers = ncores,
                                            progressbar = TRUE,
                                            tasks = length(cur_dir_spec))))

  #if error in a file
  if (any(!bpok(imported_spectra))) {
    warning(paste0("There is a problem in files for spectra: ",
                paste(sample.names[!bpok(imported_spectra)], collapse = ", "),
                ".\nFix the problem manually and re-execute the function " ,
                "otherwise these spectra are ignored."), call. = FALSE)
    sample.names <- sample.names[bpok(imported_spectra)]
    imported_spectra <- imported_spectra[bpok(imported_spectra)]
  }

  #convert in a data frame
  imported_spectra <- as.data.frame(do.call(cbind, imported_spectra),
                                    row.names = as.character(ppm.grid))
  colnames(imported_spectra) <- sample.names

  return(imported_spectra)
}


## Extract a spectrum from Bruker files
read_NMR_bruker <- function(path, ppm.grid){

  #file paths
  acq_file <- file.path(path, "acqus")
  spec_file <- file.path(path, "pdata/1/1r")
  proc_file <- file.path(path, "pdata/1/procs")

  #meta-data in acq_file
  ACQ <- readLines(acq_file)
  TD <- get_bruker_param(ACQ, "TD")
  SW <- get_bruker_param(ACQ, "SW")
  SWH <- get_bruker_param(ACQ, "SW_h")
  DTYPA <- get_bruker_param(ACQ, "DTYPA")
  BYTORDA <- get_bruker_param(ACQ, "BYTORDA")
  ENDIAN <- "little"
  SIZE <- ifelse(DTYPA == 0, 4, 8)

  #meta-data in proc_file
  PROC <- readLines(proc_file)
  OFFSET <- get_bruker_param(PROC, "OFFSET")
  SI <- get_bruker_param(PROC, "SI")

  #meta-data in spec_file
  to.read <- file(spec_file, "rb")
  maxTDSI <- max(TD, SI)
  signal <- rev(readBin(to.read, what = "int", size = SIZE, n = maxTDSI,
                        signed = TRUE, endian = ENDIAN))
  close(to.read)

  #get signal
  td <- length(signal)
  dppm <- SW / (td - 1)
  pmax <- OFFSET
  pmin <- OFFSET - SW
  ppmseq <- seq(from = pmin, to = pmax, by = dppm)
  signal <- 100 * signal / max(signal)

  #adapt the grid
  new_signal <- change_grid(signal, ppmseq, ppm.grid)

  return(new_signal)
}

# find parameter in a bruker file
get_bruker_param <- function(file, paramStr){
  regexpStr <- paste0("^...", paramStr, "=")
  as.numeric(gsub("^[^=]+= ", "" ,
                  file[which(simplify2array(regexpr(regexpStr, file)) > 0)]))
}




#' Create a pure library
#'
#' Create a new pure library from a data frame containing spectra. Noise is
#' removed by thresholding each spectrum.
#'
#' @param spectra Data frame with spectra in column and chemical shift in row.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in p.p.m).
#' @param threshold Numeric value below which pure spectrum values are
#' considered null.
#' @param nb.protons Numeric vector of the number of protons of each pure
#' metabolite spectra contained in \code{spectra} data frame.
#'
#' @return A \linkS4class{PureLibrary} object with the newly created library.
#'
#' @importFrom methods new
#' @export
#'
# nb.protons c(5,4,9,9)
create_pure_library <- function(spectra, threshold = 5, nb.protons){
  # create a new object PureLibrary with a data frame of intensities and a
  # vector of the muber of protons
  new_library <- new("PureLibrary",
                     sample.name = colnames(spectra),
                     ppm.grid = as.numeric(rownames(spectra)),
                     spectra = as.matrix(spectra),
                     nb.protons = nb.protons)
  colnames(new_library@spectra) <- NULL
  rownames(new_library@spectra) <- NULL

  # spectra thresholding
  new_library@spectra[new_library@spectra < threshold] <- 0

  return(new_library)
}



#' Preprocessing on spectra
#'
#' Create a new pure spectra object used for quantification. New spectra are
#' baseline corrected and normalized by the area under the curve.
#' @name preprocessing_spectra
#' @param spectra Data frame with spectra in column and chemical shift in row.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in p.p.m).
#' @param parallel if \code{TRUE}, apply function in parallel. Default to
#' \code{TRUE}.
#'
#' @return A \linkS4class{Spectra} object with baseline corrected and normalized
#' spectra.
#'
#' @export
#' @importFrom BiocParallel bplapply multicoreWorkers
#'
# ref à publi de correction
# spec_obj <- preprocessing_spectra(spectra_dat)
preprocessing_spectra <- function(spectra, parallel = TRUE){

  # create a new object Spectra with a data frame of intensities
  new_spectra <- new("Spectra",
                     sample.name = colnames(spectra),
                     ppm.grid = as.numeric(rownames(spectra)),
                     spectra = as.matrix(spectra))

  # number of cores
  if (parallel) {
    ncores <- multicoreWorkers()
  } else {
    ncores <- 1
  }

  new_spectra@spectra <-
    do.call(cbind,
            bplapply(as.list(seq_along(new_spectra)), internal_preprocessing,
                     new_spectra@spectra, new_spectra@ppm.grid,
                     BPPARAM = MulticoreParam(workers = ncores,
                                              progressbar = TRUE,
                                              tasks = length(new_spectra))))
  rownames(new_spectra@spectra) <- NULL

  return(new_spectra)
}

## internal function for preprocessing
internal_preprocessing <- function(idx, spectra, ppm.grid) {
  # baseline correction
  new_spectra <- baseline_corrector(spectra[, idx])
  new_spectra[new_spectra < 0] <- 0

  # normalisation by area under the curve
  new_spectra <- new_spectra / AUC(ppm.grid, new_spectra)

  return(new_spectra)
}



## Main function to load library, import a spectrum from Bruker files and
#perform the first preprocessing
#' @importFrom plyr alply
remove_areas <- function(spectrum_obj, exclusion.areas, pure.library){

  #default library or not
  if(is.null(pure.library)){
    cleaned_library <- pure_library
  } else {
    cleaned_library <- pure.library
  }

  #adapt spectrum grid to have the same than library
  cleaned_spectrum <- new("Spectra",
                          sample.name = spectrum_obj@sample.name,
                          ppm.grid = cleaned_library@ppm.grid,
                          spectra = apply(spectrum_obj@spectra, 2,
                                          change_grid,
                                          spectrum_obj@ppm.grid,
                                          cleaned_library@ppm.grid))


  #for signal in exclusion.areas, intensity is null (mixture and library)
  idx_to_remove <-
    unlist(alply(exclusion.areas, 1,
                 function(x) which(cleaned_library@ppm.grid >= x[[1]] &
                                     cleaned_library@ppm.grid <= x[[2]])),
           use.names = FALSE)

  cleaned_spectrum@spectra[idx_to_remove, ] <- 0
  cleaned_library@spectra[idx_to_remove, ] <- 0

  #remove metabolite without any signal
  with_signal <- as.numeric(which(colSums(cleaned_library@spectra) > 0))
  cleaned_library <- cleaned_library[with_signal]

  #re-normalisation of mixture and pure spectra library
  cleaned_library@spectra <- apply(cleaned_library@spectra, 2,
                                   function(x)
                                     t(x / AUC(cleaned_library@ppm.grid, x)))
  cleaned_spectrum@spectra <- apply(cleaned_spectrum@spectra, 2,
                                    function(x)
                                      t(x / AUC(cleaned_library@ppm.grid, x)))

  return(list("cleaned_spectrum" = cleaned_spectrum,
              "cleaned_library" = cleaned_library))

}

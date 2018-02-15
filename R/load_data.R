#' Import from Bruker files
#'
#' Import spectra from Bruker files contained in a single folder. This folder
#' contains subfolders for each sample. Spectra are baseline corrected
#' (optionnal) and normalized by the area under the curve.
#'
#' @param name.dir Path of the folder containing one subfolder by sample. Each
#' subfolder contains Bruker files.
#' @param which.spectra If more than one spectra by sample, spectra to import
#' (either \code{"first"}, \code{"last"} or its number). Default to
#' \code{"last"}.
#' @param baseline.correction If \code{TRUE} apply a baseline correction for
#' each spectrum (Wang et al (2013)). Default to \code{TRUE}.
#' @param ppm.grid Numeric vector of a unique grid (definition domain) for all
#' spectra (in ppm). Default to \code{NULL} (grid of default pure library is
#' used).
#' @param sample.names Character vector of sample names. Default to \code{NULL}
#' (folder names are used).
#' @param parallel If \code{TRUE}, apply function in parallel. Default to
#' \code{TRUE}.
#'
#' @references Wang, K.C., Wang, S.Y., Kuo, C.H., Tseng Y.J. (2013).
#' Distribution-based classification method for baseline correction of
#' metabolomic 1D proton nuclear magnetic resonance spectra.
#' \emph{Analytical Chemistry}, \strong{85}(2), 1231-1239.
#'
#' @return A data frame with spectra in columns and chemical shifts (in p.p.m.)
#' in rows.
#'
#' @importFrom BiocParallel bplapply bptry MulticoreParam bpok
#' @importFrom plyr llply
#' @importFrom utils head tail
#' @export
#'
#' @seealso \code{\link{baseline_correction}} \code{\link{normalisation}}
#' \code{\link{alignment}}
#'
#' @examples
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- import_spectra_bruker(current_path)
import_spectra_bruker <- function(name.dir, which.spectra = "last",
                                  baseline.correction = TRUE,
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
                                      sapply(as.list(lapply(cur_dir, function(x)
                                        sort(as.numeric(dir(x))))),
                                             head, 1)))
  } else if (length(which.spectra) == 1 && which.spectra == "last") {
    cur_dir_spec <- as.list(file.path(cur_dir,
                                      sapply(as.list(lapply(cur_dir, function(x)
                                        sort(as.numeric(dir(x))))),
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
  if (ncores > 1) {
    imported_spectra <-
      bptry(bplapply(cur_dir_spec, read_NMR_bruker, ppm.grid,
                     baseline.correction,
                     BPPARAM = MulticoreParam(workers = ncores,
                                              progressbar = TRUE,
                                              stop.on.error = FALSE,
                                              tasks = length(cur_dir_spec))))
  } else {
    suppressWarnings(
      imported_spectra <- llply(cur_dir_spec,
                                function(x, ppm, bc)
                                  try(read_NMR_bruker(x, ppm, bc), silent=TRUE),
                                ppm.grid, baseline.correction,
                                .progress = "text"))
  }

  #if error in a file
  if (any(!bpok(imported_spectra)) |
      any(sapply(imported_spectra, is, "try-error"))) {
    warning(paste0("There is a problem in files for spectra: ",
                paste(sample.names[!bpok(imported_spectra) |
                                     sapply(imported_spectra, is, "try-error")],
                      collapse = ", "),
                ".\nFix the problem manually and re-execute the function " ,
                "otherwise these spectra are ignored."), call. = FALSE)
    sample.names <- sample.names[bpok(imported_spectra) &
                                   !sapply(imported_spectra, is, "try-error")]
    imported_spectra <- imported_spectra[bpok(imported_spectra) &
                                           !sapply(imported_spectra, is,
                                                   "try-error")]
  }

  #convert in a data frame
  imported_spectra <- as.data.frame(do.call(cbind, imported_spectra),
                                    row.names = as.character(ppm.grid))
  colnames(imported_spectra) <- sample.names

  return(imported_spectra)
}


## Extract a spectrum from Bruker files
read_NMR_bruker <- function(path, ppm.grid, baseline.correction){

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

  # baseline correction
  if (baseline.correction) {
    signal <- baseline_corrector(signal)
    signal[signal < 0] <- 0
  }

  #adapt the grid
  new_signal <- change_grid(signal, ppmseq, ppm.grid)

  # normalisation by area under the curve
  new_signal <- new_signal / AUC(ppm.grid, new_signal)

  return(new_signal)
}

# find parameter in a bruker file
get_bruker_param <- function(file, paramStr){
  regexpStr <- paste0("^...", paramStr, "=")
  as.numeric(gsub("^[^=]+= ", "" ,
                  file[which(simplify2array(regexpr(regexpStr, file)) > 0)]))
}



#' Normalisation
#'
#' Normalise a data frame of spectra by the area under the curve.
#'
#' @param spectra Data frame with spectra in columns and chemical shift in rows.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in p.p.m).
#'
#' @return A data frame with normalised spectra in columns and chemical shifts
#' (in p.p.m.) in rows.
#'
#' @export
#'
#' @examples
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_norm <- normalisation(spectra_data)
normalisation <- function(spectra){

  if (!is.numeric(as.numeric(rownames(spectra)))) {
    stop("Rownames of spectra data frame don't contain ppm grid.")
  }

  spectra_norm <- data.frame(apply(spectra, 2, function(x)
    t(x / AUC(as.numeric(rownames(spectra)), x))))

  return(spectra_norm)
}



#' Baseline correction
#'
#' Apply a baseline correction to a data frame of spectra with the algorithm
#' described in Wang et al. (2013).
#'
#' @param spectra Data frame with spectra in columns and chemical shift in rows.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in p.p.m).
#' @param parallel If \code{TRUE}, apply function in parallel. Default to
#' \code{TRUE}.
#'
#' @return A data frame with baseline corrected spectra in columns and chemical
#' shifts (in p.p.m.) in rows.
#'
#' @references Wang, K.C., Wang, S.Y., Kuo, C.H., Tseng Y.J. (2013).
#' Distribution-based classification method for baseline correction of
#' metabolomic 1D proton nuclear magnetic resonance spectra.
#' \emph{Analytical Chemistry}, \strong{85}(2), 1231-1239.
#'
#' @export
#'
#' @examples
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_base_cor <- baseline_correction(spectra_data)
baseline_correction <- function(spectra, parallel = TRUE){

  if (!is.numeric(as.numeric(rownames(spectra)))) {
    stop("Rownames of spectra data frame don't contain ppm grid.")
  }

  # number of cores
  if (parallel) {
    ncores <- multicoreWorkers()
  } else {
    ncores <- 1
  }

  spectra_list <- as.list(spectra)

  #baseline correction
  if (ncores > 1) {
    spectra_bc_list <-
      bptry(bplapply(spectra_list, baseline_corrector,
                     BPPARAM = MulticoreParam(workers = ncores,
                                              progressbar = TRUE,
                                              stop.on.error = FALSE,
                                              tasks = ncol(spectra))))
  } else {
    suppressWarnings(
      spectra_bc_list <- llply(spectra_list,
                                function(x)
                                  try(baseline_corrector(x), silent=TRUE),
                                .progress = "text"))
  }

  #if error in a file
  if (any(!bpok(spectra_bc_list)) |
      any(sapply(spectra_bc_list, is, "try-error"))) {
    warning(paste0("Baseline correction algorithm can't be used for spectra: ",
                   paste(colnames(spectra)[!bpok(spectra_bc_list) |
                                             sapply(spectra_bc_list, is,
                                                    "try-error")],
                         collapse = ", ")), call. = FALSE)
  }
  spectra_list[!(!bpok(spectra_bc_list) | sapply(spectra_bc_list, is,
                                               "try-error"))] <-
    spectra_bc_list[!(!bpok(spectra_bc_list) | sapply(spectra_bc_list, is,
                                                      "try-error"))]

  #convert in a data frame
  spectra_bc <- as.data.frame(do.call(cbind, spectra_list),
                                    row.names = as.character(spectra))
  colnames(spectra_bc) <- colnames(spectra)

  spectra_bc[spectra_bc < 0] <- 0

  spectra_norm <- data.frame(apply(spectra_bc, 2, function(x)
    t(x / AUC(as.numeric(rownames(spectra)), x))))

  return(spectra_norm)
}



#' Alignment
#'
#' Align spectra of a data frame with the CluPA algorithm (Vu et al., (2011))
#'
#' @param spectra Data frame with spectra in columns and chemical shift in rows.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in p.p.m).
#' @param baseline.threshold Value of baseline threshold used to identify peaks.
#' Default to 0.02.
#' @param reference Index of reference spectra. Default to \code{NULL},
#' \emph{i.e.} reference spectrum is detect automatically with the
#' \code{findRef} function of package \code{speaq}
#' @param max.shift Maximum shift allowed for the alignment. Default to 0.002.
#'
#' @return A data frame with aligned spectra in columns and chemical shifts
#' (in p.p.m.) in rows.
#'
#' @references Vu, T. N., Valkenborg, D., Smets, K., Verwaest, K. A., Dommisse,
#' R., Lemiere, F., ... & Laukens, K. (2011). An integrated workflow for robust
#' alignment and simplified quantitative analysis of NMR spectrometry data.
#' \emph{BMC bioinformatics}, \strong{12}(1), 405.
#'
#' @export
#' @importFrom speaq findRef detectSpecPeaks dohCluster
#'
#' @examples
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_base_cor <- alignment(spectra_data)
alignment <- function(spectra, baseline.threshold = 0.02, reference = NULL,
                      max.shift = 0.002){

  if (!is.numeric(as.numeric(rownames(spectra)))) {
    stop("Rownames of spectra data frame don't contain ppm grid.")
  }

  # Maximum shift
  max_shift <- floor(max.shift / (ASICS::pure_library@ppm.grid[2] -
                                    ASICS::pure_library@ppm.grid[1]))

  # Detect peaks
  peak_list <- detectSpecPeaks(t(spectra),
                               nDivRange = 64,
                               scales = seq(1, 16, 2),
                               baselineThresh = baseline.threshold,
                               SNR.Th = -1,
                               verbose = FALSE)

  # Reference spectrum
  if (is.null(reference)) {
    reference <- findRef(peak_list)$refInd
  }

  spectra_align <- dohCluster(t(spectra),
                              peakList = peak_list,
                              refInd = reference,
                              maxShift  = max_shift,
                              acceptLostPeak = TRUE,
                              verbose = FALSE)
  spectra_align <- data.frame(t(spectra_align))

  spectra_norm <- data.frame(apply(spectra_align, 2, function(x)
    t(x / AUC(as.numeric(rownames(spectra_align)), x))))

  return(spectra_norm)
}

#' Create a pure library
#'
#' Create a new pure library from a data frame containing spectra. Noise is
#' removed by thresholding each spectrum.
#'
#' @param spectra Data frame with spectra in column and chemical shift in row.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in p.p.m).
#' @param nb.protons Numeric vector of the number of protons of each pure
#' metabolite spectra contained in \code{spectra} data frame.
#' @param threshold Numeric value below which pure spectrum values are
#' considered null. Default to 1.
#' @return A \linkS4class{PureLibrary} object with the newly created library.
#'
#' @importFrom methods new
#' @export
#'
#' @examples
#' pure_spectra <- import_spectra_bruker(system.file("extdata",
#'                                                   "example_library",
#'                                                    package = "ASICS"))
#' new_pure_library <- create_pure_library(pure_spectra,
#'                                         nb.protons = c(5, 4, 9, 9))
#'
create_pure_library <- function(spectra, nb.protons, threshold = 1){
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



#' Create a \linkS4class{Spectra} object
#'
#' Create a new spectra object used for quantification.
#'
#' @param spectra Data frame with spectra in column and chemical shift in row.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in p.p.m).
#'
#' @return A \linkS4class{Spectra} object with spectra to quantify.
#'
#' @export
#' @importFrom BiocParallel bplapply multicoreWorkers
#'
#' @seealso \linkS4class{Spectra}
#'
#' @examples
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- import_spectra_bruker(current_path)
#' spectra_obj <- create_spectra(spectra_data)
create_spectra <- function(spectra){

  # create a new object Spectra with a data frame of intensities
  new_spectra <- new("Spectra",
                     sample.name = colnames(spectra),
                     ppm.grid = as.numeric(rownames(spectra)),
                     spectra = as.matrix(spectra))

  rownames(new_spectra@spectra) <- NULL
  colnames(new_spectra@spectra) <- NULL

  return(new_spectra)
}


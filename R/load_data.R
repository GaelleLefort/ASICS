#' Import metabolomic spectra from Bruker files
#'
#' Import spectra from Bruker files contained in a single folder. This folder
#' contains one subfolder for each sample. Spectra are baseline corrected
#' (optional) and normalised by the area under the curve during the importation.
#'
#' @param name.dir Path of the folder containing one subfolder by sample. Each
#' subfolder contains the Bruker files of this sample.
#' @param which.spectra If there is more than one spectrum by sample, which is
#' the spectrum to import (either always the first one with \code{which.spectra
#' = "first"}, always the last one with \code{which.spectra = "last"} or a
#' vector of length the number of spectra that specifies the number of each
#' spectrum to import).
#' Default to \code{"last"}.
#' @param baseline.correction Logical. If \code{TRUE} a baseline correction is
#' applied for each spectrum (Wang et al (2013)). Default to \code{TRUE}.
#' @param alignment Logical. If \code{TRUE} a peak alignment is
#' applied between all spectra (Vu et al (2011)). Default to \code{FALSE}.
#' @param ppm.grid Numeric vector of a unique grid (definition domain) for all
#' spectra (in ppm). Default to \code{NULL} (in which case, the default grid of
#' the pure library is used).
#' @param sample.names Character vector of sample names. Default to \code{NULL}
#' (in which case, folder names are used).
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#' @param ... Further arguments to be passed to the function
#' \code{\link{alignment}} for specifying the parameters of the algorithm, if
#' necessary.
#'
#' @details
#' Some preprocessing steps are included during the importation. First, spectra
#' are baseline corrected if \code{baseline.correction = TRUE}. Then, all
#' spectrum definition domains are aligned to a unique one (either the one
#' specified in \code{ppm.grid} or the grid of the default library). Finally,
#' all spectra are normalised by the area under the curve.
#'
#' @references Wang, K.C., Wang, S.Y., Kuo, C.H., Tseng Y.J. (2013).
#' Distribution-based classification method for baseline correction of
#' metabolomic 1D proton nuclear magnetic resonance spectra.
#' \emph{Analytical Chemistry}, \strong{85}(2), 1231-1239.
#'
#' @references Vu, T. N., Valkenborg, D., Smets, K., Verwaest, K. A., Dommisse,
#' R., Lemiere, F., ... & Laukens, K. (2011). An integrated workflow for robust
#' alignment and simplified quantitative analysis of NMR spectrometry data.
#' \emph{BMC Bioinformatics}, \strong{12}(1), 405.
#'
#' @return A data frame with spectra in columns and chemical shifts (in ppm)
#' in rows.
#'
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SnowParam
#' @importFrom BiocParallel snowWorkers bptry bpok
#' @importFrom plyr llply
#' @importFrom utils head tail
#' @export
#'
#' @seealso \code{\link{baselineCorrection}} \code{\link{normalisation}}
#' \code{\link{alignment}}
#'
#' @examples
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- importSpectraBruker(current_path)
importSpectraBruker <- function(name.dir, which.spectra = "last",
                                baseline.correction = TRUE, alignment = FALSE,
                                ppm.grid = NULL, sample.names = NULL,
                                ncores = 1, ...){

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
                                      vapply(as.list(lapply(cur_dir, function(x)
                                        sort(as.numeric(dir(x))))),
                                        head, 1,
                                        FUN.VALUE = numeric(1))))
  } else if (length(which.spectra) == 1 && which.spectra == "last") {
    cur_dir_spec <- as.list(file.path(cur_dir,
                                      vapply(as.list(lapply(cur_dir, function(x)
                                        sort(as.numeric(dir(x))))),
                                        tail, 1,
                                        FUN.VALUE = numeric(1))))
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
  ncores <- min(ncores, length(cur_dir_spec))
  if (.Platform$OS.type == "windows") {
    para_param <- SnowParam(workers = ncores,
                            progressbar = TRUE,
                            tasks = length(cur_dir_spec),
                            stop.on.error = FALSE)
  } else {
    para_param <- MulticoreParam(workers = ncores,
                                 progressbar = TRUE,
                                 tasks = length(cur_dir_spec),
                                 stop.on.error = FALSE)
  }


  # import spectra
  if (ncores > 1) {
    imported_spectra <-
      bptry(bplapply(cur_dir_spec, .readNMRBruker, ppm.grid,
                     baseline.correction,
                     BPPARAM = para_param))
  } else {
    suppressWarnings(
      imported_spectra <- llply(cur_dir_spec,
                                function(x, ppm, bc)
                                  try(.readNMRBruker(x, ppm, bc), silent=TRUE),
                                ppm.grid, baseline.correction,
                                .progress = "text"))
  }

  # if error in a file
  if (any(!bpok(imported_spectra)) |
      any(vapply(imported_spectra, is, "try-error", FUN.VALUE = logical(1)))) {
    warning(paste0("There is a problem in files for spectra: ",
                   paste(sample.names[!bpok(imported_spectra) |
                                        vapply(imported_spectra, is, "try-error",
                                               FUN.VALUE = logical(1))],
                         collapse = ", "),
                   ".\nFix the problem manually and re-execute the function " ,
                   "otherwise these spectra are ignored."), call. = FALSE)
    sample.names <- sample.names[bpok(imported_spectra) &
                                   !vapply(imported_spectra, is, "try-error",
                                           FUN.VALUE = logical(1))]
    imported_spectra <- imported_spectra[bpok(imported_spectra) &
                                           !vapply(imported_spectra, is,
                                                   "try-error",
                                                   FUN.VALUE = logical(1))]
  }

  # convert in a data frame
  imported_spectra <- as.data.frame(do.call(cbind, imported_spectra),
                                    row.names = as.character(ppm.grid))
  colnames(imported_spectra) <- sample.names

  # peak alignment
  if (alignment) {
    imported_spectra <- alignment(imported_spectra, ...)
  }

  return(imported_spectra)
}


## Extract a spectrum from Bruker files
.readNMRBruker <- function(path, ppm.grid, baseline.correction){

  # file paths
  acq_file <- file.path(path, "acqus")
  spec_file <- file.path(path, "pdata/1/1r")
  proc_file <- file.path(path, "pdata/1/procs")

  # meta-data in acq_file
  ACQ <- readLines(acq_file)
  TD <- .getBrukerParam(ACQ, "TD")
  SW <- .getBrukerParam(ACQ, "SW")
  SWH <- .getBrukerParam(ACQ, "SW_h")
  DTYPA <- .getBrukerParam(ACQ, "DTYPA")
  BYTORDA <- .getBrukerParam(ACQ, "BYTORDA")
  ENDIAN <- "little"
  SIZE <- ifelse(DTYPA == 0, 4, 8)

  # meta-data in proc_file
  PROC <- readLines(proc_file)
  OFFSET <- .getBrukerParam(PROC, "OFFSET")
  SI <- .getBrukerParam(PROC, "SI")

  # meta-data in spec_file
  to.read <- file(spec_file, "rb")
  maxTDSI <- max(TD, SI)
  signal <- rev(readBin(to.read, what = "int", size = SIZE, n = maxTDSI,
                        signed = TRUE, endian = ENDIAN))
  close(to.read)

  # get signal
  td <- length(signal)
  dppm <- SW / (td - 1)
  pmax <- OFFSET
  pmin <- OFFSET - SW
  ppmseq <- seq(from = pmin, to = pmax, by = dppm)
  signal <- 100 * signal / max(signal)

  # baseline correction
  if (baseline.correction) {
    signal <- .baselineCorrector(signal)
    signal[signal < 0] <- 0
  }

  # adapt the grid
  new_signal <- .changeGrid(signal, ppmseq, ppm.grid)

  # normalisation by area under the curve
  new_signal <- new_signal / .AUC(ppm.grid, new_signal)

  return(new_signal)
}

## Find parameter in a bruker file
.getBrukerParam <- function(file, paramStr){
  regexpStr <- paste0("^...", paramStr, "=")
  as.numeric(gsub("^[^=]+= ", "" ,
                  file[which(simplify2array(regexpr(regexpStr, file)) > 0)]))
}



#' Normalisation
#'
#' Normalise a data frame of spectra by the area under the curve.
#'
#' @param spectra Data frame with spectra in columns and chemical shifts in
#' rows. Colnames of this data frame correspond to pure metabolite names and
#' rownames to chemical shift grid (in ppm).
#'
#' @return A data frame with normalised spectra in columns and chemical shifts
#' (in ppm) in rows.
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
    t(x / .AUC(as.numeric(rownames(spectra)), x))))
  rownames(spectra_norm) <- rownames(spectra)
  colnames(spectra_norm) <- colnames(spectra)

  return(spectra_norm)
}



#' Baseline correction
#'
#' Apply a baseline correction to a spectra with the algorithm described in
#' Wang et al. (2013).
#'
#' @param spectra Data frame with spectra in columns and chemical shifts in
#' rows. Colnames of this data frame correspond to pure metabolite names and
#' rownames to chemical shift grid (in ppm).
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#'
#' @return A data frame with baseline corrected spectra in columns and chemical
#' shifts (in ppm) in rows.
#'
#' @references Wang, K.C., Wang, S.Y., Kuo, C.H., Tseng Y.J. (2013).
#' Distribution-based classification method for baseline correction of
#' metabolomic 1D proton nuclear magnetic resonance spectra.
#' \emph{Analytical Chemistry}, \strong{85}(2), 1231-1239.
#'
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SnowParam
#' @importFrom BiocParallel snowWorkers
#' @export
#'
#' @examples
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_base_cor <- baselineCorrection(spectra_data)
baselineCorrection <- function(spectra, ncores = 1){

  if (!is.numeric(as.numeric(rownames(spectra)))) {
    stop("Rownames of spectra data frame do not contain ppm grid.")
  }

  # number of cores
  ncores <- min(ncores, ncol(spectra))
  if (.Platform$OS.type == "windows") {
    para_param <- SnowParam(workers = ncores,
                            progressbar = TRUE,
                            tasks = ncol(spectra))
  } else {
    para_param <- MulticoreParam(workers = ncores,
                                 progressbar = TRUE,
                                 tasks = ncol(spectra))
  }



  spectra_list <- as.list(spectra)

  # baseline correction
  if (ncores > 1) {
    spectra_bc_list <-
      bptry(bplapply(spectra_list, .baselineCorrector,
                     BPPARAM = para_param))
  } else {
    suppressWarnings(
      spectra_bc_list <- llply(spectra_list,
                               function(x)
                                 try(.baselineCorrector(x), silent=TRUE),
                               .progress = "text"))
  }

  # if error in a file
  if (any(!bpok(spectra_bc_list)) |
      any(vapply(spectra_bc_list, is, "try-error", FUN.VALUE = logical(1)))) {
    warning(paste0("The baseline correction algorithm can not be used for spectra: ",
                   paste(colnames(spectra)[!bpok(spectra_bc_list) |
                                             vapply(spectra_bc_list, is,
                                                    "try-error",
                                                    FUN.VALUE = logical(1))],
                         collapse = ", ")), call. = FALSE)
  }
  spectra_list[!(!bpok(spectra_bc_list) | vapply(spectra_bc_list, is,
                                                 "try-error",
                                                 FUN.VALUE = logical(1)))] <-
    spectra_bc_list[!(!bpok(spectra_bc_list) | vapply(spectra_bc_list, is,
                                                      "try-error",
                                                      FUN.VALUE = logical(1)))]

  # convert in a data frame
  spectra_bc <- as.data.frame(do.call(cbind, spectra_list))
  spectra_bc[spectra_bc < 0] <- 0

  spectra_norm <- data.frame(apply(spectra_bc, 2, function(x)
    t(x / .AUC(as.numeric(rownames(spectra)), x))))
  rownames(spectra_norm) <- rownames(spectra)
  colnames(spectra_norm) <- colnames(spectra)

  return(spectra_norm)
}



#' Alignment
#'
#' Align spectra of a data frame with the CluPA algorithm (Vu et al., (2011))
#'
#' @param spectra Data frame with spectra in columns and chemical shift in rows.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in ppm).
#' @param baseline.threshold Value of baseline threshold used to identify peaks.
#' Default to 0.02.
#' @param reference Index of the reference spectrum used for the alignment.
#' Default to \code{NULL}, \emph{i.e.} the reference spectrum is automatically
#' detected with the \code{\link[speaq]{findRef}} function of the \code{speaq}
#' package.
#' @param max.shift Maximum shift allowed for the alignment. Default to 0.002.
#'
#' @return A data frame with aligned spectra in columns and chemical shifts
#' (in ppm) in rows.
#'
#' @references Vu, T. N., Valkenborg, D., Smets, K., Verwaest, K. A., Dommisse,
#' R., Lemiere, F., ... & Laukens, K. (2011). An integrated workflow for robust
#' alignment and simplified quantitative analysis of NMR spectrometry data.
#' \emph{BMC Bioinformatics}, \strong{12}(1), 405.
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
    stop("Rownames of spectra data frame do not contain ppm grid.")
  }

  # maximum shift
  max_shift <- floor(max.shift / (ASICS::pure_library@ppm.grid[2] -
                                    ASICS::pure_library@ppm.grid[1]))

  # detect peaks
  peak_list <- detectSpecPeaks(t(spectra),
                               nDivRange = 64,
                               scales = seq(1, 16, 2),
                               baselineThresh = baseline.threshold,
                               SNR.Th = -1,
                               verbose = FALSE)

  # reference spectrum
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
    t(x / .AUC(as.numeric(rownames(spectra_align)), x))))
  rownames(spectra_norm) <- rownames(spectra)
  colnames(spectra_norm) <- colnames(spectra)

  return(spectra_norm)
}

#' Create a pure library
#'
#' Create a new pure library from a data frame containing different spectra of
#' pure metabolites. The noise is removed by thresholding each spectrum during
#' the creation of a new pure library.
#'
#' @param spectra Data frame with spectra in columns and chemical shifts in
#' rows. Colnames of this data frame correspond to pure metabolite names and
#' rownames to chemical shift grid (in ppm).
#' @param nb.protons Numeric vector of the number of protons for each pure
#' metabolite spectrum contained in \code{spectra} data frame.
#' @param threshold Numeric value below which pure spectrum values are
#' considered to be zero. Default to 1.
#' @return A \linkS4class{PureLibrary} object with the newly created library.
#'
#' @importFrom methods new
#' @export
#'
#' @examples
#' pure_spectra <- importSpectraBruker(system.file("extdata",
#'                                                 "example_library",
#'                                                 package = "ASICS"))
#' new_pure_library <- createPureLibrary(pure_spectra,
#'                                       nb.protons = c(5, 4))
#'
createPureLibrary <- function(spectra, nb.protons, threshold = 1){
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
#' @param spectra Data frame with spectra in columns and chemical shifts in
#' rows. Colnames of this data frame correspond to pure metabolite names and
#' rownames to chemical shift grid (in ppm).
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
#' spectra_data <- importSpectraBruker(current_path)
#' spectra_obj <- createSpectra(spectra_data)
createSpectra <- function(spectra){

  # create a new object Spectra with a data frame of intensities
  new_spectra <- new("Spectra",
                     sample.name = colnames(spectra),
                     ppm.grid = as.numeric(rownames(spectra)),
                     spectra = as.matrix(spectra))

  rownames(new_spectra@spectra) <- NULL
  colnames(new_spectra@spectra) <- NULL

  return(new_spectra)
}




#' Binning/Bucketing of NMR spectra
#'
#' Apply a binning function on a spectrum.
#'
#' @param spectra Data frame with spectra in columns and chemical shifts in
#' rows. Colnames of this data frame correspond to pure metabolite names and
#' rownames to chemical shift grid (in ppm).
#' @param bin Numeric value specifying the bin width.
#' @param exclusion.areas Definition domain of spectra that have to be excluded
#' of the analysis (ppm). By default, the water region is excluded
#' (4.5-5.1 ppm).
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#'
#' @return A data frame with normalised spectra in columns and buckets in rows
#' (bucket names correspond to the center of the bucket).
#'
#' @export
#' @importFrom plyr alply
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SnowParam
#' @importFrom BiocParallel snowWorkers
#'
#' @examples
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_bin <- binning(spectra_data, bin = 0.01)
binning <- function(spectra, bin = 0.01,
                    exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                    ncores = 1){

  if (!is.numeric(as.numeric(rownames(spectra)))) {
    stop("Rownames of spectra data frame don't contain ppm grid.")
  }

  if(!is.null(exclusion.areas) &&
     (!is.matrix(exclusion.areas) | ncol(exclusion.areas) != 2)){
    stop("'exclusion.areas' needs to be a matrix with 2 columns.")
  }

  old_grid <- as.numeric(rownames(spectra))

  # for signal in exclusion.areas, intensity is null (mixture and library)
  if(!is.null(exclusion.areas)){
    idx_to_remove <-
      unlist(alply(exclusion.areas, 1,
                   function(area) which(old_grid >= area[[1]] &
                                          old_grid <= area[[2]])),
             use.names = FALSE)

    spectra[idx_to_remove, ] <- 0
  }

  # number of cores
  ncores <- min(ncores, ncol(spectra))
  if (.Platform$OS.type == "windows") {
    para_param <- SnowParam(workers = ncores,
                            progressbar = TRUE,
                            tasks = ncol(spectra))
  } else {
    para_param <- MulticoreParam(workers = ncores,
                                 progressbar = TRUE,
                                 tasks = ncol(spectra))
  }

  # compute buckets
  buckets <- seq(max(0.5, min(trunc(old_grid * 10, 1) / 10)) + bin / 2,
                 min(10, max(trunc(old_grid * 10, 1) / 10)) + bin / 2,
                 by = bin)


  spectra_list <- as.list(spectra)
  # buckets values for each spectrum
  if (ncores > 1) {
    buckets_values_list <-
      bplapply(spectra_list, .binningSpectrum, old_grid, buckets, bin,
               BPPARAM = para_param)
  } else {
    buckets_values_list <- llply(spectra_list, .binningSpectrum, old_grid,
                                 buckets, bin, .progress = "text")
  }
  # convert in a data frame
  buckets_values <- as.data.frame(do.call(cbind, buckets_values_list))

  # normalisation
  spectra_norm <- data.frame(apply(buckets_values, 2, function(value)
    t(value / .AUC(buckets, value))))
  rownames(spectra_norm) <- buckets
  colnames(spectra_norm) <- colnames(spectra)

  return(spectra_norm)
}


.binningSpectrum <- function(spectrum, grid, buckets, bin){

  buckets_values <-
    vapply(buckets,
           function(x) ifelse(length(grid[grid >= x - bin / 2 &
                                            grid < x + bin / 2]) == 0,
                              0,
                              .AUC(grid[grid >= x - bin / 2 &
                                          grid < x + bin / 2],
                                   spectrum[grid >= x - bin / 2 &
                                              grid < x + bin / 2])),
           FUN.VALUE = numeric(1))

  return(buckets_values)
}

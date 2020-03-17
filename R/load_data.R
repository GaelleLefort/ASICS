#' Import metabolomic spectra
#'
#' Import spectra from text or CSV, fid or 1r (preprocessed spectrum) files.
#' (optional) Spectra are baseline corrected, aligned and normalised during the
#' importation.
#'
#' @param name.dir Path of the folder containing the spectra. Each
#' subfolder contains the fid or the 1r (preprocessed spectrum) files of this
#' sample if \code{type = "fid"} or \code{type = "1r"}.
#' @param name.file Name of the txt or csv file containing the spectra in
#' columns (spectrum names in the first line and ppm grid in the first column).
#' @param type.import Type of import. Either \code{"txt"}, \code{"csv"},
#' \code{"fid"} or \code{"1r"}.
#' @param baseline.correction Logical. If \code{TRUE} a baseline correction is
#' applied for each spectrum. Default to \code{FALSE}.
#' @param alignment Logical. If \code{TRUE} a peak alignment is
#' applied between all spectra. Default to \code{FALSE}.
#' @param normalisation Logical. If \code{TRUE} a normalisation is applied for
#' each spectrum (see \code{\link{normaliseSpectra}} for details). Default to
#' \code{TRUE}.
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#' @param verbose A boolean value to allow print out process information.
#' @param ... Further arguments to be passed to the functions
#' \code{\link{read.table}}, \code{\link{importSpectraBruker}},
#' \code{\link[PepsNMR]{Normalization}}
#' (\code{\link[PepsNMR]{PepsNMR-package}}), \code{\link{alignSpectra}} or
#' \code{\link{normaliseSpectra}} for
#' specifying the parameters of the algorithm, if necessary.
#'
#' @details
#' Some preprocessing steps are included during the importation. First, spectra
#' are baseline corrected if \code{baseline.correction = TRUE}. Then, all
#' spectrum definition domains are aligned to a unique one (either the one
#' specified in \code{ppm.grid} or the grid of the default library). Finally,
#' all spectra are normalised if \code{normalisation = TRUE} and aligned if
#' \code{alignment = TRUE}.
#'
#' @references Wang, K.C., Wang, S.Y., Kuo, C.H., Tseng Y.J. (2013).
#' Distribution-based classification method for baseline correction of
#' metabolomic 1D proton nuclear magnetic resonance spectra.
#' \emph{Analytical Chemistry}, \strong{85}(2), 1231-1239.
#'
#' @return A data frame with spectra in columns and chemical shifts (in ppm)
#' in rows.
#'
#' @import PepsNMR
#' @importFrom utils read.table capture.output
#'
#' @export
#'
#' @seealso \code{\link{importSpectraBruker}}
#' \code{\link{normaliseSpectra}} \code{\link{alignSpectra}}
#'
#' @examples
#' ## Import from txt file
#' current_path <- system.file("extdata", package = "ASICS")
#' spectra_data <- importSpectra(name.dir = current_path,
#'                      name.file = "spectra_example.txt", type.import = "txt")
#'
#' ## Import from fid file
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- importSpectra(name.dir = current_path, type.import = "fid",
#'                              subdirs = TRUE, dirs.names = TRUE)
#'
#' ## Import from txt file
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- importSpectra(name.dir = current_path, type.import = "1r")
importSpectra <- function(name.dir = NULL, name.file = NULL, type.import,
                          baseline.correction = FALSE, alignment = FALSE,
                          normalisation = TRUE, ncores = 1,
                          verbose = TRUE, ...){

  import.args <- list(...)

  if (!type.import %in% c("txt", "csv", "fid", "1r")) {
    stop("This import method is not available.")
  }

  if (type.import == "txt" | type.import == "csv") {
    if (verbose) cat("Import spectra from txt or csv files... \n")
    import.param <- c(list(file = file.path(name.dir, name.file)),
                      import.args[names(import.args) %in%
                                    names(formals(args(read.table)))])

    imported_spectra <- do.call("read.table", import.param)
  } else if (type.import == "fid") {
    if (verbose) cat("Import spectra from fid files...  \n")

    #zero-filling
    zf <- TRUE
    if(!is.null(import.args[["fn"]]) && import.args["fn"] < 1) zf <- FALSE

    import.param <- c(list(Fid_data = NULL, Fid_info = NULL,
                           data.path = name.dir, readFids = TRUE,
                           groupDelayCorr = TRUE, solventSuppression = TRUE,
                           apodization = TRUE, zerofilling = zf,
                           fourierTransform = TRUE, zeroOrderPhaseCorr = TRUE,
                           internalReferencing = TRUE,
                           baselineCorrection = FALSE, negativeValues0 = FALSE,
                           warping = FALSE, windowSelection = TRUE,
                           bucketing = FALSE, regionRemoval = FALSE,
                           zoneAggregation = FALSE, normalization = FALSE,
                           export = FALSE, writeArg = "none", verbose = verbose),
                      import.args[names(import.args) %in%
                                    names(formals(args(ReadFids)))],
                      import.args[names(import.args) %in%
                                    names(formals(args(GroupDelayCorrection)))],
                      import.args[names(import.args) %in%
                                    names(formals(args(SolventSuppression)))],
                      import.args[names(import.args) %in%
                                    names(formals(args(Apodization)))],
                      import.args[names(import.args) %in%
                                    names(formals(args(ZeroFilling)))],
                      import.args[names(import.args) %in%
                                    names(formals(args(FourierTransform)))],
                      import.args[names(import.args) %in%
                                names(formals(args(ZeroOrderPhaseCorrection)))],
                      import.args[names(import.args) %in%
                                    names(formals(args(InternalReferencing)))],
                      import.args[names(import.args) %in%
                                  names(formals(args(NegativeValuesZeroing)))],
                      import.args[names(import.args) %in%
                                    names(formals(args(WindowSelection)))])

    if (verbose) {
      imported_spectra <-
        do.call("PreprocessingChain", import.param)$Spectrum_data
    } else {
      imported_spectra <- do.call("PreprocessingChain",
                                  import.param)$Spectrum_data
    }
    imported_spectra <- as.data.frame(t(Re(imported_spectra)))
    imported_spectra <- imported_spectra[nrow(imported_spectra):1, ]
  } else if (type.import == "1r") {
    if (verbose) cat("Import spectra from 1r files... \n")
    import.param <- c(list(name.dir = name.dir,
                           ncores = ncores,
                           verbose = verbose),
                      import.args[names(import.args) %in%
                                    names(formals(args(importSpectraBruker)))])

    imported_spectra <- do.call("importSpectraBruker", import.param)
  }

  # baseline correction
  if (baseline.correction) {
    baseline.params <- c(list(Spectrum_data = t(imported_spectra),
                              verbose = verbose),
                         import.args[names(import.args) %in%
                                       c(names(formals(args(BaselineCorrection))))])

    imported_spectra <- t(do.call("BaselineCorrection", baseline.params))

    negvalue.params <- c(list(Spectrum_data = t(imported_spectra),
                              verbose = verbose))

    imported_spectra <- t(do.call("NegativeValuesZeroing", negvalue.params))
  }

  # peak alignment
  if (alignment) {
    align.param <- c(list(spectra = imported_spectra,
                          ncores = ncores,
                          verbose = verbose),
                     import.args[names(import.args) %in%
                                   c(names(formals(args(alignSpectra))))])

    imported_spectra <- do.call("alignSpectra", align.param)
  }

  # normalisation
  if (normalisation) {
    norm.param <- c(list(spectra = imported_spectra,
                         verbose = verbose),
                    import.args[names(import.args) %in%
                                  c(names(formals(args(Normalization))))])

    imported_spectra <- do.call("normaliseSpectra", norm.param)
    attr(imported_spectra, "type.norm") <-
      ifelse(!is.null(norm.param$type.norm), norm.param$type.norm, "CS")
    params <- norm.param[!(names(norm.param) %in% c("spectra", "verbose",
                                                    "type.norm"))]
    if (!is.null(norm.param$type.norm) & length(params) != 0) {
      attr(imported_spectra, "norm.param") <- params
    }
  }

  return(imported_spectra)
}




#' Import preprocessed metabolomic spectra from Bruker files
#'
#' Import preprocessed spectra from Bruker files contained in a single folder.
#' This folder contains one subfolder for each sample.
#' (optional) Spectra are baseline corrected, aligned and normalised by the area
#' under the curve during the importation.
#'
#' @param name.dir Path of the folder containing one subfolder by sample. Each
#' subfolder contains the Bruker files of this sample.
#' @param which.spectra If there is no folder with experiment number
#' (\code{all_spectra/<spectrum_name>/pdata/...}) set \code{which.spectra} to
#' \code{NULL}. Else if there is more than one spectrum by sample
#' (\code{all_spectra/<spectrum_name>/<experiment_number>/pdata/...}), which is
#' the spectrum to import (either always the first one with \code{which.spectra
#' = "first"}, always the last one with \code{which.spectra = "last"} or a
#' vector of length the number of spectra that specifies the number of each
#' spectrum to import).
#' Default to \code{"first"}.
#' @param ppm.grid Numeric vector of a unique grid (definition domain) for all
#' spectra (in ppm). Default to \code{NULL} (in which case, the default grid of
#' the pure library is used).
#' @param sample.names Character vector of sample names. Default to \code{NULL}
#' (in which case, folder names are used).
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#' @param verbose A boolean value to allow print out process information.
#'
#' @return A data frame with spectra in columns and chemical shifts (in ppm)
#' in rows.
#'
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SerialParam
#' @importFrom BiocParallel bptry bpok
#' @export
#'
#' @seealso \code{\link{normaliseSpectra}}
#' \code{\link{alignSpectra}}
#'
#' @examples
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- importSpectraBruker(current_path)
importSpectraBruker <- function(name.dir, which.spectra = "first",
                                ppm.grid = NULL, sample.names = NULL,
                                ncores = 1, verbose = TRUE){

  if(!dir.exists(name.dir)){
    stop("Path of the Bruker files doesn't exist!")
  }

  if (!is.null(which.spectra) && (length(which.spectra) != 1 |
                                  (!(which.spectra %in% c("first", "last", NULL)) &
                                    !is.numeric(which.spectra)))) {
    stop(paste("'which.spectra' must be of length 1 and either 'first', 'last'",
               "'NULL' or a number."))
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
                                        utils::head, 1,
                                        FUN.VALUE = numeric(1))))
  } else if (length(which.spectra) == 1 && which.spectra == "last") {
    cur_dir_spec <- as.list(file.path(cur_dir,
                                      vapply(as.list(lapply(cur_dir, function(x)
                                        sort(as.numeric(dir(x))))),
                                        utils::tail, 1,
                                        FUN.VALUE = numeric(1))))
  } else if (is.null(which.spectra)) {
    cur_dir_spec <- cur_dir
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
  if (.Platform$OS.type == "windows" | ncores == 1) {
    para_param <- SerialParam(progressbar = verbose,
                              stop.on.error = FALSE)
  } else {
    para_param <- MulticoreParam(workers = ncores,
                                 progressbar = verbose,
                                 tasks = length(cur_dir_spec),
                                 manager.hostname = "localhost",
                                 stop.on.error = FALSE)
  }

  # import spectra
  imported_spectra <-
    bptry(bplapply(cur_dir_spec, .readNMRBruker, ppm.grid,
                   BPPARAM = para_param))

  # if error in a file
  if (any(!bpok(imported_spectra))) {
    warning(paste0("There is a problem in files for spectra: ",
                   paste(sample.names[!bpok(imported_spectra)],
                         collapse = ", "),
                   ".\nFix the problem manually and re-execute the function ",
                   "otherwise these spectra are ignored."), call. = FALSE)
    sample.names <- sample.names[bpok(imported_spectra)]
    imported_spectra <- imported_spectra[bpok(imported_spectra)]
  }

  # convert in a data frame
  imported_spectra <- as.data.frame(do.call(cbind, imported_spectra),
                                    row.names = as.character(ppm.grid))
  colnames(imported_spectra) <- sample.names

  return(imported_spectra)
}


## Extract a spectrum from Bruker files
.readNMRBruker <- function(path, ppm.grid){

  # file paths
  acq_file <- file.path(path, "acqus")
  spec_file <- file.path(path, "pdata/1/1r")
  proc_file <- file.path(path, "pdata/1/procs")

  # meta-data in acq_file
  ACQ <- readLines(acq_file)
  TD <- .getBrukerParam(ACQ, "TD")
  SW <- .getBrukerParam(ACQ, "SW")
  DTYPA <- .getBrukerParam(ACQ, "DTYPA")
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

  # adapt the grid
  new_signal <- .changeGrid(signal, ppmseq, ppm.grid)

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
#' Normalise a data frame of spectra to a constant sum (CS) or with a
#' method of \code{PepsNMR} package (see \code{\link[PepsNMR]{Normalization}}).
#'
#' @param spectra Data frame with spectra in columns and chemical shifts in
#' rows. Colnames of this data frame correspond to pure metabolite names and
#' rownames to chemical shift grid (in ppm).
#' @param type.norm Type of normalisation : "CS", "mean", "pqn", "median",
#' "firstquartile" or "peak". Default to "CS".
#' @param verbose A boolean value to allow print out process information.
#' @param ... other arguments to be passed to
#' \code{\link[PepsNMR]{Normalization}}
#'
#' @return A data frame with normalised spectra in columns and chemical shifts
#' (in ppm) in rows.
#'
#' @importFrom utils capture.output
#' @export
#'
#' @examples
#' current_path <- system.file("extdata", package = "ASICS")
#' spectra_data <- importSpectra(name.dir = current_path,
#'                      name.file = "spectra_example.txt", type.import = "txt")
#' spectra_norm <- normaliseSpectra(spectra_data, type.norm = "pqn")
normaliseSpectra <- function(spectra, type.norm = "CS", verbose = TRUE, ...){

  if (!type.norm %in% c("CS", "mean", "pqn", "median", "firstquartile",
                        "peak")) {
    stop("This normalisation method is not available.")
  }

  if (!is.numeric(as.numeric(rownames(spectra)))) {
    stop("Rownames of spectra data frame don't contain ppm grid.")
  }

  if (verbose) cat("Normalisation... \n")

  if (type.norm == "CS") {
    spectra_norm <- data.frame(apply(spectra, 2, function(x)
      t(x / .AUC(as.numeric(rownames(spectra)), x))))
  } else {
    capture.output(spectra_norm  <-
                     data.frame(t(Normalization(t(spectra),
                                                type.norm = type.norm, ...))))
  }

  rownames(spectra_norm) <- rownames(spectra)
  colnames(spectra_norm) <- colnames(spectra)

  return(spectra_norm)
}


#' Alignment
#'
#' Align spectra of a data frame by a method based on the CluPA algorithm
#' (Vu et al., (2011))
#'
#' @param spectra Data frame with spectra in columns and chemical shift in rows.
#' Colnames of this data frame correspond to pure metabolite names and rownames
#' to chemical shift grid (in ppm).
#' @param reference Index of the reference spectrum used for the alignment.
#' Default to \code{NULL}, \emph{i.e.} the reference spectrum is automatically
#' detected.
#' @param max.shift Maximum shift allowed for the alignment. Default to 0.002.
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#' @param verbose A boolean value to allow print out process information.
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
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SerialParam
#'
#' @examples
#' current_path <- system.file("extdata", package = "ASICS")
#' spectra_data <- importSpectra(name.dir = current_path,
#'                      name.file = "spectra_example.txt", type.import = "txt")
#' spectra_align <- alignSpectra(spectra_data)
alignSpectra <- function(spectra, reference = NULL, max.shift = 0.02,
                         ncores = 1, verbose = TRUE){

  if (!is.numeric(as.numeric(rownames(spectra)))) {
    stop("Rownames of spectra data frame do not contain ppm grid.")
  }

  if (verbose) cat("Alignment... \n")

  # maximum shift
  max_shift <- floor(max.shift / (ASICS::pure_library@ppm.grid[2] -
                                    ASICS::pure_library@ppm.grid[1]))

  # detect peaks
  if (verbose) cat("Peak detection \n")
  peak_list <- bplapply(as.list(seq_len(ncol(spectra))),
                        .findPeaks, spectra,
                        BPPARAM = .createEnv(ncores, ncol(spectra), verbose))


  # reference spectrum
  if (is.null(reference)) {
    if (verbose) cat("Finding reference spectrum \n")
    reference <- .findReference(spectra, ncores, verbose)

    if (verbose) cat(paste("The reference spectrum is the number", reference,
                           ":", colnames(spectra)[reference],"\n"))
  }

  # import spectra
  if (verbose) cat("Align spectra \n")
  align_spectra <-
    bplapply(as.list(seq_len(ncol(spectra))),
             .alignmentIdx, spectra, peak_list, reference, max_shift,
             BPPARAM = .createEnv(ncores, ncol(spectra), verbose))

  align_spectra_df <- as.data.frame(do.call(cbind, align_spectra),
                                    row.names = rownames(spectra))
  colnames(align_spectra_df) <- colnames(spectra)

  return(align_spectra_df)
}


.alignmentIdx <- function(idx, spectra, peak_list, reference, max_shift) {
  spectrum_align <- as.numeric(spectra[, idx])
  if (idx != reference) {
    spectrum_align <- .doAlignment(spectra[, idx], peak_list[[idx]],
                                   spectra[, reference],
                                   peak_list[[reference]],
                                   max_shift,
                                   acceptLostPeak = TRUE)
  }

  names(spectrum_align) <- NULL
  return(spectrum_align)
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
#' @param threshold Numeric value or numeric vector of length
#' \code{ncol(spectra)} below which pure spectrum values are considered to be
#' zero. Default to 1.
#' @return A \linkS4class{PureLibrary} object with the newly created library.
#'
#' @importFrom methods new
#' @importFrom Matrix Matrix
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
                     spectra = Matrix(as.matrix(spectra)),
                     nb.protons = nb.protons)
  colnames(new_library@spectra) <- NULL
  rownames(new_library@spectra) <- NULL

  # spectra thresholding
  if(length(threshold) == 1) {
    new_library@spectra[new_library@spectra < threshold] <- 0
  } else {
    for(i in seq_along(threshold)) {
      new_library@spectra[new_library@spectra[, i] < threshold[i], i] <- 0
    }
  }

  return(new_library)
}



#' Create a \linkS4class{Spectra} object
#'
#' Create a new spectra object used for quantification.
#'
#' @param spectra Data frame with spectra in columns and chemical shifts in
#' rows. Colnames of this data frame correspond to pure metabolite names and
#' rownames to chemical shift grid (in ppm).
#' @param norm.method Character specifying the normalisation method to use on
#' spectra ONLY if the \code{\link{importSpectra}} function was not used.
#' @param norm.params List containing normalisation parameteres (see
#' \code{\link{normaliseSpectra}} for details) ONLY if the
#' \code{\link{importSpectra}} function was not used.
#'
#' @return A \linkS4class{Spectra} object with spectra to quantify.
#'
#' @export
#' @importFrom BiocParallel bplapply multicoreWorkers
#' @importFrom Matrix Matrix
#'
#' @seealso \linkS4class{Spectra}
#'
#' @examples
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- importSpectraBruker(current_path)
#' spectra_obj <- createSpectra(spectra_data)
createSpectra <- function(spectra, norm.method = NULL, norm.params = NULL){

  if (!is.null(attr(spectra, "type.norm")) && !is.null(norm.method) &&
               attr(spectra, "type.norm") != norm.method) {
    warning(paste("Normalisation method is not the same than the one use",
                  "during data import."))
  }

  if (is.null(norm.method) & !is.null(attr(spectra, "type.norm"))) {
    norm.method <- attr(spectra, "type.norm")
  } else if (is.null(norm.method) & is.null(attr(spectra, "type.norm"))) {
    norm.method <- "CS"
  }

  if (is.null(norm.params) & !is.null(attr(spectra, "norm.param"))) {
    norm.params <- attr(spectra, "norm.param")
  } else if (is.null(norm.params) & is.null(attr(spectra, "norm.param"))) {
    norm.params <- list()
  }

  # create a new object Spectra with a data frame of intensities
  new_spectra <- new("Spectra",
                     sample.name = colnames(spectra),
                     ppm.grid = as.numeric(rownames(spectra)),
                     spectra = Matrix(as.matrix(spectra)),
                     norm.method = norm.method,
                     norm.params = norm.params)

  rownames(new_spectra@spectra) <- NULL
  colnames(new_spectra@spectra) <- NULL

  return(new_spectra)
}




#' Binning/Bucketing of NMR spectra
#'
#' Apply a binning function on a spectrum.
#'
#' @param spectra Data frame with spectra in columns and chemical shifts in
#' rows. Colnames of this data frame correspond to sample names and
#' rownames to chemical shift grid (in ppm).
#' @param bin Numeric value specifying the bin width.
#' @param exclusion.areas Definition domain of spectra that have to be excluded
#' of the analysis (ppm). By default, the water region is excluded
#' (4.5-5.1 ppm).
#' @param normalisation Logical. If \code{TRUE} a normalisation is applied for
#' each spectrum (see \code{\link{normaliseSpectra}} for details). Default to
#' \code{TRUE}.
#' @param ncores Number of cores used in parallel evaluation. Default to
#' \code{1}.
#' @param verbose A boolean value to allow print out process information.
#' @param ... Further arguments to be passed to the function
#' \code{\link{normaliseSpectra}}
#'
#' @return A data frame with normalised spectra in columns and buckets in rows
#' (bucket names correspond to the center of the bucket).
#'
#' @export
#' @importFrom plyr alply
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers SerialParam
#'
#' @examples
#' current_path <- system.file("extdata", package = "ASICS")
#' spectra_data <- importSpectra(name.dir = current_path,
#'                      name.file = "spectra_example.txt", type.import = "txt")
#' spectra_bin <- binning(spectra_data, bin = 0.01, type.norm = "pqn")
binning <- function(spectra, bin = 0.01,
                    exclusion.areas = matrix(c(4.5, 5.1), ncol = 2),
                    normalisation= TRUE, ncores = 1, verbose = TRUE, ...){

  extra.args <- list(...)

  if (!is.numeric(as.numeric(rownames(spectra)))) {
    stop("Rownames of spectra data frame don't contain ppm grid.")
  }

  if(!is.null(exclusion.areas) &&
     (!is.matrix(exclusion.areas) | ncol(exclusion.areas) != 2)){
    stop("'exclusion.areas' needs to be a matrix with 2 columns.")
  }

  if (verbose) cat("Binning \n")

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

  # compute buckets
  buckets <- seq(max(0.5, min(trunc(old_grid * 10, 1) / 10)) + bin / 2,
                 min(10, max(trunc(old_grid * 10, 1) / 10)) + bin / 2,
                 by = bin)


  spectra_list <- as.list(spectra)
  # buckets values for each spectrum
  buckets_values_list <-
    bplapply(spectra_list, .binningSpectrum, old_grid, buckets, bin,
             BPPARAM = .createEnv(ncores, ncol(spectra), verbose))

  # convert in a data frame
  buckets_values <- as.data.frame(do.call(cbind, buckets_values_list))

  # normalisation
  # normalisation
  spectra_norm <- buckets_values
  if (normalisation) {
    norm.param <- c(list(spectra = buckets_values,
                         verbose = verbose),
                    extra.args[names(extra.args) %in%
                                 c(names(formals(args(Normalization))))])

    spectra_norm <- do.call("normaliseSpectra", norm.param)
  }

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

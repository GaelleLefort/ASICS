## Main function to load library, import a spectrum from Bruker files and
#perform the first preprocessing
load_data <- function(path, library.metabolites = NULL,
                      exclusion.areas = matrix(c(4.5,5.1), ncol = 2, nrow = 1),
                      max.shift = 0.02, which.spectra = "last"){

  #load default or user metabolites library
  if(!is.null(library.metabolites)) {
    if(!file.exists(library.metabolites)){
      stop("Pure library file doesn't exist !")
    }

    load(library.metabolites)
  }

  ## checking the validity of the library
  if(sum(names(pure_library) %in% c("grid", "name", "spectra", "nb_protons"))
     != 4){
    stop("grid, name, spectra and nb_protons must be elements of pure library
         list.")
  }

  if(length(pure_library$grid) != nrow(pure_library$spectra) |
     ncol(pure_library$spectra) != length(pure_library$name) |
     ncol(pure_library$spectra) != length(pure_library$nb_protons)){
    stop("Pure library is not well formatted.")
  }

  #get the mixture
  mixture <- get_mixture(path, pure_library$grid, which.spectra)

  #for signal in exclusion.areas, intensity is null (mixture and library)
  idx_to_remove <-
    unlist(plyr::alply(exclusion.areas, 1, function(x)
      which(pure_library$grid >= x[[1]] & pure_library$grid <= x[[2]])),
      use.names = FALSE)

  mixture[idx_to_remove] <- 0
  pure_library$spectra[idx_to_remove, ] <- 0

  #remove metabolite without any signal
  with_signal <- as.numeric(which(colSums(pure_library$spectra) > 0))
  cleaned_library <- subset_library(pure_library, with_signal)

  #re-normalisation of mixture and pure spectra library
  cleaned_library$spectra <- apply(cleaned_library$spectra, 2,
                                 function(x)
                                   t(x / AUC(cleaned_library$grid, x)))
  mixture <- mixture / AUC(cleaned_library$grid, mixture)

  return(list("mixture" = mixture, "pure_library" = cleaned_library))

}


## Get the mixture from files, correct baseline and normalize it
get_mixture <- function(path, library.grid, which.spectra = "last") {

  #path of the chosen spectrum
  if(which.spectra == "first"){
    path_spectrum <- file.path(path, "1")
  } else if(which.spectra == "last"){
    path_spectrum <- file.path(path, dir(path)[length(dir(path))])
  } else {
    path_spectrum <- file.path(path, which.spectra)
  }

  #import mixture from bruker files
  imported_mixture <- read_NMR_bruker(path_spectrum)

  #baseline correction
  imported_mixture[, 2] <- baseline_corrector(imported_mixture[, 2])

  #baseline improvement of complex mixture : remove negative values
  imported_mixture[, 2][imported_mixture[, 2] < 0] <- 0

  #normalisation by area under the curve
  imported_mixture[, 2] <- imported_mixture[, 2] /
    AUC(imported_mixture[, 1], imported_mixture[, 2])

  #obtain the same grid for mixture than for metabolites library
  mixture <- f_o_phi(imported_mixture[, 1], imported_mixture[, 2], library.grid)

  return(mixture)
}


## Extract a spectrum from Bruker files
read_NMR_bruker <- function(path){

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

  sample_spectrum <- cbind(ppmseq, signal)
  return(sample_spectrum)
}

# find parameter in a bruker file
get_bruker_param <- function(file, paramStr){
  regexpStr <- paste0("^...", paramStr, "=")
  as.numeric(gsub("^[^=]+= ", "" ,
                  file[which(simplify2array(regexpr(regexpStr, file)) > 0)]))
}


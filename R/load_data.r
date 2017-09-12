## Main function to load library, import a spectrum from Bruker files and
#perform the first preprocessing
load_data <- function(path, library.metabolites = NULL,
                      exclusion.areas = matrix(c(4.5,5.1), ncol = 2, nrow = 1),
                      max.shift = 0.02, which.spectra = "last"){

  #load default or user metabolites library
  if(is.null(library.metabolites)) {
    data(Library)
  } else {
    load(library.metabolites)
  }

  #get the mixture
  mixture <- get_mixture(path, pure_library$grid, which.spectra)

  #for signal in exclusion.areas, intesity is null (mixture and library)
  idx_to_remove <- unlist(alply(exclusion.areas, 1,
                                function(x) which(pure_library$grid >= x[[1]] &
                                                    pure_library$grid <= x[[2]])),
                          use.names = FALSE)

  mixture[idx_to_remove] <- 0
  pure_library$spectra[idx_to_remove, ] <- 0

  #remove metabolite without any signal
  with_signal <- as.numeric(which(colSums(pure_library$spectra) > 0))
  cleaned_library <- subset_library(pure_library, with_signal)

  #re-normalisation of mixture and pure spectra library
  cleaned_library$spectra <- apply(cleaned_library$spectra, 2,
                                 function(x) t(x / AUC(cleaned_library$grid, x)))
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

  #normalisation by area under the curve
  imported_mixture[, 2] <- imported_mixture[, 2] /
    AUC(imported_mixture[, 1], imported_mixture[, 2])

  #obtain the same grid for mixture than for metabolites library
  mixture <- f_o_phi(imported_mixture[, 1], imported_mixture[, 2], library.grid)

  return(mixture)
}


## Extract a spectrum from Bruker files
###### !! A revoir sauf si ça vient d'un script déjà écrit
read_NMR_bruker <- function(path){

  ACQFILE <- file.path(path, "acqus")
  SPECFILE <- file.path(path, "pdata/1/1r")
  PROCFILE <- file.path(path, "pdata/1/procs")

  ACQ <- readLines(ACQFILE)
  TD      <- get_bruker_param(ACQ, "TD")
  SW      <- get_bruker_param(ACQ, "SW")
  SWH     <- get_bruker_param(ACQ, "SW_h")
  DTYPA   <- get_bruker_param(ACQ, "DTYPA")
  BYTORDA <- get_bruker_param(ACQ, "BYTORDA")
  ENDIAN <- "little"
  SIZE <- ifelse(DTYPA == 0, 4, 8)

  PROC <- readLines(PROCFILE)
  OFFSET <- get_bruker_param(PROC, "OFFSET")
  SI <- get_bruker_param(PROC, "SI")

  to.read <- file(SPECFILE, "rb")
  maxTDSI <- max(TD, SI)
  signal <- rev(readBin(to.read, what = "int", size = SIZE, n = maxTDSI,
                        signed = TRUE, endian = ENDIAN))
  close(to.read)

  td <- length(signal)

  dppm <- SW / (td - 1)
  pmax <- OFFSET
  pmin <- OFFSET - SW
  ppmseq <- seq(from = pmin, to = pmax, by = dppm)
  signal <- 100 * signal / max(signal)

  SampleSpectrum <- cbind(ppmseq, signal)
  return(SampleSpectrum)
}


get_bruker_param <- function(ACQ, paramStr){
  regexpStr <- paste0("^...", paramStr, "=")
  as.numeric(gsub("^[^=]+= ", "" ,
                  ACQ[which(simplify2array(regexpr(regexpStr, ACQ)) > 0)]))
}


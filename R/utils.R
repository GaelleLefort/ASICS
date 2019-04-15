## Compute the area under the curve (y = f(x)) using the trapezoidal rule
.AUC <- function(x, y) {
  auc <- sum(diff(x) * (y[seq_len(length(y) - 1)] + y[2:length(y)]) / 2)
  return(auc)
}

## Linear interpolation to adapt old_spectrum on the old grid to the new grid
.changeGrid <- function(old_spectrum, old_grid, new_grid) {

  # new_grid must be included in old_grid
  if (new_grid[1] < old_grid[1]) new_grid[1] <- old_grid[1]
  n <- length(old_grid)
  nphi <- length(new_grid)
  if (new_grid[nphi] > old_grid[n]) new_grid[nphi] <- old_grid[n]

  # combine both grid (1 = old and 0 = new)
  grid_indicator <- rep(c(1, 0), c(n, nphi))
  old_grid[1] <- old_grid[1] - 1e-06
  old_grid[n] <- old_grid[n] + 1e-06
  combined_grid <- c(old_grid, new_grid)

  # sort combined grid and get points on old grid just before points on new grid
  idx_sorted_grid <- order(combined_grid)
  u <- cumsum(grid_indicator[idx_sorted_grid])
  lower_bound <- u[grid_indicator[idx_sorted_grid] == 0]

  # spectrum on new grid
  new_spectrum <- old_spectrum[lower_bound] +
    (new_grid - old_grid[lower_bound]) *
    (old_spectrum[(lower_bound + 1)] - old_spectrum[lower_bound]) /
    (old_grid[(lower_bound + 1)] - old_grid[lower_bound])

  if (is.na(new_spectrum[nphi])) new_spectrum[nphi] <- old_spectrum[n]
  if (is.na(new_spectrum[1])) new_spectrum[1] <- old_spectrum[1]

  return(new_spectrum)
}



## Main function to load library, import a spectrum from Bruker files and
#perform the first preprocessing
#' @importFrom plyr alply
#' @importFrom Matrix colSums
.removeAreas <- function(spectrum_obj, exclusion.areas, cleaned_library){

  # remove extremities of library grid if the spectrum is shorter
  cleaned_library@spectra <-
    cleaned_library@spectra[cleaned_library@ppm.grid >=
                              min(spectrum_obj@ppm.grid) &
                              cleaned_library@ppm.grid <=
                              max(spectrum_obj@ppm.grid), ]
  cleaned_library@ppm.grid <-
    cleaned_library@ppm.grid[cleaned_library@ppm.grid >=
                               min(spectrum_obj@ppm.grid) &
                               cleaned_library@ppm.grid <=
                               max(spectrum_obj@ppm.grid)]


  # adapt spectrum grid to have the same than library
  spectrum_obj@spectra <- Matrix(.changeGrid(spectrum_obj@spectra,
                                             spectrum_obj@ppm.grid,
                                             cleaned_library@ppm.grid))
  spectrum_obj@ppm.grid <- cleaned_library@ppm.grid


  # for signal in exclusion.areas, intensity is null (mixture and library)
  if(!is.null(exclusion.areas)){
    idx_to_remove <-
      unlist(alply(exclusion.areas, 1,
                   function(x) which(cleaned_library@ppm.grid >= x[[1]] &
                                       cleaned_library@ppm.grid <= x[[2]])),
             use.names = FALSE)

    if (length(idx_to_remove) != 0) {
      spectrum_obj@spectra[idx_to_remove, ] <- 0
      cleaned_library@spectra[idx_to_remove, ] <- 0
    }
  }

  # remove metabolite without any signal
  with_signal <- as.numeric(which(colSums(cleaned_library@spectra) > 0))
  cleaned_library <- cleaned_library[with_signal]

  # re-normalisation of mixture and pure spectra library
  cleaned_library@spectra <- Matrix(apply(cleaned_library@spectra, 2,
                                   function(x)
                                     t(x / .AUC(cleaned_library@ppm.grid, x))))

  spectra_to_norm <- as.data.frame(as.matrix(spectrum_obj@spectra))
  rownames(spectra_to_norm) <- spectrum_obj@ppm.grid
  norm.param <- c(list(spectra = spectra_to_norm,
                       verbose = FALSE,
                       type.norm = spectrum_obj@norm.method),
                  spectrum_obj@norm.params)
  spectrum_obj@spectra <-
    Matrix(as.matrix(do.call("normaliseSpectra", norm.param)))

  return(list("cleaned_spectrum" = spectrum_obj,
              "cleaned_library" = cleaned_library))

}


## Remove metabolites that cannot belong to the mixture
#' @importFrom plyr aaply
.cleanLibrary <- function(obj, threshold.noise, nb_points_shift){

  #support of mixture superior to threshold
  signal_mixture <- 1 * (obj[["cleaned_spectrum"]]@spectra > threshold.noise)
  signal_mixture_shift <- .signalWithShift(signal_mixture, nb_points_shift)

  #normalize each spectra with the maximum of mixture
  norm_library <- t(aaply(obj[["cleaned_library"]]@spectra, 2, function(x) x *
                            max(obj[["cleaned_spectrum"]]@spectra) / max(x)))
  signal_library <- 1 * (norm_library > threshold.noise)

  #keep metabolites for which signal is included in mixture signal
  metab_to_keep <- which(apply(signal_library, 2,
                               function(x)
                                 sum(x - signal_mixture_shift == 1) == 0))

  obj[["cleaned_library"]] <- obj[["cleaned_library"]][metab_to_keep]
  return(obj)
}


## Extend each peak in signal of nb_points_shift
.signalWithShift <- function(signal, nb_points_shift) {

  # get indexes with a signal
  idx_signal <- which(signal == 1)

  # get extremities of peaks
  min_extremities <- idx_signal[!((idx_signal - 1) %in% idx_signal)]
  max_extremities <- idx_signal[!((idx_signal + 1) %in% idx_signal)]

  # new signal indexes with shift
  shift_before <- unlist(lapply(min_extremities,
                                function(x) max(x - nb_points_shift, 1):x))
  shift_after <- unlist(lapply(max_extremities,
                               function(x) x:min(x + nb_points_shift,
                                                 length(signal))))

  # new signal
  signal[c(shift_before, shift_after)] <- 1

  return(signal)
}


## Non-negative least square (formula : y~x with weights w)
#' @importFrom quadprog solve.QP
.lmConstrained <- function(y, x, w = 1, precision = 1){
  # constraints
  Amat <- diag(rep(1, ncol(x)))
  b0 <- rep(0, ncol(x))

  # quadratic function to be minimized
  D <- t(x) %*% (w * x) * precision
  dlittle <- t(x) %*% (w * y) * precision

  # optimisation
  res_lm <- solve.QP(D, dlittle, Amat, b0)

  # coefficients
  coefficients <- res_lm$solution
  coefficients[coefficients < 0] <- 0

  # residuals
  residuals <- y - x %*% coefficients

  return(list(coefficients = coefficients, residuals = residuals))
}


## To create a parallel environment
#' @importFrom BiocParallel SerialParam MulticoreParam register
.createEnv <- function(ncores, ntasks, verbose){
  ncores <- min(ncores, ntasks)

  para_param <- NULL
  if (.Platform$OS.type == "windows" | ncores == 1) {
    para_param <- SerialParam(progressbar = verbose)
  } else {
    para_param <- MulticoreParam(workers = ncores,
                                 progressbar = verbose,
                                 tasks = ntasks)
  }
  return(para_param)
}


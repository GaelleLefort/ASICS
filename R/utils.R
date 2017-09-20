## Compute the area under the curve (y = f(x)) using the trapezoidal rule
AUC <- function(x, y) {
  auc <- sum(diff(x) * zoo::rollmean(y, 2))
  return(auc)
}


## Create a subset of metabolites library from metabolites indexes
subset_library <- function(pure_library, idx){
  sub_library <- list()
  sub_library$name <- pure_library$name[idx]
  sub_library$grid <- pure_library$grid
  sub_library$spectra <- pure_library$spectra[, idx]
  sub_library$nb_protons <- pure_library$nb_protons[idx]
  return(sub_library)
}


## remove metabolites that cannot belong to the mixture
clean_library <- function(mixture, pure_library, threshold.noise,
                          nb_points_shift){
  #support of mixture superior to threshold
  signal_mixture <- 1 * (mixture > threshold.noise)
  signal_mixture_shift <- signal_with_shift(signal_mixture, nb_points_shift)

  #normalize each spectra with the maximum of mixture
  norm_library <- t(plyr::aaply(pure_library$spectra, 2, function(x) x *
                                  max(mixture) / max(x)))
  signal_library <- 1 * (norm_library > threshold.noise)

  #keep metabolites for which signal is included in mixture signal
  metab_to_keep <- which(apply(signal_library, 2,
                               function(x)
                                 sum(x - signal_mixture_shift == 1) == 0))

  pure_lib_clean <- subset_library(pure_library, metab_to_keep)
  return(pure_lib_clean)
}


## Extend each peak in signal of nb_points_shift
signal_with_shift <- function(signal, nb_points_shift) {

  #Get indexes with a signal
  idx_signal <- which(signal == 1)

  #Get extremities of peaks
  min_extremities <- idx_signal[!((idx_signal - 1) %in% idx_signal)]
  max_extremities <- idx_signal[!((idx_signal + 1) %in% idx_signal)]

  #New signal indexes with shift
  shift_before <- unlist(lapply(min_extremities,
                                function(x) max(x - nb_points_shift, 1):x))
  shift_after <- unlist(lapply(max_extremities,
                               function(x) x:min(x + nb_points_shift,
                                                 length(signal))))

  #New signal
  signal[c(shift_before, shift_after)] <- 1

  return(signal)
}


## Non-negative least square (formula : y~x with weights w)
lm_constrained <- function(y, x, w = 1, precision = 1){
  #constraints
  Amat <- diag(rep(1, ncol(x)))
  b0 <- rep(0, ncol(x))

  #quadratic function to be minimized
  D <- t(x) %*% (w * x) * precision
  dlittle <- t(x) %*% (w * y) * precision

  #optimisation
  res.lm <- quadprog::solve.QP(D, dlittle, Amat, b0)

  #coefficients
  coefficients <- res.lm$solution
  coefficients[coefficients < 0] <- 0

  #residuals
  residuals <- y - x %*% coefficients

  return(list(coefficients = coefficients, residuals = residuals))
}


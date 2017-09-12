## Compute the area under the curve (y = f(x)) using the trapezoidal rule
AUC <- function(x, y) {
  auc <- sum(diff(x) * rollmean(y, 2))
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





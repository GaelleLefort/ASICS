## Find the best translation between each pure spectra and mixture ####
#and sort metabolites by regression residuals
#' @importFrom stats lm
#' @importFrom speaq findShiftStepFFT
#' @importFrom plyr laply
#' @keywords internal
translate_library <- function(cleaned_spectrum, cleaned_library,
                              mixture_weights, nb_points_shift){

  #find best shift with FFT cross-correlation
  ref <- as.numeric(cleaned_spectrum@spectra)*mixture_weights
  which_shift <-
    apply(cleaned_library@spectra, 2,
          function(x) findShiftStepFFT(ref, x,
                                       maxShift = nb_points_shift)$stepAdj)
  shift <- which_shift * (cleaned_library@ppm.grid[2] -
                            cleaned_library@ppm.grid[1])


  #create a matrix of all pure spectra with shift
  spectra_all_shift <- rbind(matrix(0, nrow = nb_points_shift,
                                    ncol = ncol(cleaned_library@spectra)),
                             cleaned_library@spectra,
                             matrix(0, nrow = nb_points_shift,
                                    ncol = ncol(cleaned_library@spectra)))

  #shift library according to FFT cross-correlation
  shifted_library <- cleaned_library
  shifted_library@spectra <-
    t(laply(as.list(1:ncol(cleaned_library@spectra)),
            function(i)
              spectra_all_shift[(nb_points_shift - which_shift[i] + 1):
                                  (nb_points_shift - which_shift[i] +
                                     nrow(shifted_library@spectra)), i]))

  #optimal residuals
  residuals_opti <-
    unlist(lapply(1:length(shifted_library),
                  function(i) sum((lm(cleaned_spectrum@spectra~
                                        shifted_library@spectra[, i] -
                                        1)$residuals) ^ 2)))

  #metabolites sorted by decreasing residual sum
  order <- sort(residuals_opti, decreasing = FALSE, index.return = TRUE)$ix
  shift <- shift[order]
  sorted_library <- shifted_library[order]

  return(list(sorted_library = sorted_library, shift = shift))
}




##Localized deformations of pure spectra
#' @importFrom methods is
#' @importFrom plyr laply
#' @keywords internal
deform_library <- function(cleaned_spectrum, sorted_library,
                           mixture_weights, nb_points_shift,
                           max.shift, shift){

  #Shifted library
  deformed_library <- sorted_library

  #Algorithm parameters
  nb_iter_by_library <- 1

  nb_iter_lib <- 0

  while(nb_iter_lib < nb_iter_by_library){
    #Linear regression between mixture and each pure spectra
    least_square <- try(lm_constrained(as.numeric(cleaned_spectrum@spectra),
                                       deformed_library@spectra,
                                       mixture_weights), silent = TRUE)

    if(is(least_square, "try-error")){
      least_square <- lm_constrained(as.numeric(cleaned_spectrum@spectra),
                                     deformed_library@spectra,
                                     mixture_weights, 10e-3)
    }

    #Deform each spectrum
    deformed_library@spectra <-
      t(laply(as.list(1:length(deformed_library)),
              deform_spectra, deformed_library, least_square,
              mixture_weights, nb_points_shift, max.shift, shift))

    nb_iter_lib <- nb_iter_lib + 1

  }


  return(deformed_library)
}




## Deforme each peak of a pure spectrum to align it on the complex mixture
deform_spectra <- function(idx_to_deform, pure_lib, least_square,
                           mixture_weights, nb_points_shift, max.shift, shift) {


  #Algorithm parameters
  peak_threshold <- 1
  nb_iter_by_peak <- 2

  #Linear regression coefficient of spectrum idx_to_deform
  LS_coeff <- max(0, as.numeric(least_square$coefficients[idx_to_deform]))
  LS_residuals <- as.numeric(least_square$residuals)

  #Expanded connected components
  signal_peak_lib <- which(pure_lib@spectra[, idx_to_deform] > peak_threshold)

  min_extremities <- signal_peak_lib[!((signal_peak_lib - 1) %in%
                                         signal_peak_lib)] -
    floor(nb_points_shift / 10)
  max_extremities <- signal_peak_lib[!((signal_peak_lib + 1) %in%
                                         signal_peak_lib)] +
    floor(nb_points_shift / 10)
  peaks_extremities <- cbind(min_extremities, max_extremities)

  #Remove overlapping
  long_signal <- unique(unlist(apply(peaks_extremities, 1,
                                     function(x) x[[1]]:x[[2]])))

  min_extremities <- long_signal[!((long_signal - 1) %in% long_signal)]
  max_extremities <- long_signal[!((long_signal + 1) %in% long_signal)]
  peaks_extremities <- cbind(min_extremities, max_extremities)

  #Deform on each connected component
  for(peak in 1:nrow(peaks_extremities)){

    #area of peak to deforme
    peak_area <- peaks_extremities[peak, 1]:peaks_extremities[peak, 2]
    #peak to deform
    to_deform <- pure_lib@spectra[peak_area, idx_to_deform]
    #grid to deform
    grid_to_deform <- pure_lib@ppm.grid[peak_area]

    #parameter of deformation
    range_a <- seq(- (max.shift / 5) / (max(grid_to_deform) -
                                          min(grid_to_deform)) / 0.5^2,
                   (max.shift / 5) / (max(grid_to_deform) -
                                        min(grid_to_deform)) / 0.5^2,
                   length.out = 20)
    range_a <- range_a[range_a < 1 & range_a > -1]

    #residuals without those corresponding to the current spectrum
    residuals_without_idx <- LS_residuals[peak_area] + LS_coeff * to_deform

    #repliacate nb_iter_by_peak times
    iter <- 0
    while(iter < nb_iter_by_peak){
      #Optimisation criterion
      opti_criterion <- abs(sum(to_deform * residuals_without_idx *
                                  mixture_weights[peak_area]) /
                              sqrt(sum((to_deform ^ 2) *
                                         mixture_weights[peak_area])))

      for(a in range_a){
        #for a shift of a:
        deformed_spectrum_small <- deforme(grid_to_deform, to_deform, a)
        deformed_grid <- phi(grid_to_deform, a)

        new_opti <- abs(sum(deformed_spectrum_small * residuals_without_idx *
                              mixture_weights[peak_area]) /
                          sqrt(sum((deformed_spectrum_small ^ 2) *
                                     mixture_weights[peak_area])))

        if(new_opti > opti_criterion &
           max(abs(deformed_grid - pure_lib@ppm.grid[peak_area] +
                   shift[idx_to_deform])) < max.shift){
          to_deform <- deformed_spectrum_small
          opti_criterion <- new_opti
        }
      }
      iter <- iter + 1
    }
    pure_lib@spectra[peak_area, idx_to_deform] <- to_deform
  }

  #Normalize deformed spectrum
  pure_lib@spectra[, idx_to_deform] <- pure_lib@spectra[, idx_to_deform] /
    AUC(pure_lib@ppm.grid, pure_lib@spectra[, idx_to_deform])



  return(pure_lib@spectra[, idx_to_deform])
}


## Deforme grid x with parameter a and compute new spectrum on deformed grid
#(y = old spectrum)
deforme <- function(x, y, a) {
  phix <- phi(x, a)
  return(change_grid(y, x, phix))
}

## Deforme grid x with parameter a
phi <- function(x, a) {
  # Put the old grid on [0,1]
  u <- min(x)
  v <- max(x) - min(x)
  z <- (x - u)/v

  # New grid with deformation : phi(x) = ax(1-x) + x
  tt <- z + a*z*(1 - z)
  return(u + tt*v)
}





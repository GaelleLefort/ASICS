## Find the best translation between each pure spectra and mixture ####
#and sort metabolites by regression residuals
#' @importFrom stats lm
#' @keywords internal
translate_library <- function(mixture, pure_lib_clean, nb_points_shift){
  #create a matrix with all possible shifted mixture in column
  mixture_all_shift <- c(rep(0, nb_points_shift),
                         mixture,
                         rep(0, nb_points_shift))

  mixture_shift <- t(plyr::laply(.data = as.list(1:(nb_points_shift * 2 + 1)),
                                 .fun = function (x)
                                   mixture_all_shift[x:(x + length(mixture) -
                                                          1)]))

  #create a matrix of all pure spectra with shift
  spectra_all_shift <- rbind(matrix(0, nrow = nb_points_shift,
                                    ncol = ncol(pure_lib_clean$spectra)),
                             pure_lib_clean$spectra,
                             matrix(0, nrow = nb_points_shift,
                                    ncol = ncol(pure_lib_clean$spectra)))

  #compute least square between all possible shifted mixture and all metabolites
  XtX <- diag(colSums(pure_lib_clean$spectra^2))
  least_square_shift <- solve(XtX) %*% t(pure_lib_clean$spectra) %*%
    mixture_shift
  max_least_square <- apply(least_square_shift, 1, which.max)

  #shift between original spectrum and translate spectrum
  shift <- (max_least_square - (ncol(mixture_shift) - 1) / 2) *
    diff(pure_lib_clean$grid)[1]

  #shift library according to least square
  pure_lib_shifted <- pure_lib_clean
  pure_lib_shifted$spectra <-
    t(plyr::laply(as.list(1:ncol(pure_lib_clean$spectra)),
                  function(i)
                    spectra_all_shift[(ncol(mixture_shift) - max_least_square[i]
                                       + 1):(ncol(mixture_shift) -
                                               max_least_square[i] +
                                           nrow(pure_lib_shifted$spectra)), i]))

  #optimal residuals
  residuals_opti <- unlist(lapply(1:length(pure_lib_shifted$name),
                                  function(x)
                                    sum((lm(mixture_shift[, max_least_square[x]]~
                                              pure_lib_shifted$spectra[,x] -
                                              1)$residuals) ^ 2)))

  #metabolites sorted by decreasing residual sum
  order <- sort(residuals_opti, decreasing = FALSE, index.return = TRUE)$ix
  shift <- shift[order]
  pure_lib_sorted <- subset_library(pure_lib_shifted, order)

  return(list(pure_lib_sorted = pure_lib_sorted, shift = shift))
}



##Localized deformations of pure spectra
#' @importFrom methods is
#' @keywords internal
deform_library <- function(mixture, pure_lib_sorted, nb_points_shift,
                           max.shift, shift){
  #noises and weights
  s1 <- 0.172 #standard deviation of multiplicative noise
  s2 <- 0.15 #standard deviation of additive noise
  noises <- abs(mixture) * s1 ^ 2 + s2 ^ 2
  mixture_weights <- 1 / noises

  #Shifted library
  pure_lib_deformed <- pure_lib_sorted

  #Algorithm parameters
  nb_iter_by_library <- 5

  nb_iter_lib <- 0

  while(nb_iter_lib < nb_iter_by_library){

    #Linear regression between mixture and each pure spectra
    least_square <- try(lm_constrained(mixture, pure_lib_deformed$spectra,
                                       mixture_weights), silent = TRUE)

    if(is(least_square, "try-error")){
      least_square <- lm_constrained(mixture, pure_lib_deformed$spectra,
                                         mixture_weights, 10e-3)
    }

    #Deform each spectrum
    pure_lib_deformed$spectra <-
      t(plyr::laply(as.list(1:ncol(pure_lib_deformed$spectra)),
                    deform_spectra, pure_lib_deformed, least_square,
                    mixture_weights, nb_points_shift, max.shift, shift))

    nb_iter_lib <- nb_iter_lib + 1
  }

  return(pure_lib_deformed)
}



## Deforme each peak of a pure spectrum to align it on the complex mixture
deform_spectra <- function(idx_to_deform, pure_lib, least_square,
                           mixture_weights, nb_points_shift, max.shift, shift) {
  #Algorithm parameters
  peak_threshold <- 1
  nb_iter_by_peak <- 4

  #Linear regression coefficient of spectrum idx_to_deform
  LS_coeff <- max(0, as.numeric(least_square$coefficients[idx_to_deform]))
  LS_residuals <- as.numeric(least_square$residuals)

  #Expanded connected components
  signal_peak_lib <- which(pure_lib$spectra[, idx_to_deform] > peak_threshold)

  min_extremities <- signal_peak_lib[!((signal_peak_lib - 1) %in%
                                         signal_peak_lib)] -
    floor(nb_points_shift / 5)
  max_extremities <- signal_peak_lib[!((signal_peak_lib + 1) %in%
                                         signal_peak_lib)] +
    floor(nb_points_shift / 5)
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
    to_deform <- pure_lib$spectra[peak_area, idx_to_deform]
    #grid to deform
    grid_to_deform <- pure_lib$grid[peak_area]

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
           max(abs(deformed_grid - pure_lib$grid[peak_area])) +
           shift[idx_to_deform] < max.shift){
          to_deform <- deformed_spectrum_small
          opti_criterion <- new_opti
        }
      }
      iter <- iter + 1
    }
    pure_lib$spectra[peak_area, idx_to_deform] <- to_deform
  }

  #Normalize deformed spectrum
  pure_lib$spectra[, idx_to_deform] <- pure_lib$spectra[, idx_to_deform] /
    AUC(pure_lib$grid, pure_lib$spectra[, idx_to_deform])

  return(pure_lib$spectra[, idx_to_deform])
}


## Deforme grid x with parameter a and compute new spectrum on deformed grid
#(y = old spectrum)
deforme <- function(x, y, a) {
  phix <- phi(x, a)
  return(f_o_phi(x, y, phix))
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


## Linear interpolation to adapt f on the old grid x to the new grid phix
f_o_phi <- function(x, f, phix) {

  #phix must be included in x
  if (phix[1] < x[1]) phix[1] <- x[1]
  n <- length(x)
  nphi <- length(phix)
  if (phix[nphi] > x[n]) phix[nphi] <- x[n]

  #Combine both grid (1 = old and 0 = new)
  grid_indicator <- rep(c(1, 0), c(n, nphi))
  x[1] <- x[1] - 1e-06
  x[n] <- x[n] + 1e-06
  combined_grid <- c(x, phix)

  #Sort combined grid and get points on old grid just before points on new grid
  idx_sorted_grid <- order(combined_grid)
  u <- cumsum(grid_indicator[idx_sorted_grid])
  lower_bound <- u[grid_indicator[idx_sorted_grid] == 0]

  #Spectrum on new grid
  fophi <- f[lower_bound] + (phix - x[lower_bound]) *
    (f[(lower_bound + 1)] - f[lower_bound]) /
    (x[(lower_bound + 1)] - x[lower_bound])

  if (is.na(fophi[nphi])) fophi[nphi] <- f[n]
  if (is.na(fophi[1])) fophi[1] <- f[1]

  return(fophi)
}



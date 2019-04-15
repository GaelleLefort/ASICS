## Find the best translation between each pure spectra and mixture and sort
# metabolites by regression residuals
#' @importFrom stats lm
#' @importFrom plyr laply
#' @keywords internal
.translateLibrary <- function(spectrum_obj, nb_points_shift, max.shift){

  # find best shift with FFT cross-correlation
  ref <- as.numeric(spectrum_obj[["cleaned_spectrum"]]@spectra) *
    spectrum_obj[["mixture_weights"]]
  which_shift <-
    apply(spectrum_obj[["cleaned_library"]]@spectra, 2,
          function(x) .findBestShift(ref, x,
                                       maxShift = nb_points_shift)$stepAdj)
  spectrum_obj[["shift"]] <-
    which_shift * (spectrum_obj[["cleaned_library"]]@ppm.grid[2] -
                     spectrum_obj[["cleaned_library"]]@ppm.grid[1])


  # create a matrix of all pure spectra with shift
  spectra_all_shift <-
    rbind(matrix(0, nrow = nb_points_shift,
                 ncol = ncol(spectrum_obj[["cleaned_library"]]@spectra)),
          spectrum_obj[["cleaned_library"]]@spectra,
          matrix(0, nrow = nb_points_shift,
                 ncol = ncol(spectrum_obj[["cleaned_library"]]@spectra)))

  # shift library according to FFT cross-correlation
  spectrum_obj[["cleaned_library"]]@spectra <-
    Matrix(t(laply(as.list(seq_len(ncol(
      spectrum_obj[["cleaned_library"]]@spectra))),
      function(i)
        as.matrix(
          spectra_all_shift[(nb_points_shift - which_shift[i] + 1):
                              (nb_points_shift - which_shift[i] +
                              nrow(spectrum_obj[["cleaned_library"]]@spectra)),
                            i]))))

  # optimal residuals
  residuals_opti <-
    unlist(lapply(as.list(seq_along(spectrum_obj[["cleaned_library"]])),
                  function(i) sum((lm(
                    as.matrix(spectrum_obj[["cleaned_spectrum"]]@spectra)~
                    as.matrix(spectrum_obj[["cleaned_library"]]@spectra)[, i] -
                      1)$residuals) ^ 2)))

  # metabolites sorted by decreasing residual sum
  order <- sort(residuals_opti, decreasing = FALSE, index.return = TRUE)$ix
  spectrum_obj[["shift"]] <- spectrum_obj[["shift"]][order]
  sorted_library <- spectrum_obj[["cleaned_library"]][order]

  spectrum_obj[["max_shift"]] <-
    rep(max.shift, length(spectrum_obj[["cleaned_library"]]))
  spectrum_obj[["max_nb_points_shift"]] <-
    rep(nb_points_shift, length(spectrum_obj[["cleaned_library"]]))

  return(spectrum_obj)
}


## Find the best translation between each pure spectra and mixture (test of
## various max.shift to find the best one) then sort
# metabolites by regression residuals
#' @importFrom stats lm
#' @importFrom plyr laply adply
#' @keywords internal
.translateLibrary_combineVersion <- function(spectra_obj, max.shift,
                                             nb_points_shift, spec_bin,
                                             pure.library, ncores, ntasks,
                                             verbose){

  # compute all shifts and relative concentrations
  if (verbose) cat("Compute shifts for all maximum shift values \n")
  all_res <- bplapply(spectra_obj, .allShift, nb_points_shift,
                      BPPARAM = .createEnv(ncores, ntasks, verbose))

  all_shift <- bplapply(seq_along(nb_points_shift), function(y)
    join_all(lapply(all_res, function(x, y) x[[y]][["shift"]], y),
             by = "metab", type = "full"),
    BPPARAM = .createEnv(ncores, ntasks, verbose))
  all_shift <- bplapply(all_shift,
                        function(x) {rownames(x) <- x$metab; x$metab <- NULL;
                        x[is.na(x)] <- 0; return(x)},
                        BPPARAM = .createEnv(ncores, ntasks, verbose))

  # same shift for all sample
  if (verbose) cat("Put the median shift for extreme shift values \n")
  all_shift <- bplapply(all_shift, .same_shift,
                        BPPARAM = .createEnv(ncores, ntasks, verbose))


  all_conc <- bplapply(seq_along(nb_points_shift), function(y)
    join_all(lapply(all_res, function(x, y) x[[y]][["relative_conc"]], y),
             by = "metab", type = "full"),
    BPPARAM = .createEnv(ncores, ntasks, verbose))
  all_conc <- bplapply(all_conc,
                        function(x) {rownames(x) <- x$metab; x$metab <- NULL;
                        x[is.na(x)] <- 0; return(x)},
                       BPPARAM = .createEnv(ncores, ntasks, verbose))

  # compute correlation
  if (verbose)
    cat("Compute correlations between buckets and quantifications \n")
  list_all_cor <- bplapply(all_conc, .compute_correlation, spec_bin,
                           pure.library, nb_points_shift,
                           BPPARAM = .createEnv(ncores, ntasks, verbose))
  list_all_cor <- lapply(seq_along(list_all_cor),
                         function(x) {colnames(list_all_cor[[x]])[2] <- x;
                         return(list_all_cor[[x]])})
  all_cor <- join_all(list_all_cor, by = "metab", type = "full")

  # best shift by metabolite and sample
  all_cor <- adply(all_cor, 1,
                   function(x) {if (all(is.na(x[2:length(x)]))) x[2] <- 1;
                   return(x)})
  which_corr_max <- apply(all_cor[, 2:ncol(all_cor)], 1, which.max)
  names(which_corr_max) <- all_cor[, 1]

  final_shift <-
    do.call("rbind", lapply(seq_along(which_corr_max),
                            function(i) all_shift[[which_corr_max[i]]][i,]))

  # shift all spectra according to final shift
  if (verbose) cat("Shift all spectra according to the best shift \n")
  spectra_obj <- bplapply(spectra_obj, .do_shift_corr, final_shift, max.shift,
                          nb_points_shift, which_corr_max,
                          BPPARAM = .createEnv(ncores, ntasks, verbose))

  return(spectra_obj)
}

# Compute all shifs for all maximum shifts and all spectra
.allShift <- function(spectrum_obj, nb_points_shift) {
  ## reference spectrum
  ref <- as.numeric(spectrum_obj[["cleaned_spectrum"]]@spectra) *
    spectrum_obj[["mixture_weights"]]

  # Find the best translation between each pure spectra and mixture for each
  #max shift
  shift_and_conc <- lapply(nb_points_shift, .allShift_internal, spectrum_obj,
                           ref)

  return(shift_and_conc)
}


#' @importFrom plyr ldply
.allShift_internal <- function(nb_points_shift_i, spectrum_obj, ref) {
  # find best shift with FFT cross-correlation
  which_shift <-
    apply(spectrum_obj[["cleaned_library"]]@spectra, 2,
          function(x) .findBestShift(ref, x,
                                     maxShift = nb_points_shift_i)$stepAdj)

  # data-frame of shifts
  shift_i <- data.frame(D = which_shift)
  colnames(shift_i) <- spectrum_obj[["cleaned_spectrum"]]@sample.name
  shift_i <- cbind(metab = spectrum_obj[["cleaned_library"]]@sample.name,
                   shift_i)

  # perform shit
  spectrum_obj[["cleaned_library"]]@spectra <-
    Matrix(t(ldply(seq_len(nrow(shift_i)), function(x)
      .doShift(spectrum_obj[["cleaned_library"]]@spectra[, x],
                     shift_i[x, 2], 0, 0))))


  # pseudo MLE estimation
  B2 <- try(.lmConstrained(spectrum_obj[["cleaned_spectrum"]]@spectra,
                           spectrum_obj[["cleaned_library"]]@spectra,
                           spectrum_obj[["mixture_weights"]])$coefficients,
            silent = TRUE)
  if(is(B2, "try-error")){
    B2 <- .lmConstrained(spectrum_obj[["cleaned_spectrum"]]@spectra,
                         spectrum_obj[["cleaned_library"]]@spectra,
                         spectrum_obj[["mixture_weights"]], 10e-3)$coefficients
  }

  # data-frame of relative concentrations
  relative_conc_i <-
    data.frame(RC = B2 / spectrum_obj[["cleaned_library"]]@nb.protons)
  colnames(relative_conc_i) <- spectrum_obj[["cleaned_spectrum"]]@sample.name
  relative_conc_i <-
    cbind(metab = spectrum_obj[["cleaned_library"]]@sample.name,
          relative_conc_i)

  return(list(shift = shift_i, relative_conc = relative_conc_i))
}


# Put the median shift for extreme shift values
.same_shift <- function(shift) {
  shift_med <- round(apply(shift, 1, median), 0)

  change_noncommom_shift <-
    as.data.frame(t(vapply(seq_len(nrow(shift)),
           function(i) unlist(ifelse(shift[i, ] > shift_med[i] - 5 &
                                       shift[i, ] < shift_med[i] + 5,
                                     shift[i, ], shift_med[i])),
           FUN.VALUE = numeric(ncol(shift)))))
  rownames(change_noncommom_shift) <- rownames(shift)
  colnames(change_noncommom_shift) <- colnames(shift)

  return(change_noncommom_shift)
}


# Compute the correlation with buckets
#' @importFrom stats cor
.compute_correlation <- function(concentration_i, spec_bin,
                                 pure.library, nb_points_shift) {
  corr_i <- suppressWarnings(cor(t(concentration_i),
                                 t(spec_bin), use = "complete.obs"))

  corr_max <- unlist(lapply(seq_len(nrow(corr_i)), .correlation_buckets, corr_i,
                            pure.library, nb_points_shift))

  return(data.frame(metab = rownames(corr_i), corr_max))
}


# Find buckets near pure spectra peaks and compute the correlation
.correlation_buckets <- function(idx, corr_i, pure.library, nb_points_shift) {
  pure_spectrum <- pure.library@spectra[, pure.library@sample.name ==
                                          rownames(corr_i)[idx]]

  # expanded connected components
  signal_peak_lib <- which(pure_spectrum != 0)

  min_extremities <- signal_peak_lib[!((signal_peak_lib - 1) %in%
                                         signal_peak_lib)] -
    floor(nb_points_shift[length(nb_points_shift)])
  max_extremities <- signal_peak_lib[!((signal_peak_lib + 1) %in%
                                         signal_peak_lib)] +
    floor(nb_points_shift[length(nb_points_shift)])
  peaks_extremities <-
    cbind(vapply(min_extremities, max, 1,
                 FUN.VALUE = numeric(1)),
          vapply(max_extremities, min, length(pure.library@ppm.grid),
                 FUN.VALUE = numeric(1)))

  # remove overlapping
  long_signal <- unique(unlist(apply(peaks_extremities, 1,
                                     function(x) x[[1]]:x[[2]])))

  min_extremities <- long_signal[!((long_signal - 1) %in% long_signal)]
  max_extremities <- long_signal[!((long_signal + 1) %in% long_signal)]
  peaks_extremities <- cbind(pure.library@ppm.grid[min_extremities],
                             pure.library@ppm.grid[max_extremities])

  # which buckets have signal in pure spectrum
  ppm_bucket <- as.numeric(names(corr_i[idx, ]))
  which_buck <- unlist(apply(peaks_extremities, 1,
                             function(x) which(ppm_bucket > x[1] &
                                                 ppm_bucket < x[2])))

  # return the maximum correlation
  return(max(corr_i[idx, which_buck]))
}


# Perform the best library shift for each spectra
.do_shift_corr <- function(spectrum_obj, final_shift, max.shift,
                           nb_points_shift, which_corr_max) {
  obj <- spectrum_obj[["cleaned_spectrum"]]
  lib <- spectrum_obj[["cleaned_library"]]
  spectrum_obj[["cleaned_library"]]@spectra <-
    Matrix(do.call("cbind",
            lapply(seq_len(ncol(lib@spectra)),
                   function(idx)
                     .doShift(lib@spectra[, idx],
                                    final_shift[lib@sample.name[idx],
                                                obj@sample.name], 0, 0))))

  spectrum_obj[["shift"]] <- unlist(
    lapply(seq_len(ncol(lib@spectra)),
                   function(idx)
                     final_shift[lib@sample.name[idx], obj@sample.name] *
             (lib@ppm.grid[2] - lib@ppm.grid[1])))

  spectrum_obj[["max_shift"]] <- unlist(
    lapply(seq_len(ncol(lib@spectra)),
           function(idx)
             max.shift[which_corr_max[lib@sample.name[idx]]]))

  spectrum_obj[["max_nb_points_shift"]] <- unlist(
    lapply(seq_len(ncol(lib@spectra)),
           function(idx)
             nb_points_shift[which_corr_max[lib@sample.name[idx]]]))

  return(spectrum_obj)
}


## Localize deformations of pure spectra
#' @importFrom methods is
#' @importFrom plyr laply
#' @keywords internal
.deformLibrary <- function(spectrum_obj){

  # linear regression between mixture and each pure spectra
  least_square <-
    try(.lmConstrained(as.numeric(spectrum_obj[["cleaned_spectrum"]]@spectra),
                       spectrum_obj[["cleaned_library"]]@spectra,
                       spectrum_obj[["mixture_weights"]]), silent = TRUE)

  if(is(least_square, "try-error")){
    least_square <-
      .lmConstrained(as.numeric(spectrum_obj[["cleaned_spectrum"]]@spectra),
                     spectrum_obj[["cleaned_library"]]@spectra,
                     spectrum_obj[["mixture_weights"]], 10e-3)
  }

  # deform each spectrum
  spectrum_obj[["cleaned_library"]]@spectra <-
    Matrix(t(laply(as.list(seq_along(spectrum_obj[["cleaned_library"]])),
            .deformSpectra, spectrum_obj[["cleaned_library"]], least_square,
            spectrum_obj[["mixture_weights"]],
            spectrum_obj[["max_nb_points_shift"]],
            spectrum_obj[["max_shift"]],
            spectrum_obj[["shift"]])))

  return(spectrum_obj)
}




## Deforme each peak of a pure spectrum to align it on the complex mixture
.deformSpectra <- function(idx_to_deform, pure_lib, least_square,
                           mixture_weights, nb_points_shift, max.shift, shift) {

  # algorithm parameters
  peak_threshold <- 1
  nb_iter_by_peak <- 2

  # linear regression coefficient of spectrum idx_to_deform
  LS_coeff <- max(0, as.numeric(least_square$coefficients[idx_to_deform]))
  LS_residuals <- as.numeric(least_square$residuals)

  # expanded connected components
  signal_peak_lib <- which(pure_lib@spectra[, idx_to_deform] > peak_threshold)

  min_extremities <- signal_peak_lib[!((signal_peak_lib - 1) %in%
                                         signal_peak_lib)] -
    floor(nb_points_shift[idx_to_deform] / 10)
  max_extremities <- signal_peak_lib[!((signal_peak_lib + 1) %in%
                                         signal_peak_lib)] +
    floor(nb_points_shift[idx_to_deform] / 10)
  peaks_extremities <-
    cbind(vapply(min_extremities, max, 1,
                 FUN.VALUE = numeric(1)),
          vapply(max_extremities, min, length(pure_lib@ppm.grid),
                 FUN.VALUE = numeric(1)))

  # remove overlapping
  long_signal <- unique(unlist(apply(peaks_extremities, 1,
                                     function(x) x[[1]]:x[[2]])))

  min_extremities <- long_signal[!((long_signal - 1) %in% long_signal)]
  max_extremities <- long_signal[!((long_signal + 1) %in% long_signal)]
  peaks_extremities <- cbind(min_extremities, max_extremities)

  # remove small peak (less than 3 points)
  peaks_extremities <- t(peaks_extremities[peaks_extremities[, 2] -
                                            peaks_extremities[, 1] >= 3, ])

  if (nrow(peaks_extremities) > 0) {
    # deform on each connected component
    for(peak in seq_len(nrow(peaks_extremities))){

      # area of peak to deforme
      peak_area <- peaks_extremities[peak, 1]:peaks_extremities[peak, 2]
      # peak to deform
      to_deform <- pure_lib@spectra[peak_area, idx_to_deform]
      # grid to deform
      grid_to_deform <- pure_lib@ppm.grid[peak_area]

      # parameter of deformation
      range_a <- seq(-(max.shift[idx_to_deform] / 10) /
                       (max(grid_to_deform) - min(grid_to_deform)) / 0.5^2,
                     (max.shift[idx_to_deform] / 10) /
                       (max(grid_to_deform) - min(grid_to_deform)) / 0.5^2,
                     length.out = 20)
      range_a <- range_a[range_a < 1 & range_a > -1]

      # residuals without those corresponding to the current spectrum
      residuals_without_idx <- LS_residuals[peak_area] + LS_coeff * to_deform

      # repliacate nb_iter_by_peak times
      iter <- 0
      while(iter < nb_iter_by_peak){
        # optimisation criterion
        opti_criterion <- abs(sum(to_deform * residuals_without_idx *
                                    mixture_weights[peak_area]) /
                                sqrt(sum((to_deform ^ 2) *
                                           mixture_weights[peak_area])))

        for(a in range_a){
          # for a shift of a:
          deformed_spectrum_small <- .deforme(grid_to_deform, to_deform, a)
          deformed_grid <- .phi(grid_to_deform, a)

          new_opti <- abs(sum(deformed_spectrum_small * residuals_without_idx *
                                mixture_weights[peak_area]) /
                            sqrt(sum((deformed_spectrum_small ^ 2) *
                                       mixture_weights[peak_area])))

          if(new_opti > opti_criterion &
             max(abs(deformed_grid - pure_lib@ppm.grid[peak_area] +
                     shift[idx_to_deform])) < max.shift[idx_to_deform]){
            to_deform <- deformed_spectrum_small
            opti_criterion <- new_opti
          }
        }
        iter <- iter + 1
      }
      pure_lib@spectra[peak_area, idx_to_deform] <- to_deform
    }
  }


  # normalize deformed spectrum
  pure_lib@spectra[, idx_to_deform] <- pure_lib@spectra[, idx_to_deform] /
    .AUC(pure_lib@ppm.grid, pure_lib@spectra[, idx_to_deform])



  return(pure_lib@spectra[, idx_to_deform])
}


## Deform grid x with parameter a and compute new spectrum on deformed grid
#(y = old spectrum)
.deforme <- function(x, y, a) {
  phix <- .phi(x, a)
  return(.changeGrid(y, x, phix))
}

## Deforme grid x with parameter a
.phi <- function(x, a) {
  # put the old grid on [0,1]
  u <- min(x)
  v <- max(x) - min(x)
  z <- (x - u)/v

  # new grid with deformation: phi(x) = ax(1-x) + x
  tt <- z + a*z*(1 - z)
  return(u + tt*v)
}





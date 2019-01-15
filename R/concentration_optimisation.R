#### Threshold and concentration optimisation for each metabolite ####
#' @importFrom stats rnorm runif
#' @importFrom methods is
#' @keywords internal
.concentrationOpti <- function(spectrum_obj){

  # variance matrix of maximum likelihood
  A <- as.numeric(1 / sqrt(1/spectrum_obj[["mixture_weights"]])) *
    spectrum_obj[["cleaned_library"]]@spectra
  VMLE <- solve(t(A)%*%A)

  # first threshold minimisation
  N <- 1000
  C <- t(chol(VMLE))
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N, rnorm(nrow(C)*N))

  # regularization paramater of lasso under positive constraints
  se <- sqrt(diag(VMLE))

  # L1 norm of threshold optimisation
  nb_draw <- 30
  p <- ncol(VMLE)
  W <- matrix(nrow = nb_draw, ncol = p)
  u <- numeric(nb_draw)

  delta0 <- rep(0.5, p) # 0.5 by metabolite
  a_min <- sum(.computeThreshold(delta0, ZMLE, se)) #threshold sum by metabolite
  err <- 0.4

  for(i in seq_len(400)) {
    err <- 0.99 * err
    W <- matrix(delta0 + err * c(runif(p * nb_draw, -1, 1)), nrow = nb_draw,
                byrow = TRUE)
    W[W > 1]  <-  1
    W[W < 0] <- 0.0001
    u <- apply(W, 1, .toMinimize, ZMLE, se)
    if (min(u) < a_min)
    {
      delta0 <- W[which.min(u), ]
      a_min <- min(u)
    }
  }

  # pseudo MLE estimation
  B2 <- try(.lmConstrained(spectrum_obj[["cleaned_spectrum"]]@spectra,
                           spectrum_obj[["cleaned_library"]]@spectra,
                           spectrum_obj[["mixture_weights"]])$coefficients,
            silent = TRUE)
  if(is(B2, "try-error")){
    B2 <- .lmConstrained(spectrum_obj[["cleaned_spectrum"]]@spectra,
                         spectrum_obj[["cleaned_library"]]@spectra,
                         spectrum_obj[["mixture_weights"]],
                         10e-3)$coefficients
  }

  # compute all thresholds
  N <- 10000
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N, rnorm(nrow(C)*N))

  # concentration lasso estimation with positive constraints
  identified_metab <- (B2 > .tuning(delta0, ZMLE) / delta0) & (B2 > 0)

  identified_library <- spectrum_obj[["cleaned_library"]][which(identified_metab)]

  B_final_tot <- try(.lmConstrained(spectrum_obj[["cleaned_spectrum"]]@spectra,
                                    identified_library@spectra,
                                    spectrum_obj[["mixture_weights"]])$coefficients,
                     silent = TRUE)
  if(is(B_final_tot,"try-error")){
    B_final_tot <- .lmConstrained(spectrum_obj[["cleaned_spectrum"]]@spectra,
                                  identified_library@spectra,
                                  spectrum_obj[["mixture_weights"]],
                                  10e-3)$coefficients
  }

  # test of coefficients positivity
  B_final <- B_final_tot[B_final_tot > 0]
  identified_metab[identified_metab][B_final_tot <= 0] <- FALSE
  spectrum_obj[["cleaned_library"]] <- spectrum_obj[["cleaned_library"]][which(identified_metab)]

  # reconstituted mixture with estimated coefficents
  spectrum_obj[["est_mixture"]] <-
    spectrum_obj[["cleaned_library"]]@spectra %*% B_final
  spectrum_obj[["cleaned_library"]]@spectra <-
    spectrum_obj[["cleaned_library"]]@spectra %*% diag(B_final)

  # compute relative concentration of identified metabolites
  spectrum_obj[["relative_concentration"]] <-
    B_final / spectrum_obj[["cleaned_library"]]@nb.protons
  # sort library according to relative concentration
  sorted_idx <- sort(spectrum_obj[["relative_concentration"]], decreasing = TRUE,
                     index.return = TRUE)$ix
  spectrum_obj[["cleaned_library"]] <-
    spectrum_obj[["cleaned_library"]][sorted_idx]

  spectrum_obj[["relative_concentration"]] <-
    data.frame(spectrum_obj[["relative_concentration"]][sorted_idx])
  rownames(spectrum_obj[["relative_concentration"]]) <-
    spectrum_obj[["cleaned_library"]]@sample.name
  colnames(spectrum_obj[["relative_concentration"]]) <-
    spectrum_obj[["cleaned_spectrum"]]@sample.name

  # change format of pure library
  temp_df <- as.data.frame(as.matrix(spectrum_obj[["cleaned_library"]]@spectra))
  rownames(temp_df) <- spectrum_obj[["cleaned_library"]]@ppm.grid
  colnames(temp_df) <- spectrum_obj[["cleaned_library"]]@sample.name
  pure_lib_format <-
    reshape(temp_df, idvar = "ppm_grid", ids = row.names(temp_df),
            times = names(temp_df), timevar = "metabolite_name",
            varying = list(names(temp_df)), direction = "long",
            v.names = "intensity")
  rownames(pure_lib_format) <- NULL
  pure_lib_format <- pure_lib_format[pure_lib_format$intensity != 0, ]
  spectrum_obj[["format_library"]] <-
    cbind(sample = rep(spectrum_obj[["cleaned_spectrum"]]@sample.name,
                       nrow(pure_lib_format)), pure_lib_format)

  return(spectrum_obj)
}



##### Threshold optimisation functions
.toMinimize <- function(x, ZMLE, se) {
  sum(.computeThreshold(x, ZMLE, se))
}

#' @importFrom stats qnorm
#' @keywords internal
.computeThreshold <- function(x, ZMLE, se) {
  return(.tuning(x, ZMLE)/x + se*qnorm(0.95))
}

#' @importFrom stats quantile
#' @keywords internal
.tuning <- function(delta, ZMLE) {
  observation <- apply(as.double(delta) * (ZMLE), 2, max)
  return(quantile(observation, 0.95))
}

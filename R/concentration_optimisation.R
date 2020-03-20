#### Threshold and concentration optimisation for each metabolite ####
#' @importFrom methods is
#' @keywords internal
.concentrationOpti <- function(spectrum_obj){

  # variance matrix of maximum likelihood
  A <- as.numeric(1 / sqrt(1/spectrum_obj[["mixture_weights"]])) *
    spectrum_obj[["cleaned_library"]]@spectra
  VMLE <- as.matrix(solve(t(A)%*%A))

  th <- optimal_thres_uni(VMLE)

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

  # concentration lasso estimation with positive constraints
  identified_metab <- (B2 > th) & (B2 > 0)

  identified_library <-
    spectrum_obj[["cleaned_library"]][which(identified_metab)]

  B_final_tot <-
    try(.lmConstrained(spectrum_obj[["cleaned_spectrum"]]@spectra,
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
  spectrum_obj[["cleaned_library"]] <-
    spectrum_obj[["cleaned_library"]][which(identified_metab)]

  # reconstituted mixture with estimated coefficents
  spectrum_obj[["est_mixture"]] <-
    spectrum_obj[["cleaned_library"]]@spectra %*% B_final
  spectrum_obj[["cleaned_library"]]@spectra <-
    spectrum_obj[["cleaned_library"]]@spectra %*% diag(B_final)

  # compute relative concentration of identified metabolites
  spectrum_obj[["relative_concentration"]] <-
    B_final / spectrum_obj[["cleaned_library"]]@nb.protons
  # sort library according to relative concentration
  sorted_idx <- sort(spectrum_obj[["relative_concentration"]],
                     decreasing = TRUE, index.return = TRUE)$ix
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
#' @importFrom mvtnorm rmvnorm
optimal_thres_uni <- function(V){

  n <- 1000 # number of simulation for the quantile computation
  p <- ncol(V)

  se <- sqrt(diag(V)) # standard deviation
  C <- diag(1/se)%*%V%*%diag(1/se) # correlation matrix associated with V

  Y <- rmvnorm(n, mean = rep(0, nrow(C)), sigma = C)

  seuil <- newton_volume(C, Y) # numerical computation of optimal thresholds
  return(se * seuil)
}

newton_volume <- function(C, Y){

  # initial threshold
  sup <- apply(Y, 1, max)
  q95 <- quantile(sup, probs = 0.95, names = FALSE)
  lim1 <- rep(q95, nrow(C))

  lim_temp <- lim1
  flag <- TRUE
  while (flag == T) {
    grad <- compute_gradient(lim1, C, Y)
    inc <- 1 / grad / max(1 / grad)

    flag2 <- TRUE
    lambda <- 0.4
    while (flag2 == T) {
      lim_temp <- shrink(lim1 - lambda * inc, Y)
      if (sum(log(lim_temp)) < sum(log(lim1))) {
        lim1 <- lim_temp
        flag2 <- FALSE
      } else {
        lambda <- lambda * 0.6666
        if (lambda < 0.005) {
          flag2 <- FALSE
          flag <- FALSE
        }
      }
    }
  }

  return(lim1)
}

#' @importFrom stats dnorm
compute_gradient <- function(lim1, C, Y) {

  grad <- numeric(nrow(C))

  for (i in 1:nrow(C))  {
    Sig22_ <- solve(C[-i,-i])
    s2 <- C[i, i] - C[i, -i] %*% Sig22_ %*% C[-i,i]
    sel <- which(rowSums((Y[, -i] > lim1[-i])) == 0)
    m <- t(C[i, -i] %*% Sig22_ %*% t(Y[sel, -i]))
    grad[i] <- mean(dnorm(as.numeric(Y[sel,i]), as.numeric(m),
                          rep(s2, length(m))))
  }

  return(grad)
}

#' @importFrom stats quantile
shrink <- function(lim, Y) {
  ZZ <- Y * matrix(rep(1 / lim, nrow(Y)), nrow = nrow(Y), byrow = TRUE)
  sup <- apply(ZZ, 1, max)
  q95 <- quantile(sup, probs = 0.95, names = FALSE)
  return(lim * rep(q95, ncol(Y)))
}












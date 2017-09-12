## Deforme each peak of a pure spectrum to align it on the complex mixture
deform_spectra <- function(idx_to_deform, pure_lib, least_square,
                           mixture_weights, nb_points_shift, max.shift) {
  #Algorithm parameters
  peak_threshold <- 1
  nb_iter_deform <- 4
  range_a <- -9:9/10

  #Deformed spetrum
  deform_spectrum <- pure_lib$spectra[, idx_to_deform]

  #Expanded connected components
  signal_peak_lib <- which(pure_lib$spectra[, idx_to_deform] > peak_threshold)

  min_extremities <- signal_peak_lib[!((signal_peak_lib - 1) %in%
                                         signal_peak_lib)] - nb_points_shift
  max_extremities <- signal_peak_lib[!((signal_peak_lib + 1) %in%
                                         signal_peak_lib)] + nb_points_shift
  peaks_extremities <- cbind(min_extremities, max_extremities)

  #Remove overlapping
  long_signal <- unique(unlist(apply(peaks_extremities, 1,
                                     function(x) x[[1]]:x[[2]])))

  min_extremities <- long_signal[!((long_signal - 1) %in% long_signal)]
  max_extremities <- long_signal[!((long_signal + 1) %in% long_signal)]
  peaks_extremities <- cbind(min_extremities, max_extremities)

  #Residuals without those corresponding to the current spectrum
  residuals_without_idx <- as.numeric(least_square$residuals) +
    max(0, as.numeric(least_square$coefficients[idx_to_deform])) *
    pure_lib$spectra[, idx_to_deform]

  #Deform on each connected component
  for(peak in 1:nrow(peaks_extremities)){
    #area of peak to deforme
    peak_area <- peaks_extremities[peak, 1]:peaks_extremities[peak, 2]
    #peak to deform
    to_deform <- deform_spectrum[peak_area]
    #grid to deform
    grid_to_deform <- pure_lib$grid[peak_area]

    #repliacate nb_iter_deform times
    iter <- 0
    while(iter < nb_iter_deform){
      #Optimisation criterion
      opti_criterion <- abs(sum(to_deform * residuals_without_idx[peak_area] *
                                  mixture_weights[peak_area]) /
                              sqrt(sum((to_deform ^ 2) *
                                         mixture_weights[peak_area])))

      for(a in range_a){
        #for a shift of a:
        deformed_grid <- phi(grid_to_deform, a)
        deformed_spectrum_small <- deforme(pure_lib$grid[peak_area],
                                           to_deform, a)

        new_opti <- abs(sum(deformed_spectrum_small * residuals_without_idx[peak_area] *
                              mixture_weights[peak_area]) /
                          sqrt(sum((deformed_spectrum_small ^ 2) *
                                     mixture_weights[peak_area])))

        if(new_opti > opti_criterion &
           max(abs(deformed_grid - pure_lib$grid[peak_area])) < max.shift){
          to_deform <- deformed_spectrum_small
          opti_criterion <- new_opti
          grid_to_deform <- deformed_grid
        }
      }
      iter <- iter + 1
    }

    deform_spectrum[peak_area] <- to_deform
  }
  deform_spectrum <- deform_spectrum/AUC(pure_lib$grid, deform_spectrum)

  return(deform_spectrum)
}









f_o_phi <- function(x, f, phix) {
  #il y a des pbs avec les ex aequo dans u et v
  # Attention, il faut que inf phi(x) >= inf x ET sup phi(x)<= sup x
  #phix contient le vecteur phi(x), f contient f(x)
  #renvoie f(phi(x))

  #Le minimum et le maximum de phix doivent être compris dans x
  if (phix[1] < x[1]) phix[1] <- x[1]
  n <- length(x)
  nphi <- length(phix)
  if (phix[nphi] > x[n]) phix[nphi] <- x[n]

  #1 pour la grille d'origine et 0 pour la nouvelle
  orde <- rep(c(1, 0), c(n, nphi))

  x[1] <- x[1] - 1e-06
  x[n] <- x[n] + 1e-06

  #Concaténation des 2 grilles
  zx <- c(x, phix)

  #On ordonne les 2 grilles et on récupère les indices
  o <- order(zx)

  #On augmente de 1 si on est sur la grille d'origine
  #On a donc un vecteur d'indice allant de 1 à n où la valeur ne change pas quand
  #un point de la nouvelle grille est compris dans l'ancienne
  u <- cumsum(orde[o])

  #On récupère les points de l'ancienne grille juste inférieurs à ceux de la nouvelle
  v <- u[orde[o] == 0]

  #Nouvelle fonction par interpolation linéaire
  fophi <- f[v] + (phix - x[v])*(f[(v + 1)] - f[v])/(x[(v + 1)] - x[v])

  if (is.na(fophi[nphi])) fophi[nphi] <- f[n]
  if (is.na(fophi[1])) fophi[1] <- f[1]

  return(fophi)
}



phi <- function(x, a) {
  u <- min(x)
  v <- max(x) - min(x)
  z <- (x - u)/v
  tt <- z + a*z*(1 - z)
  return(u + tt*v)
}

deforme <- function(x, y, a) {
  phix <- phi(x, a)
  return(f_o_phi(x, y, phix)) # phix deformation de l'axe des absisses et f_o_phi est la valeur de la spectre deforme sur les points
  #de discretisation de l'axe ds absisses.
}






##### fonctions optimisation
Tuning <- function(delta, ZMLE) {
  Observation <- apply(as.double(delta)*(ZMLE), 2, max)
  return(quantile(Observation, 0.95))
}

seuil <- function(x, ZMLE, se) {
  return(Tuning(x, ZMLE)/x + se*qnorm(0.95))
}


aminimiser <- function(x, ZMLE, se) {
  sum(seuil(x, ZMLE, se))
}










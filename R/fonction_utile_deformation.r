recupere_intervalle <- function(X, e)
{
  #X est le spectre d'un m?tabolite
  # e permet de d?finir la notion de support (qui n'est pas d?ini pour le spectre bruit? d'un m?tabolite)
  supp <- which(X > e) #points du support
  u <- which((supp[-1] - supp[1:(length(supp) - 1)] - 1) != 0)#recupere les composantes connexes du support
  nb_int <- length(u) + 1
  v <- matrix(nrow = nb_int, ncol = 2)
  v[, 2] <- c(supp[u], max(supp))
  v[, 1] <- c(min(supp), supp[u + 1])
  return(v)
}


grossir_ensemble <- function(S, e)
{
  #Le vecteur S est un vecteur bool?ens dont on suppose que les extr?mit? ont suffisament de z?ros.
  #l'entier e donne dle grossissement de l'ensemble S
  M <- matrix(nrow = (2*e+1), ncol = length(S))
  u <- S[(e + 1):(length(S) - e)]
  M[1, ] <- c(u, rep(0, (2*e))) #d?calage de S le plus ? gauche
  for(i in 1:(2*e))
  {
    M[1+i, ] <- c(rep(0, i), u, rep(0, (2*e - i))) #On d?place le vecteur S vers la gauche
  }
  R <- 1*(colSums(M) != 0) #grossisement du vecteur S
  return(R)
}



recal_1_pic <- function(X, Y, d1, d2)
{
  #d1=d?calage vers la gauche et d2 d?calage vers la droite
  M <- matrix(nrow = d1 + d2 + 1, ncol = length(X))
  u <- X[(1 + d1):(length(X) - d2)]
  M[1, ] <- c(u, rep(0, (d1 + d2))) #D?calage du pic le plus ? gauche
  for(i in 1:(d1 + d2))
  {
    M[1+i, ] <- c(rep(0, i), u, rep(0, (d1 + d2 - i))) #On d?place le vecteur X vers la gauche
  }
  i0 <- which.max(Y%*%t(M))
  return(c(i0, M[i0, ]))
}


translation_intervalle <- function(Y, X, d, e)
{
  # X est la matrice des spectres m?tabolites
  # Y est le spectre du m?lange
  # d est la d?formation maximale
  mc <- lm(Y~X - 1)
  res <- mc$residuals
  for(i in 1:ncol(X))
  {
    Nexp <- res + mc$coefficients[i]*X[, i] #r?sidues plus partie non expliqu? par X[,i]
    v <- recupere_intervalle(X[, i], e)
    for(j in 1:nrow(v))
    {
      if(j == 1){b1 <- v[j, 1] - d} else {b1 <- max(v[j, 1] - d, mean(c(v[j - 1, 2], v[j, 1])))}
      d1 <- abs(v[j, 1] - b1)
      if(j == nrow(v)){b2 <- v[j, 2] + d} else {b2 <- min(v[j, 2] + d, mean(c(v[j, 2], v[j + 1, 1])))}
      d2 <- abs(b2 - v[j, 2])
      M <- recal_1_pic(X[b1:b2, i], Y[b1:b2], d1, d2)
      i0 <- M[1]
      X[b1:b2, i] <- M[-1]
      b <- lm(Nexp[b1:b2]~X[b1:b2, i] - 1)$coefficients[1]
      A <- rep(0, length(res))
      A[b1:b2] <- b*X[b1:b2, i]
      Nexp <- Nexp - A
    }
  }
  return(X)
}

deformation <- function(biblio_deformee,y)

  for (i in 1:ncol(biblio_deformee))
  {
    yy <- residu + coeff[i]*biblio[, i]
    adeformer <- biblio_deformee[, i]
    sel <- which(biblio_deformee[, i] > seuil1)
    u <- which((sel[-1] - sel[1:(length(sel) - 1)] - 1) != 0)
    nb_int <- length(u) + 1
    v <- matrix(nrow = nb_int, ncol = 2)
    v[, 2] <- c(sel[u], max(sel))
    v[, 1] <- c(min(sel), sel[u + 1])
    v[, 1] <- v[, 1] - h
    v[, 2] <- v[, 2] + h
    for (j in 1:nb_int)
    {
      xx <- x[v[j, 1]:v[j, 2]]
      zz <- adeformer[v[j, 1]:v[j, 2]]
      yyy <- yy[v[j, 1]:v[j, 2]]
      for (k in 1:nb_iter_deform)
      {
        opti <- abs(sum(yyy*zz)/sum(zz))
        flag <- F
        for (a in rangea)
        {
          candidat <- deforme(xx, zz, a)
          if (abs(sum(candidat*yyy)/sum(zz)) > opti)
          {
            zz <- candidat
            opti <- abs(sum(candidat*yyy)/sum(zz))
            flag <- T
          }
        }
      }
      if (flag==T) biblio_deformee[v[j, 1]:v[j, 2], i] <- zz
      else k <- nb_iter_deform + 1
    }
    #cat(" i =",i,"\n")
    biblio_deformee[, i] <- biblio_deformee[, i]/AUC(x, biblio_deformee[, i])
  }


phi <- function(x, a)
{
  set.seed(12345)
  u <- min(x)
  v <- max(x) - min(x)
  z <- (x - u)/v
  tt <- z + a*z*(1 - z)
  return(u + tt*v)
}

deforme <- function(x, y, a)
{
  set.seed(12345)
  phix <- phi(x, a)
  return(f_o_phi(x, y, phix)) # phix deformation de l'axe des absisses et f_o_phi est la valeur de la spectre deforme sur les points
  #de discretisation de l'axe ds absisses.
}


deform_spectra <- function(idx_to_deform, pure_lib, least_square,
                           mixture_weights, nb_points_shift, max.shift){
  #Algorithm parameters
  peak_threshold <- 1
  nb_iter_deform <- 4
  range_a <- -9:9/10
  set.seed(12345)

  #Deformed spetrum
  deform_spectrum <- pure_lib$spectra[, idx_to_deform]

  #Expanded connected components
  signal_peak_lib <- which(pure_lib$spectra[, idx_to_deform] > peak_threshold)

  min_extremities <- signal_peak_lib[!((signal_peak_lib - 1) %in%
                                         signal_peak_lib)] - nb_points_shift
  max_extremities <- signal_peak_lib[!((signal_peak_lib + 1) %in%
                                         signal_peak_lib)] + nb_points_shift
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
    #print(deform_spectrum[peak_area])
  }
  deform_spectrum <- deform_spectrum/AUC(pure_lib$grid, deform_spectrum)

  return(deform_spectrum)
}



##### fonctions optimisation
Tuning <- function(delta, ZMLE)
{
  Observation <- apply(as.double(delta)*(ZMLE), 2, max)
  return(quantile(Observation, 0.95))
}

seuil <- function(x, ZMLE, se)
{
  return(Tuning(x, ZMLE)/x + se*qnorm(0.95))
}


aminimiser <- function(x, ZMLE, se)
{
  sum(seuil(x, ZMLE, se))
}










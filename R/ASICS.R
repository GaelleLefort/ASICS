#' Automatic Statistical Identification in Complex Spectra
#'
#' Description
#' @param path folder path of the Bruker files
#' @param exclusion.areas exclusion areas of the quantification
#' @param max.shift maximum chemical shift allowed (in ppm)
#' @param which.spectra if more than one spectra by sample, spectra to choose
#' (either "first", "last" or its number)
#' @param library.metabolites path of the library of standard if not the default
#' one
#' @param threshold.noise threshold noise
#' @keywords NMR quantification metabolites
#' @export
#' @examples
#' result <- ASICS(path = "./spectres_exemple/AG_faq_Beck01",
#'  exclusion.areas = matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE),
#'  max.shift = 0.02, which.spectra = "last", library.metabolites = NULL)
ASICS <- function(path, exclusion.areas = matrix(c(4.5, 5.1), ncol = 2, nrow = 1),
                  max.shift = 0.02, which.spectra = "last",
                  library.metabolites = NULL, threshold.noise = 0.02)
{

  ##Seed and variables declaration
  set.seed(12345)


  #-----------------------------------------------------------------------------
  #Import complexe mixture and pure spectra library
  import <- Total_Bruker(path = path, exclusion.areas = exclusion.areas,
                         max.shift = max.shift, which.spectra = which.spectra,
                         library.metabolites = library.metabolites)

  mixture <- import$mixture
  pure_library <- import$pure_library

  #number of points on library grid corresponding to maximum shift
  nb_points_shift <- floor(max.shift / (pure_library$grid[2] -
                                          pure_library$grid[1]))



  #-----------------------------------------------------------------------------
  #Baseline improvement of complex mixture

  #Minimum and minimum index for nb_interval of same length
  nb_interval <- 100
  interval_length <- floor(length(pure_library$grid) / nb_interval)

  minima <- rollapply(mixture, FUN = min,
                      width = interval_length, by = interval_length)

  minima_idx <- rollapply(mixture, FUN = which.min,
                          width = interval_length, by = interval_length) +
    (1:nb_interval - 1) * interval_length

  #Smooth grid with linear interpolation between minimum
  minima_idx <- c(1, minima_idx, length(pure_library$grid))
  minima <- c(0, minima, 0)
  fitted_baseline <- approx(minima_idx, minima, method = "linear",
                            n = length(pure_library$grid))

  #Corrected mixture
  mixture <- mixture - fitted_baseline$y * (fitted_baseline$y < 0)



  #-----------------------------------------------------------------------------
  #Cleaning step: remove metabolites that cannot belong to the mixture

  signal_mixture <- 1 * (mixture > threshold.noise)
  signal_mixture_shift <- grossir_ensemble(signal_mixture, nb_points_shift)
  signal_library <- 1 * (pure_library$spectra > 0.1)

  #keep metabolites for which signal is included in mixture signal
  metab_to_keep <- which(apply(signal_library, 2,
                               function(x)
                                 sum(x - signal_mixture_shift == 1) == 0))

  pure_lib_clean <- subset_library(pure_library, metab_to_keep)



  #-----------------------------------------------------------------------------
  #Sort metabolites by regression residuals

  #create a matrix with all possible shifted mixture in column
  max_shift_left <- c(tail(head(mixture, -nb_points_shift), -nb_points_shift),
                      rep(0, 2 * nb_points_shift))
  mixture_shift <- t(laply (.data = as.list(1:(nb_points_shift * 2 - 1)),
                            .fun = function (x) c(tail(max_shift_left, x),
                                                  head(max_shift_left, -x))))
  mixture_shift <- cbind(max_shift_left, mixture_shift)

  #compute least square between all possible shifted mixture and all metabolites
  XtX <- diag(colSums(pure_lib_clean$spectra^2))
  least_square_shift <- solve(XtX)%*%t(pure_lib_clean$spectra)%*%mixture_shift
  max_least_square <- apply(least_square_shift, 1, which.max)

  #optimal shift
  residuals_opti <- unlist(lapply(1:length(pure_lib_clean$name),
                                  function(x)
                                    sum((lm(mixture_shift[, max_least_square[x]]~
                                              pure_lib_clean$spectra[,x] -
                                              1)$residuals) ^ 2)))

  # metabolites sorted by decreasing residual sum
  order <- sort(residuals_opti, decreasing = FALSE, index.return = TRUE)$ix
  pure_lib_sorted <- subset_library(pure_lib_clean, order)



  #-----------------------------------------------------------------------------
  #Deformations of pure spectra

  #Noises and weights
  s1 <- 0.172 #standard deviation of multiplicative noise
  s2 <- 0.15 #standard deviation of additive noise
  noises <- abs(mixture) * s1 ^ 2 + s2 ^ 2
  mixture_weights <- 1 / noises

  #Linear regression between mixture and each pure spectra
  least_square <- lm(mixture~pure_lib_sorted$spectra - 1,
                     weights = mixture_weights)

  #Shifted library
  pure_lib_shifted <- pure_lib_sorted

  #Deform each spectrum
  pure_lib_shifted$spectra <- t(laply(as.list(1:ncol(pure_lib_sorted$spectra)),
                                      deform_spectra, pure_lib_sorted, least_square,
                                      mixture_weights, nb_points_shift, max.shift))


  #-----------------------------------------------------------------------------
  #Localized deformations of pure spectra
  nb_iter_deform_loc <- 5
  nb_iter_deform <- 4
  range_a <- -9:9/10
  seuil1 <- 1

  nb_iter <- 0
  while(nb_iter < nb_iter_deform_loc){

    #####################!! script d'origine : pour moi il faut réutiliser la
    #fonction précédente mais il y quelques petites choses qui ne sont pas
    #faite pareilles
    MC <- lm(mixture~pure_lib_shifted$spectra - 1, weights = mixture_weights)
    r <- as.numeric(MC$residuals)


    for (i in 1:ncol(pure_lib_shifted$spectra))
    {

      c <- as.numeric(MC$coefficients[i])
      ##recupere les extremites des composantes connexes##
      sel <- which(pure_lib_shifted$spectra[, i] > seuil1)
      u <- which((sel[-1] - sel[1:(length(sel)-1)]-1) != 0)
      nb_int <- length(u) + 1
      v <- matrix(nrow = nb_int, ncol = 2)
      v[,2] <- c(sel[u], max(sel))
      v[,1] <- c(min(sel), sel[u+1])
      ####
      for(j in 1:nb_int)
      {
        x_res <- pure_lib_shifted$grid[v[j, 1]:v[j, 2]]
        r_res <- r[v[j, 1]:v[j, 2]]
        i0 <- which.max(abs(r_res)) ######## Pourquoi partir de là ??
        amin <- x_res[i0] - 0.002
        bmax <- x_res[i0] + 0.002
        r_opti <- r[(pure_lib_shifted$grid > amin) & (pure_lib_shifted$grid < bmax)]
        x_def <- pure_lib_shifted$grid[(pure_lib_shifted$grid > amin) & (pure_lib_shifted$grid < bmax)]


        adeformer <- pure_lib_shifted$spectra[(pure_lib_shifted$grid > amin) & (pure_lib_shifted$grid < bmax), i]
        pondere_res <- mixture_weights[(pure_lib_shifted$grid > amin) & (pure_lib_shifted$grid < bmax)]
        reste_res <- r_opti + c*pure_lib_shifted$spectra[(pure_lib_shifted$grid > amin) & (pure_lib_shifted$grid < bmax), i]
        for (k in 1:nb_iter_deform)
        {
          opti <- abs(sum(adeformer*reste_res*pondere_res)/sqrt(sum((adeformer^2)*pondere_res)))
          flag <- F
          for (a in range_a) ## même range_a que précédemment
          {
            candidat <- deforme(x_def, adeformer, a)

            if (abs(sum(candidat*reste_res*pondere_res)/sqrt(sum((candidat^2)*pondere_res))) > opti)
            {
              adeformer <- candidat
              opti <- abs(sum(candidat*reste_res*pondere_res)/sqrt(sum((candidat)^2*pondere_res)))
              flag <- T
            }
          }
          if (flag == T) pure_lib_shifted$spectra[(pure_lib_shifted$grid > amin) & (pure_lib_shifted$grid < bmax), i] <- adeformer
        }

      }

      pure_lib_shifted$spectra[, i] <- pure_lib_shifted$spectra[, i]/AUC(pure_lib_shifted$grid, pure_lib_shifted$spectra[, i])
    }
    ####################################
    nb_iter <- nb_iter + 1
  }





  #-----------------------------------------------------------------------------
  #Lasso optimisation
  ### Je n'ai pas changé les noms des variables ou essayer d'optimiser cette partie
  #du code

  #Construction de la matrice de variance du maximum de vraisemblance
  U <- 1 / sqrt(noises)
  A <- as.numeric(U) * pure_lib_shifted$spectra
  VMLE <- solve(t(A)%*%A) #La matrice du MLE est (X^T %*% D^{-1} %*% X)^{-1}

  N <- 1000 # Pour la minimisation (grossière des seuils), on prend peu d'observation de ZMLE
  C <- t(chol(VMLE))
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N, rnorm(nrow(C)*N))

  ##### parametre de regulariusation du lasso sous contraintes positives #######
  se <- sqrt(diag(VMLE))


  #Optimisation de la norme 1 des seuils

  Nbtirage <- 30
  p <- ncol(VMLE) #nb de metab
  W <- matrix(nrow = Nbtirage, ncol = p)
  u <- numeric(Nbtirage)

  delta0 <- rep(0.5, p) # 0.5 par metab
  a_min <- sum(seuil(delta0, ZMLE, se)) # somme des seuil par metab
  err <- 0.4

  for(i in 1:400) ## pk 400 ?
  {
    err <- 0.99 * err
    W <- matrix(delta0 + err * c(runif(p * Nbtirage, -1, 1)), nrow = Nbtirage,
                byrow = T)
    W[W > 1]  <-  1
    W[W < 0] <- 0.0001
    u <- apply(W, 1, aminimiser, ZMLE, se)
    if (min(u) < a_min)
    {
      delta0 <- W[which.min(u), ]
      a_min <- min(u)
    }
  }


  #Estimation du pseudo MLE
  B2 <- as.numeric(lm(mixture~pure_lib_shifted$spectra - 1,
                      weights = mixture_weights)$coefficients)

  N <- 10000
  #Calcul des seuils on prend beaucoup de realisations de ZMLE
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N,rnorm(nrow(C)*N))

  ####### Estimation lasso de l'active avec contraintes de positivite ##########
  identified_metab <- (B2 > Tuning(delta0, ZMLE) / delta0) & (B2 > 0)
  pure_lib_identified <- subset_library(pure_lib_shifted, which(identified_metab))

  B_final_tot <- as.numeric(lm(mixture~pure_lib_identified$spectra - 1, weights = mixture_weights)$coefficients)

  #Test de la positivité des coefficients
  B_final <- B_final_tot[B_final_tot > 0]

  identified_metab[identified_metab][B_final_tot < 0] <- FALSE

  pure_lib_final <- subset_library(pure_lib_shifted, which(identified_metab))







  #-----------------------------------------------------------------------------
  #Results

  #Reconstituted mixture with estimated coefficents
  est_mixture <- pure_lib_final$spectra%*%B_final

  #Compute relative concentration of identified metabolites
  relative_concentration <- B_final / pure_lib_final$nb_protons
  #Sort library according to relative concentration
  sorted_idx <- sort(relative_concentration, decreasing = TRUE,
                     index.return = TRUE)$ix
  pure_lib_final_sorted <- subset_library(pure_lib_final, sorted_idx)

  present_metab <- data.frame(Name = pure_lib_final_sorted$name,
                              Relative_Concentration =
                                relative_concentration[sorted_idx])


  #Compute identification threshold of non identified metabolites
  identification_threshold <- round((seuil(delta0, ZMLE, se)[!identified_metab] *
                                       pure_lib_final_sorted$nb_protons[1]) /
                                      (pure_lib_shifted$nb_protons[!identified_metab]),
                                    digits = 4)

  non_identified_metab <- data.frame(Name = pure_lib_shifted$name[!identified_metab],
                                     Threshold = identification_threshold)

  #List to return
  L <- list("Present" = present_metab, "Non_identified" = non_identified_metab,
            "Mixture" = mixture, "Estimated_mixture" = est_mixture,
            "Grid" = pure_lib_final$grid)
  return(L)
}

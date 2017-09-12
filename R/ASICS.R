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
#' @param threshold.noise threshold noised
#' @export
#' @examples
#' result <- ASICS(path = "./spectres_exemple/AG_faq_Beck01",
#'  exclusion.areas = matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE),
#'  max.shift = 0.02, which.spectra = "last", library.metabolites = NULL)
ASICS <- function(path, exclusion.areas = matrix(c(4.5, 5.1), ncol = 2, nrow = 1),
                  max.shift = 0.02, which.spectra = "last",
                  library.metabolites = NULL, threshold.noise = 0.02){

  T1 <- Sys.time()

  ##Seed and variables declaration
  set.seed(12345)


  #-----------------------------------------------------------------------------
  #### Import complexe mixture and pure spectra library ####
  import <- load_data(path = path, exclusion.areas = exclusion.areas,
                      max.shift = max.shift, which.spectra = which.spectra,
                      library.metabolites = library.metabolites)

  mixture <- import$mixture
  pure_library <- import$pure_library

  #number of points on library grid corresponding to maximum shift
  nb_points_shift <- floor(max.shift / (pure_library$grid[2] -
                                          pure_library$grid[1]))

  T2 <- Sys.time()


  #-----------------------------------------------------------------------------
  #### Baseline improvement of complex mixture : remove negative values --> dans script d'import ####
  mixture[mixture < 0] <- 0
  mixture <- mixture / AUC(pure_library$grid, mixture)


  #-----------------------------------------------------------------------------
  #### Cleaning step: remove metabolites that cannot belong to the mixture ####
  signal_mixture <- 1 * (mixture > threshold.noise)
  signal_mixture_shift <- signal_with_shift(signal_mixture, nb_points_shift)

  #normalize each spectra with the maximum of mixture
  norm_library <- t(aaply(pure_library$spectra, 2, function(x) x * max(mixture) /
                            max(x)))
  signal_library <- 1 * (norm_library > threshold.noise)

  #keep metabolites for which signal is included in mixture signal
  metab_to_keep <- which(apply(signal_library, 2,
                               function(x)
                                 sum(x - signal_mixture_shift == 1) == 0))

  pure_lib_clean <- subset_library(pure_library, metab_to_keep)


  T3 <- Sys.time()

  #-----------------------------------------------------------------------------
  #### Find the best translation between each pure spectra and mixture ####
  #and sort metabolites by regression residuals
#### !! attention au shif maximal des 2 déformations
  #create a matrix with all possible shifted mixture in column
  mixture_all_shift <- c(rep(0, nb_points_shift),
                         mixture,
                         rep(0, nb_points_shift))

  mixture_shift <- t(laply(.data = as.list(1:(nb_points_shift * 2 + 1)),
                           .fun = function (x)
                             mixture_all_shift[x:(x + length(mixture) - 1)]))

  #create a matrix of all pure spectra with shift
  spectra_all_shift <- rbind(matrix(0, nrow = nb_points_shift,
                                    ncol = ncol(pure_lib_clean$spectra)),
                             pure_lib_clean$spectra,
                             matrix(0, nrow = nb_points_shift,
                                    ncol = ncol(pure_lib_clean$spectra)))

  #compute least square between all possible shifted mixture and all metabolites
  XtX <- diag(colSums(pure_lib_clean$spectra^2))
  least_square_shift <- solve(XtX)%*%t(pure_lib_clean$spectra)%*%mixture_shift
  max_least_square <- apply(least_square_shift, 1, which.max)

  #shift library according to least square
  pure_lib_shifted <- pure_lib_clean
  pure_lib_shifted$spectra <- t(laply(as.list(1:ncol(pure_lib_clean$spectra)),
                                      function(i)
                                        spectra_all_shift[
                                          (ncol(mixture_shift) - max_least_square[i] +
                                             1):
                                            (ncol(mixture_shift) - max_least_square[i] +
                                               nrow(pure_lib_shifted$spectra)), i]))


  #optimal residuals
  residuals_opti <- unlist(lapply(1:length(pure_lib_shifted$name),
                                  function(x)
                                    sum((lm(mixture_shift[, max_least_square[x]]~
                                              pure_lib_shifted$spectra[,x] -
                                              1)$residuals) ^ 2)))


  #metabolites sorted by decreasing residual sum
  order <- sort(residuals_opti, decreasing = FALSE, index.return = TRUE)$ix
  pure_lib_sorted <- subset_library(pure_lib_shifted, order)


  #-----------------------------------------------------------------------------
  #### Localized deformations of pure spectra ####

  #noises and weights
  s1 <- 0.172 #standard deviation of multiplicative noise
  s2 <- 0.15 #standard deviation of additive noise
  noises <- abs(mixture) * s1 ^ 2 + s2 ^ 2
  mixture_weights <- 1 / noises

### !!! tester en mettant à jour les résidus

  #Shifted library
  pure_lib_deformed <- pure_lib_sorted

  #Algorithm parameters
  peak_threshold <- 1
  nb_iter_by_peak <- 4
  range_a <- -9:9/10
  nb_iter_by_library <- 5

  #Linear regression between mixture and each pure spectra
  updated_least_square <- lm(mixture~pure_lib_deformed$spectra - 1,
                             weights = mixture_weights)
  LS_residuals <- as.numeric(residuals(updated_least_square))

  nb_iter_lib <- 0
  while(nb_iter_lib < nb_iter_by_library){

    #Deforme each peak of a pure spectrum to align it on the complex mixture
    for (i in 1:ncol(pure_lib_deformed$spectra)) {

      #Linear regression coefficient of spectrum i
      LS_coeff <- max(0, as.numeric(updated_least_square$coefficients[i]))

      #Expanded connected components
      signal_peak_lib <- which(pure_lib_deformed$spectra[, i] > peak_threshold)

      min_extremities <- signal_peak_lib[!((signal_peak_lib - 1) %in%
                                             signal_peak_lib)] - floor(nb_points_shift / 10)
      max_extremities <- signal_peak_lib[!((signal_peak_lib + 1) %in%
                                             signal_peak_lib)] + floor(nb_points_shift / 10)
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
        to_deform <- pure_lib_deformed$spectra[peak_area, i]
        #grid to deform
        grid_to_deform <- pure_lib_deformed$grid[peak_area]

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
               max(abs(deformed_grid - pure_lib_sorted$grid[peak_area])) < max.shift){
              to_deform <- deformed_spectrum_small
              opti_criterion <- new_opti
            }
          }
          iter <- iter + 1
        }
        pure_lib_deformed$spectra[peak_area, i] <- to_deform
      }

      pure_lib_deformed$spectra[, i] <- pure_lib_deformed$spectra[, i] /
        AUC(pure_lib_deformed$grid, pure_lib_deformed$spectra[, i])
    }

    nb_iter_lib <- nb_iter_lib + 1
  }


  T4 <- Sys.time()








  #-----------------------------------------------------------------------------
  #### Lasso optimisation ????????????? ####
  ### Je n'ai pas changé les noms des variables ou essayer d'optimiser cette partie
  #du code

  #Construction de la matrice de variance du maximum de vraisemblance
  U <- 1 / sqrt(noises)
  A <- as.numeric(U) * pure_lib_shifted$spectra
  VMLE <- solve(t(A)%*%A) #La matrice du MLE est (X^T %*% D^{-1} %*% X)^{-1}

  N <- 1000 # Pour la minimisation (grossière des seuils), on prend peu d'observation de ZMLE
  C <- t(chol(VMLE))
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N, rnorm(nrow(C)*N))

  ##### parametre de regulariusation du lasso sous contraintes positives
  se <- sqrt(diag(VMLE))


  #Optimisation de la norme 1 des seuils

  Nbtirage <- 30
  p <- ncol(VMLE) #nb de metab
  W <- matrix(nrow = Nbtirage, ncol = p)
  u <- numeric(Nbtirage)

  delta0 <- rep(0.5, p) # 0.5 par metab
  a_min <- sum(seuil(delta0, ZMLE, se)) # somme des seuil par metab
  err <- 0.4

  for(i in 1:400)
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

  # Estimation lasso de l'active avec contraintes de positivite
  identified_metab <- (B2 > Tuning(delta0, ZMLE) / delta0) & (B2 > 0)
  pure_lib_identified <- subset_library(pure_lib_shifted, which(identified_metab))

  B_final_tot <- as.numeric(lm(mixture~pure_lib_identified$spectra - 1, weights = mixture_weights)$coefficients)

  #Test de la positivité des coefficients
  B_final <- B_final_tot[B_final_tot > 0]

  identified_metab[identified_metab][B_final_tot < 0] <- FALSE

  pure_lib_final <- subset_library(pure_lib_shifted, which(identified_metab))







  #-----------------------------------------------------------------------------
  #### Results ####

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

  T5 <- Sys.time()
  return(L)
}

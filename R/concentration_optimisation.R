#### Threshold and concentration optimisation for each metabolites ####
#' @importFrom stats rnorm runif
#' @importFrom methods is
#' @keywords internal
concentration_opti <- function(mixture, pure_lib_deformed, nb.iter.signif){
  #noises and weights
  s1 <- 0.172 #standard deviation of multiplicative noise
  s2 <- 0.15 #standard deviation of additive noise
  noises <- abs(mixture) * s1 ^ 2 + s2 ^ 2
  mixture_weights <- 1 / noises

  #Variance matrix of maximum likelihood
  A <- 1 / sqrt(noises) * pure_lib_deformed$spectra
  VMLE <- solve(t(A)%*%A)

  #First threshold minimisation
  N <- 1000
  C <- t(chol(VMLE))
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N, rnorm(nrow(C)*N))

  #Regularization paramater of lasso under positive constraints
  se <- sqrt(diag(VMLE))

  #L1 norm of threshold optimisation
  nb_draw <- 30
  p <- ncol(VMLE)
  W <- matrix(nrow = nb_draw, ncol = p)
  u <- numeric(nb_draw)

  delta0 <- rep(0.5, p) # 0.5 by metabolite
  a_min <- sum(compute_threshold(delta0, ZMLE, se)) #threshold sum by metabolite
  err <- 0.4

  for(i in 1:nb.iter.signif) {
    err <- 0.99 * err
    W <- matrix(delta0 + err * c(runif(p * nb_draw, -1, 1)), nrow = nb_draw,
                byrow = TRUE)
    W[W > 1]  <-  1
    W[W < 0] <- 0.0001
    u <- apply(W, 1, to_minimize, ZMLE, se)
    if (min(u) < a_min)
    {
      delta0 <- W[which.min(u), ]
      a_min <- min(u)
    }
  }

  #Pseudo MLE estimation
  B2 <- try(lm_constrained(mixture, pure_lib_deformed$spectra,
                           mixture_weights)$coefficients, silent = TRUE)
  if(is(B2, "try-error")){
    B2 <- lm_constrained(mixture, pure_lib_deformed$spectra,
                             mixture_weights, 10e-3)$coefficients
  }

  #Compute all thresholds
  N <- 10000
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N, rnorm(nrow(C)*N))

  #Concentration lasso estimation with positive constraints
  identified_metab <- (B2 > tuning(delta0, ZMLE) / delta0) & (B2 > 0)
  pure_lib_identified <- subset_library(pure_lib_deformed,
                                        which(identified_metab))

  B_final_tot <- try(lm_constrained(mixture, pure_lib_identified$spectra,
                                    mixture_weights)$coefficients,
                     silent = TRUE)
  if(is(B_final_tot,"try-error")){
    B_final_tot <- lm_constrained(mixture, pure_lib_identified$spectra,
                                 mixture_weights, 10e-3)$coefficients
  }

  #Test of coefficients positivity
  B_final <- B_final_tot[B_final_tot > 0]
  identified_metab[identified_metab][B_final_tot < 0] <- FALSE
  pure_lib_final <- subset_library(pure_lib_deformed, which(identified_metab))

  return(list(pure_lib_final = pure_lib_final, B_final = B_final,
              threshold = compute_threshold(delta0, ZMLE, se),
              identified_metab = identified_metab))
}



##### Threshold optimisation functions
to_minimize <- function(x, ZMLE, se) {
  sum(compute_threshold(x, ZMLE, se))
}

#' @importFrom stats qnorm
#' @keywords internal
compute_threshold <- function(x, ZMLE, se) {
  return(tuning(x, ZMLE)/x + se*qnorm(0.95))
}

#' @importFrom stats quantile
#' @keywords internal
tuning <- function(delta, ZMLE) {
  observation <- apply(as.double(delta) * (ZMLE), 2, max)
  return(quantile(observation, 0.95))
}

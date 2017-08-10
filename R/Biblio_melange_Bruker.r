f_o_phi <- function(x, f, phix)
{
  #il y a des pbs avec les ex aequo dans u et v
  # Attention, il faut que inf phi(x) >= inf x ET sup phi(x)<= sup x
  #phix contient le vecteur phi(x), f contient f(x)
  #renvoie f(phi(x))
  if (phix[1] < x[1]) phix[1] <- x[1]
  n <- length(x)
  nphi <- length(phix)
  if (phix[nphi] > x[n]) phix[nphi] <- x[n]
  orde <- rep(c(1, 0), c(n, nphi))
  x[1] <- x[1] - 1e-06
  x[n] <- x[n] + 1e-06
  zx <- c(x, phix)
  o <- order(zx)
  u <- cumsum(orde[o])
  v <- u[orde[o] == 0]
  fophi <- f[v] + (phix - x[v])*(f[(v + 1)] - f[v])/(x[(v + 1)] - x[v])
  if (is.na(fophi[nphi])) fophi[nphi] <- f[n]
  if (is.na(fophi[1])) fophi[1] <- f[1]
  return(fophi)
}



AUC <- function(x, y) # calcul l'AUC de y=f(x) avec la methode des trapezes
{
  a <- diff(x)
  auc <- 0.5*sum((y[2:length(y)] + y[1:length(y) - 1])*a)
  return(auc)
}


# crunch <- function(x, nb_interval=6000)
# {
#   re1 <- sum(x[,1] < 4.5)
#   re2 <- sum(x[,1] > 5)
#   a <- floor((re1 + re2)/(nb_interval))
#   #attention au trou entre 4.5 et 5#
#   bb1 <- rep(1:round(re1/a + 1), length.out = re1, each = a)
#   bb2 <- rep((max(bb1) + 1):(max(bb1) + 1 + round(re2/a + 1)), length.out = re2, each = a)
#   bb <- c(bb1, bb2)
#
#   yy <- matrix(nrow = max(bb), ncol = 2)
#   bb <- as.factor(bb)
#   yy[,1] <- tapply(x[, 1], bb, mean)
#   yy[,2] <- tapply(x[, 2], bb, sum)
#   return(as.matrix(yy))
# }

subset_library <- function(pure_library, idx){
  sub_library <- list()
  sub_library$name <- pure_library$name[idx]
  sub_library$grid <- pure_library$grid
  sub_library$spectra <- pure_library$spectra[, idx]
  sub_library$nb_protons <- pure_library$nb_protons[idx]
  return(sub_library)
}








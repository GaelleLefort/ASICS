#' Automatic Statistical Identification in Complex Spectra
#'
#' Description
#' @param name folder path of the Bruker files
#' @param ZoneAEnlever exclusion areas to remove before the quantification
#' @param DecalageMax maximum chemical shift allowed
#' @param DecalageLib chemical shift between library of pure spectra and complex mixture
#' @param sample if more than one spectra by sample, spectra to choose (either "first",
#' "last" or its number
#' @param seed seed
#' @param libraryMetab path of the library of standard if not the default one
#' @param seuilBruit threshold noise
#' @keywords NMR quantification metabolites
#' @export
#' @examples
#' result <- ASICS(name = "./spectres_exemple/AG_faq_Beck01",
#'  ZoneAEnlever = matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE),
#'  DecalageMax = 0.02, DecalageLib = 0, sample = "last", seed = 12345,
#'  libraryMetab = NULL)
ASICS <- function(name, ZoneAEnlever = matrix(c(4.5, 5.1), ncol = 2, nrow = 1),
                  DecalageMax = 0.02, DecalageLib = 0, sample = "last", seed = 12345,
                  libraryMetab = NULL, seuilBruit = 0.02)
{

  s1 <- 0.172 #ecart type du bruit multiplicatif
  s2 <- 0.15 #ecart type du bruit additif1
  set.seed(seed)

  ##On va chercher le melange et la bibliotheque de metabolites
  a_analyser <- Total_Bruker(name = name, ZoneAEnlever = ZoneAEnlever,
                             DecalageMax = DecalageMax, DecalageLib = DecalageLib,
                             sample = sample, libraryMetab = libraryMetab)
  melange <- a_analyser$Melange
  nameMeta <- a_analyser$Library$Name
  x <- a_analyser$Library$Grid
  biblio <- a_analyser$Library$Metab
  dep <- DecalageMax
  h <- floor(dep/(x[2] - x[1]))
  proton <- a_analyser$Library$Protons


  #####################Amelioration de la ligne de base ##############################
  N <- floor(length(x)/100)
  m <- numeric(100)
  r <- numeric(100)
  for(i in 1:100)
  {
    m[i] <- min(melange[(1+(i-1)*N):(i*N)])
    r[i] <- which.min(melange[(1+(i-1)*N):(i*N)])+(i-1)*N
  }

  cuts <- c(1, r, length(x))
  vals <- c(0, m, 0)

  f <- approx(cuts, vals, method = "linear", n = length(x))

  melange <- melange - f$y*(f$y<0)

  #####################################################################################
  seuil_melange <- seuilBruit
  Smelange <- 1*(melange > seuil_melange)
  Saugmente <- grossir_ensemble(Smelange, h)
  Test <- numeric(ncol(biblio))
  for(i in 1:ncol(biblio))
  {
    Smetab <- 1*(biblio[,i] > 0.1) #Support du metabolite i (on doit etre sur de ne selectionner que le support
    Test[i] <- 1*(sum(Smetab*(1 - Saugmente)) == 0) #Test si le support du metabolite est inclus dans le support du melange augmente
  }
  v0 <- (1:ncol(biblio))[Test == 1] #Metabolites non elimines apres l'etape de nettoyage.
  nameMeta1 <- nameMeta[v0]


  ##################### Etape pour ordonner les metabolites #####################

  u <- melange[(h+1):(length(x)-h)]
  M <- matrix(nrow = length(x),ncol = (2*h))
  #matrice ou chaque colonne est le spectre du melange deforme
  for(i in (0:(2*h-1)))
  {
    M[, (i+1)] <- c(rep(0, i),u,rep(0, (2*h-i)))
  }
  D <- diag(colSums(biblio[, v0]^2))
  A <- solve(D)%*%t(biblio[, v0])%*%M #matrice des coefficients des moindres carres ; chaque ligne correspond au coefficient des
  #moindres carres entre un metabolite et les melanges deformes

  R <- apply(A, 1, which.max)
  Resultat <- numeric(length(v0))

  for(i in (1:length(v0)))
  {
    Resultat[i] <- sum((lm(M[, R[i]]~biblio[, v0[i]]-1)$residuals)^2)
  }

  order <- sort(Resultat, decreasing = FALSE, index.return = TRUE)$ix
  v1 <- v0[order]
  proton1 <- proton[v1]

  ###############################################################################

  biblio_a_conserver <- biblio[,v1]
  nameMeta_final <- nameMeta[v1]

  phi <- function(x, a)
  {
    u <- min(x)
    v <- max(x) - min(x)
    z <- (x - u)/v
    tt <- z + a*z*(1 - z)
    return(u + tt*v)
  }
  deforme <- function(x, y, a)
  {
    phix <- phi(x, a)
    return(f_o_phi(x, y, phix)) # phix deformation de l'axe des absisses et f_o_phi est la valeur de la spectre deforme sur les points
    #de discretisation de l'axe ds absisses.
  }


  ############################ Reprend l'etape de deformation #########################

  var <- abs(melange)*s1^2 + s2^2
  pondere <- 1/var
  MC <- lm(melange~biblio_a_conserver - 1, weights = pondere)
  biblio_deformee <- biblio_a_conserver
  rangea <- -9:9/10
  nb_iter_deform <- 4
  seuil1 <- 1

  for (i in 1:ncol(biblio_deformee))
  {
    ##recupere les extremites des composantes connexes##
    control <- x
    sel <- which(biblio_deformee[, i] > seuil1)
    u <- which((sel[-1] - sel[1:(length(sel) - 1)] - 1) != 0)
    nb_int <- length(u) + 1
    v <- matrix(nrow = nb_int, ncol = 2)
    v[, 2]=c(sel[u], max(sel))
    v[, 1]=c(min(sel), sel[u + 1])
    ####
    ##agrandi les composantes connexes##
    v[, 1]=v[, 1] - h
    v[, 2]=v[, 2] + h
    ####

    ###si elles se chevauchent on les fusionne
    # iii<-1
    #while(iii<nrow(matrix(v,ncol=2)))
    #{   if (v[iii,2]>v[iii+1,1])
    #    {v[iii,2]<-v[iii+1,2]
    #      v<-v[-(iii+1),]
    #     }
    #      else
    #      {iii<-iii+1}
    #}
    #v<-matrix(v,ncol=2)
    #nb_int<-nrow(v)

    #reste <- as.numeric(MC$residuals) + as.numeric(MC$coefficients[i])*biblio_deformee[, i]
    reste <- as.numeric(MC$residuals) + max(0, as.numeric(MC$coefficients[i]))*biblio_deformee[, i]

    ## deforme sur chaque composante connexe##
    for (j in 1:nb_int)
    {
      #_res signifie que l'on se restreint à un intervalle
      reste_res <- reste[v[j, 1]:v[j, 2]]
      biblio_res <- biblio_deformee[v[j, 1]:v[j, 2],] #restriction des spectres de metabolite sur la composante connexe
      adeformer <- biblio_res[, i] #partie du spectre du metabolite a deformer
      x_res <- x[v[j, 1]:v[j, 2]]
      pondere_res <- pondere[v[j, 1]:v[j, 2]]
      tt <- x_res
      for (k in 1:nb_iter_deform)
      {
        opti <- abs(sum(adeformer*reste_res*pondere_res)/sqrt(sum((adeformer^2)*pondere_res)))
        flag <- F
        for (a in rangea)
        {
          control[v[j, 1]:v[j, 2]] <- phi(tt, a)
          candidat <- deforme(x_res, adeformer, a)

          if ((abs(sum(candidat*reste_res*pondere_res)/sqrt(sum((candidat^2)*pondere_res))) > opti) & (max(abs(control - x)[v[j, 1]:v[j, 2]]) < dep))
          {
            adeformer <- candidat
            opti <- abs(sum(candidat*reste_res*pondere_res)/sqrt(sum((candidat)^2*pondere_res)))
            flag <- T
            tt <- phi(tt, a)
          }
        }
        if (flag == T) biblio_deformee[v[j, 1]:v[j, 2], i] <- adeformer
      }

    }

    biblio_deformee[, i] <- biblio_deformee[, i]/AUC(x, biblio_deformee[, i])
  }


  ###################### Deformation localisee #######################
  biblio_deformeee <- biblio_deformee
  for(l in 1:5)
  {

    MC <- lm(melange~biblio_deformeee - 1, weights = pondere)
    r <- as.numeric(MC$residuals)


    for (i in 1:ncol(biblio_deformeee))
    {

      c <- as.numeric(MC$coefficients[i])
      ##recupere les extremites des composantes connexes##
      sel <- which(biblio_deformeee[, i] > seuil1)
      u <- which((sel[-1] - sel[1:(length(sel)-1)]-1) != 0)
      nb_int <- length(u) + 1
      v <- matrix(nrow = nb_int, ncol = 2)
      v[,2] <- c(sel[u], max(sel))
      v[,1] <- c(min(sel), sel[u+1])
      ####
      for(j in 1:nb_int)
      {
        x_res <- x[v[j, 1]:v[j, 2]]
        r_res <- r[v[j, 1]:v[j, 2]]
        i0 <- which.max(abs(r_res))
        amin <- x_res[i0] - 0.002
        bmax <- x_res[i0] + 0.002
        r_opti <- r[(x > amin) & (x < bmax)]
        x_def <- x[(x > amin) & (x < bmax)]
        adeformer <- biblio_deformeee[(x > amin) & (x < bmax), i]
        pondere_res <- pondere[(x > amin) & (x < bmax)]
        reste_res <- r_opti + c*biblio_deformeee[(x > amin) & (x < bmax), i]
        for (k in 1:nb_iter_deform)
        {
          opti <- abs(sum(adeformer*reste_res*pondere_res)/sqrt(sum((adeformer^2)*pondere_res)))
          flag <- F
          for (a in rangea)
          {
            candidat <- deforme(x_def, adeformer, a)

            if (abs(sum(candidat*reste_res*pondere_res)/sqrt(sum((candidat^2)*pondere_res))) > opti)
            {
              adeformer <- candidat
              opti <- abs(sum(candidat*reste_res*pondere_res)/sqrt(sum((candidat)^2*pondere_res)))
              flag <- T
            }
          }
          if (flag == T) biblio_deformeee[(x > amin) & (x < bmax), i] <- adeformer
        }

      }

      biblio_deformeee[, i] <- biblio_deformeee[, i]/AUC(x, biblio_deformeee[, i])
    }
  }

  ####################### Optimisation lasso ##########################

  #Construction de la matrice de variance du maximum de vraisemblance
  U <- 1/sqrt(var)
  A <- as.numeric(U)*biblio_deformeee
  VMLE <- solve(t(A)%*%A) #La matrice du MLE est (X^T %*% D^{-1} %*% X)^{-1}

  N <- 1000 # Pour la minimisation (grossière des seuils), on prend peu d'observation de ZMLE
  C <- t(chol(VMLE))
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N, rnorm(nrow(C)*N))

  ################## parametre de regulariusation du lasso sous contraintes positives #############

  Tuning <- function(delta)
  {
    Observation <- apply(as.double(delta)*(ZMLE), 2, max)
    return(quantile(Observation, 0.95))
  }

  #########################################################

  se <- sqrt(diag(VMLE))

  seuil <- function(x)
  {
    return(Tuning(x)/x + se*qnorm(0.95))
  }


  aminimiser <- function(x)
  {
    sum(seuil(x))
  }

  #Optimisation de la norme 1 des seuils

  Nbtirage <- 30
  p <- ncol(VMLE)
  W <- matrix(nrow = Nbtirage, ncol = p)
  u <- numeric(Nbtirage)

  delta0 <- rep(0.5, p)
  a_min <- sum(seuil(delta0))
  err <- 0.4

  for(i in 1:400)
  {
    err <- 0.99*err
    W <- matrix(delta0 + err*c(runif(p*Nbtirage, -1, 1)), nrow = Nbtirage, byrow = T)
    W[W >1 ]  <-  1
    W[W < 0] <- 0.0001
    u <- apply(W, 1, aminimiser)
    if (min(u) < a_min)
    {
      delta0 <- W[which.min(u), ]
      a_min <- min(u)
    }
  }


  #Estimation du pseudo MLE
  B2 <- as.numeric(lm(melange~biblio_deformeee - 1, weights = pondere)$coefficients) #estimation du pseudo MLE

  N <- 10000
  ZMLE <- C%*%matrix(nrow = nrow(C), ncol = N,rnorm(nrow(C)*N)) #Calcul des seuils on prend beaucoup de realisations de ZMLE

  ################## Estimation lasso de l'active avec contraintes de positivite #############
  v <- 1*(B2 > Tuning(delta0)/delta0)*1*(B2 > 0)
  B_final_tot <- as.numeric(lm(melange~biblio_deformeee[, v==1] - 1, weights = pondere)$coefficients)
  Positif <- B_final_tot > 0
  B_final <- B_final_tot[Positif]
  if (any(Positif == FALSE))
  {
    Compo_Negative <- v1[v == 1][Positif == FALSE]
    for (j in 1:length(Compo_Negative))
    {
      v[v1 == Compo_Negative[j]] <- 0
    }   ##on met les negatifs dans non identifies
  }

  Metamelange <- nameMeta_final[v == 1]

  ############ Reconstitution du melange avec les metabolites observes ##########
  Melange_rec <- biblio_deformeee[, v == 1]%*%B_final #spectre du melange reconstitue

  ############## Concentrations relatives ###########################
  v2 <- v1[v == 1] #Metabolite conserver à l'issue de l'etape d'elimination et du lasso
  proton_final <- proton[v2]
  concentration_final <- B_final/proton_final #concentration est proportionnelle à l'AUC et est inversement proportionnelle au nombre de proton
  ord <- sort(concentration_final, decreasing = TRUE, index.return = TRUE)$ix
  Metaordonnee <- Metamelange[ord]

  Metabolite_ref <- Metaordonnee[1] #Metabolite de reference
  B_ord <- B_final[ord]
  proton_ord <- proton_final[ord]
  concentration_ord <- concentration_final[ord]

  #concentration relatives à celle de Metabolite_ref
  #C_final=(B_ord[-1]*proton_ord[1])/(B_ord[1]*proton_ord[-1])

  ##!!!!!!
  #C_final=(concentration_ord[-1])/(concentration_ord[1])
  C_final <- concentration_ord

  # Seuil d'identification relatif
  Meta_non_identifie <- nameMeta_final[v == 0]
  seuil_identification_relatif <- round((seuil(delta0)[v == 0]*proton_ord[1])/(proton[v1[v == 0]]), digits = 4)

  ##A renvoyer
  ##!!!
  Present <- data.frame(Name = Metaordonnee, Relative_Concentration = C_final)
  Non_identified <- data.frame(Name = Meta_non_identifie, Threshold = seuil_identification_relatif)
  L <- list("Present" = Present, "Non_identified" = Non_identified, "Mixture" = melange,
            "Estimated_mixture" = Melange_rec, "Grid" = a_analyser$Library$Grid)
  return(L)
}

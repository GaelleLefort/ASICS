######Exemple de programme ##############################
#Pour le parallelisme
Total_Bruker <- function(name, libraryMetab = NULL, extensionMelange = "Bruker",
                         ZoneAEnlever = matrix(c(4.5,5.1), ncol = 2, nrow = 1),
                         DecalageMax = 0.02, DecalageLib = 0, sample = "last")
{
  #name est le nom du dossier dans lequel on a le melange
  #DecalageMax est le decalage qu'on permet sur le spectre en ppm
  # [borneinfZoneEnlever;bornesupZoneEnlever] zone qu'on enlève, souvent zone où il y a l'eau
  # en général [4.4;5.1] si on est dans l'urine ça peut etre [4.4;6.3] ou autre
  #extensionMelange = "txt" ou "Bruker"
  #libraryMetab est le chemin pour trouver la bibliotheque a partir de l'endroit ou les prog sont
  #library(linprog)
  #  library(doParallel)
  #  nbcoeurs<-trunc(detectCores()*0.75)
  #  cl<-makeCluster(nbcoeurs)
  #  registerDoParallel(cl)
  #  here<-getwd()
  #
  # #on dit ou est la bibliotheque
  # if (is.null(libraryMetab)) {libraryMetab<-paste(here,"/LibraryMetab",sep="")}
  # else {libraryMetab<-paste(here,"/",libraryMetab,sep="")}
  #
  #  #on va chercher les meta
  #   source("Biblio_Totale_Meta_Brucker.r")
  #   libraryMeta<-recup_biblio_meta(libraryMetab,borneinfZoneEnlever,bornesupZoneEnlever)
  #   data<-libraryMeta$data
  #   nameMeta<-libraryMeta$name
  #   pas<-data[[1]][2,1]-data[[1]][1,1]
  #

  if(is.null(libraryMetab)) {
    data(Library)
  } else {
    load(libraryMetab)
  }

  #Decaler la librairie
  Library$Grid <- Library$Grid + DecalageLib

  ##on va chercher le melange
  a_analyser <- recup_melange_bru(name, Library$Grid, sample = sample)

  ##on enleve la zone a enlever
  Zone <- vector(mode = "logical", length = length(Library$Grid))
  for (i in 1:nrow(ZoneAEnlever))
  {
    Zone <- Zone | ((Library$Grid > ZoneAEnlever[i, 1]) & (Library$Grid < ZoneAEnlever[i, 2]))
  }
  a_analyser[Zone, 2] <- 0
  Library$Metab[Zone, ] <- 0
  ###on teste si un metabolite n'a plus que des valeurs = 0
  testZero <- colSums(Library$Metab) > 0
  Library$Metab <- Library$Metab[, testZero]
  Library$Name <- Library$Name[testZero]
  Library$Protons <- Library$Protons[testZero]

  ##on renormalise les metabolites et le melange
  for (i in 1:ncol(Library$Metab))
  {
    Library$Metab[, i] <- Library$Metab[, i]/AUC(Library$Grid, Library$Metab[, i])
  }
  a_analyser[, 2] <- a_analyser[, 2]/AUC(Library$Grid, a_analyser[, 2])

  return(list("Melange" = a_analyser[, 2], "Library" = Library))
  #return(list("Metabolites"=data,"Noms_Metab"=nameMeta,"Melange"=a_analyser))

}


recup_melange_bru <- function(nameMel, grilleMeta, borneinfZoneEnlever = 0,
                              bornesupZoneEnlever = 0.02, sample = "last")
{
  # melange a analyser
  if(sample == "first"){
    nameSample <- file.path(nameMel, "1")
  } else if(sample == "last"){
    nameSample <- file.path(nameMel, dir(nameMel)[length(dir(nameMel))])
  } else {
    nameSample <- file.path(nameMel, sample)
  }


  a_analyserT <- NmRBrucker_read(nameSample, a_analyserT)
  a_analyserT[, 2] <- B_Corrector(a_analyserT[, 2])[[2]]

  a_analyserT[, 2] <- a_analyserT[, 2]/(AUC(a_analyserT[, 1], a_analyserT[, 2]))
  #on se remet sur la meme grille pour tous
  dataM <- matrix(0, nrow = length(grilleMeta), ncol = 2)
  dataM[, 2] <- f_o_phi(a_analyserT[, 1], a_analyserT[, 2], grilleMeta)
  dataM[, 1] <- grilleMeta
  return(dataM)
}

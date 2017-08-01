#' Extract concentrations
#'
#' Description
#' @param resASICS result of ASICS or ASICS_multiFiles functions
#' @export
#' @examples
#' quantification <- extractConcentrations(result_all)
extractConcentrations <- function(resASICS){
  relConc_list <- list()

  #Nom des échantillons dans chaque data.frame
  pos <- 1
  for(i in 1:length(resASICS)){
    if(class(resASICS[[i]]) != "try-error"){
      relConc_list[[pos]] <- resASICS[[i]]$Present
      colnames(relConc_list[[pos]])[2] <- names(resASICS)[i]
      pos <- pos + 1
    }

  }

  #Obtention d'un tableau à partir d'une liste
  relConc <- join_all(relConc_list, by = "Name", type = "full")
  colnames(relConc)[1] <- "metabolites"

  return(relConc)
}

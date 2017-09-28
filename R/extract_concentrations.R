#' Extract concentrations
#'
#' Combine results of multiple ASICS function to obtain quantified relative
#' concentration of each spectrum in one dataset
#'
#' @param res_ASICS result of ASICS_multiFiles function
#' @return A data frame containing relative concentrations of identified
#' metabolites for each spectrum
#' @export
#' @examples
#' \dontrun{
#' #Compute quantification on all spectra
#' res_multi <- ASICS_multiFiles(name.dir = system.file("extdata",
#'                                                      "example_spectra",
#'                                                      package = "ASICS"),
#'                            exclusion.areas = matrix(c(4.5,5.1,5.5,6.5),
#'                                                     ncol = 2, byrow = TRUE))
#'
#' #Extract relative concentrations
#' quantification <- extract_concentrations(res_multi)
#' }

extract_concentrations <- function(res_ASICS){
  concentration_list <- list()

  # Name of each sample
  pos <- 1
  for(i in 1:length(res_ASICS)){
    if(class(res_ASICS[[i]]) != "try-error"){
      concentration_list[[pos]] <- res_ASICS[[i]]@present_metabolites
      colnames(concentration_list[[pos]])[2] <- names(res_ASICS)[i]
      pos <- pos + 1
    }

  }

  #Get a data-frame from a list
  concentration <- plyr::join_all(concentration_list, by = "Name",
                                  type = "full")
  colnames(concentration)[1] <- "metabolites"

  return(concentration)
}

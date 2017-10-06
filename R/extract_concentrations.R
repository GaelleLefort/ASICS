#' Extract concentrations
#'
#' Combine results of the multiple file ASICS function to obtain quantified
#' relative concentrations of metabolites for each spectrum in one dataset
#'
#' @param res_ASICS result of the \code{\link{ASICS_multiFiles}} function
#' @return A data frame containing relative concentrations of identified
#' metabolites for each spectrum
#' @export
#' @examples
#' \dontshow{
#' lib_file <- system.file("extdata", "library_for_examples.rda",
#'                         package = "ASICS")
#' cur_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' res_multi <- ASICS_multiFiles(name.dir = cur_path,
#'                               exclusion.areas = to_exclude,
#'                               nb.iter.signif = 10, which.spectra = 2,
#'                               library.metabolites = lib_file, ncores = 1)
#' quantification <- extract_concentrations(res_multi)
#' }
#' \dontrun{
#' cur_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' res_multi <- ASICS_multiFiles(name.dir = cur_path,
#'                               exclusion.areas = to_exclude,
#'                               ncores = 2)
#'
#' # extract relative concentrations
#' quantification <- extract_concentrations(res_multi)
#' }

extract_concentrations <- function(res_ASICS){
  concentration_list <- list()

  # sample name
  pos <- 1
  for(i in 1:length(res_ASICS)){
    if(class(res_ASICS[[i]]) != "try-error"){
      concentration_list[[pos]] <- res_ASICS[[i]]@present_metabolites
      colnames(concentration_list[[pos]])[2] <- names(res_ASICS)[i]
      pos <- pos + 1
    }

  }

  # get a data-frame from a list
  concentration <- plyr::join_all(concentration_list, by = "Name",
                                  type = "full")
  colnames(concentration)[1] <- "metabolites"
  concentration[is.na(concentration)] <- 0


  return(concentration)
}

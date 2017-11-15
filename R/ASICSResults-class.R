#' Class \code{ASICSResults}
#'
#' Objects of class \code{ASICSResults} represent results of ASICS for a set of
#' spectra.
#' It an extension of class \linkS4class{Spectra} with an additional slots for
#' quantification results.
#'
#'
#' @name ASICSResults-class
#' @exportClass ASICSResults
#'
#' @slot sample.name Character vector of sample names.
#' @slot ppm.grid Numeric vector of a unique grid (definition domain) for all
#' spectra (in ppm).
#' @slot spectra Numeric matrix of original spectra in column.
#' Columns are in the same order as for \code{sample.name} and rows correspond
#' to points of \code{ppm.grid}.
#' @slot recomposed.spectra Numeric matrix of recomposed spectra (in column)
#' with estimated concentrations. Columns are in the same order as for
#' \code{sample.name} and rows correspond to points of \code{ppm.grid}.
#' @slot quantification Data-frame with identified metabolites and their
#' relative concentrations.
#' @slot deformed.library a list of \linkS4class{PureLibrary} containing the
#' deformed libraries of each sample.
#'
#' @include Spectra-class.R
#' @seealso \linkS4class{Spectra}

setClass(
  Class = "ASICSResults",
  slots = list(
    recomposed.spectra = "matrix",
    quantification = "data.frame",
    deformed.library = "list"
  ),
  contains = "Spectra"
)



#### Accessors

setGeneric("get_recomposed_spectra",
           function(object) standardGeneric("get_recomposed_spectra")
)
setGeneric("get_quantification",
           function(object) standardGeneric("get_quantification")
)
setGeneric("get_deformed_library",
           function(object) standardGeneric("get_deformed_library")
)


#' @export
#' @describeIn ASICSResults extract recomposed spectra matrix of the
#' \code{ASICSResults} object.
setMethod("get_recomposed_spectra", "ASICSResults",
          function(object) return(object@recomposed.spectra)
)

#' @export
#' @describeIn ASICSResults extract a data-frame with metabolite quantifications
#' computed with ASICS for each sample of the \code{ASICSResults} object.
setMethod("get_quantification", "ASICSResults",
          function(object) return(object@quantification)
)

#' @describeIn ASICSResults extract a list of all deformed libraries of the
#' \code{ASICSResults} object.
#' @export
setMethod("get_deformed_library", "ASICSResults",
          function(object) return(object@deformed.library)
)





#### Basic methods

#' @describeIn ASICSResults show a summary of the \code{ASICSResults} object.
#' @param object an object of class \code{ASICSResults}
#' @aliases show.ASICSResults
#' @importFrom utils head
#' @export
setMethod(
  f = "show",
  signature = "ASICSResults",
  definition = function(object){
    cat("An object of class", class(object), "\n")
    cat("It contains", length(object@sample.name), "spectra of",
        length(object@ppm.grid), "points. \n\n")
    cat("ASICS results:", nrow(object@quantification),
        "metabolites are identified for this set of spectra. \n")
    cat("Most concentrated metabolites are \n")
    print(head(object@quantification$Metabolites))
  }
)




#' @describeIn ASICSResults extract some samples from a \code{ASICSResults}
#' object.
#' @param i indices specifying samples to extract
#' @aliases [.ASICSResults
#' @export
setMethod(
  f = "[",
  signature(x = "ASICSResults", i = "ANY"),
  function (x, i){
    return(new("ASICSResults",
               sample.name = x@sample.name[i],
               ppm.grid = x@ppm.grid,
               spectra = x@spectra[, i],
               recomposed.spectra = x@recomposed.spectra[, i],
               quantification = x@quantification[, i + 1],
               deformed.library = x@deformed.library[[i]]))
  }
)


#' @describeIn ASICSResults combine \code{ASICSResults} objects.
#' @aliases c.ASICSResults
#' @importFrom plyr join_all
#' @export
setMethod(
  "c",
  signature(x = "ASICSResults"),
  function(x, ...){
    elements <- list(x, ...)

    # first grid for all objects
    for(i in 2:length(elements)){
      if(!any(elements[[1]]@ppm.grid == elements[[i]]@ppm.grid)){
        elements[[i]]@spectra <- apply(elements[[i]]@spectra, 2, change_grid,
                                       elements[[i]]@ppm.grid,
                                       elements[[1]]@ppm.grid)
      }
    }

    # merge quantification
    extract_quantif_list <- lapply(elements, get_quantification)

    all_quantification <- join_all(extract_quantif_list,
                                   by = "Metabolites", type = "full")

    return(new("ASICSResults",
               sample.name = do.call("c", lapply(elements, get_sample_name)),
               ppm.grid = x@ppm.grid,
               spectra = do.call("cbind", lapply(elements, get_spectra)),
               recomposed.spectra =
                 do.call("cbind", lapply(elements, get_recomposed_spectra)),
               quantification = all_quantification,
               deformed.library =
                 do.call("c", lapply(elements, get_deformed_library))))
  }
)



#### Other methods

#' @aliases diagnostics.ASICSResults
#' @export
diagnostics <- function(object){
  NULL
}













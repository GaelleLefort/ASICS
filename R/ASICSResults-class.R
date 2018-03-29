#' Class \linkS4class{ASICSResults}
#'
#' Objects of class \linkS4class{ASICSResults} contains results of ASICS
#' quantification method for a set of spectra. This object is an extension of
#' the class \linkS4class{Spectra}, with additional slots for quantification
#' results, reconstructed spectra and deformed library.
#'
#'
#' @name ASICSResults-class
#' @exportClass ASICSResults
#'
#' @slot sample.name Character vector of sample names.
#' @slot ppm.grid Numeric vector of a unique grid (definition domain) for all
#' spectra (in ppm).
#' @slot spectra Numeric matrix of original spectra. Columns contain the spectra
#' and are in the same order than \code{sample.name}. Rows correspond to points
#' of \code{ppm.grid}.
#' @slot reconstructed.spectra Numeric matrix of reconstructed spectra (in
#' columns) with estimated concentrations. Columns are in the same order than
#' \code{sample.name} and rows correspond to points of \code{ppm.grid}.
#' @slot quantification Data-frame with identified metabolites and their
#' relative concentrations.
#' @slot deformed.library A data frame containing the deformed library of each
#' sample.
#'
#' @section Methods:
#'   Multiple methods can be applied to \linkS4class{Spectra} objects.
#'   \itemize{
#'     \item As usual for S4 object, show and summary methods are available, see
#'     \link[=summary-methods]{Object summary}
#'     \item All slots have an accessor \code{get_slot name}, see
#'     \link[=accessors-methods]{Accessors}
#'     \item Two objects can be combined or a subset can be extracted, see
#'     \link[=combineAndSubset-methods]{Combine and subset methods}
#'     \item All spectra contained in an object can be represented in a plot,
#'     see \link[=visualisation-methods-spectra]{Visualisation methods}
#'   }
#'
#' @include Spectra-class.R
#' @seealso \linkS4class{Spectra}

setClass(
  Class = "ASICSResults",
  slots = list(
    reconstructed.spectra = "matrix",
    quantification = "data.frame",
    deformed.library = "data.frame"
  ),
  contains = "Spectra"
)



#### Accessors

setGeneric("getReconstructedSpectra",
           function(object) standardGeneric("getReconstructedSpectra")
)
setGeneric("getQuantification",
           function(object) standardGeneric("getQuantification")
)
setGeneric("getDeformedLibrary",
           function(object) standardGeneric("getDeformedLibrary")
)


#' @export
#' @aliases getReconstructedSpectra
#' @rdname accessors-methods
setMethod("getReconstructedSpectra", "ASICSResults",
          function(object) return(object@reconstructed.spectra)
)

#' @export
#' @aliases getQuantification
#' @rdname accessors-methods
setMethod("getQuantification", "ASICSResults",
          function(object) return(object@quantification)
)


#' @export
#' @aliases getDeformedLibrary
#' @rdname accessors-methods
setMethod("getDeformedLibrary", "ASICSResults",
          function(object) return(object@deformed.library)
)





#### Basic methods

#' @rdname summary-methods
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
    cat("ASICS results: \n", nrow(object@quantification),
        "metabolites are identified for this set of spectra. \n")
    cat("Most concentrated metabolites are:",
        paste(head(rownames(object@quantification)), collapse = ", "))
  }
)

#' @rdname summary-methods
#' @export
#' @aliases dim.Spectra
setMethod(f = "dim", signature(x = "ASICSResults"),
          function(x) return(c(nrow(x@quantification), length(x@sample.name)))
)


#' @rdname combineAndSubset-methods
#' @aliases [.ASICSResults
#' @export
setMethod(
  f = "[",
  signature(x = "ASICSResults", i = "ANY"),
  function (x, i){

    quantification <- as.data.frame(x@quantification[, i])
    colnames(quantification) <- x@sample.name[i]
    rownames(quantification) <- rownames(x@quantification)

    return(new("ASICSResults",
               sample.name = x@sample.name[i],
               ppm.grid = x@ppm.grid,
               spectra = as.matrix(x@spectra[, i]),
               reconstructed.spectra = as.matrix(x@reconstructed.spectra[, i]),
               quantification = quantification,
               deformed.library = x@deformed.library[x@deformed.library$sample
                                                     %in% x@sample.name[i], ]))
  }
)


#' @rdname combineAndSubset-methods
#' @aliases c.ASICSResults
#' @importFrom plyr join_all
#' @export
setMethod(
  "c",
  signature(x = "ASICSResults"),
  function(x, ...){
    elements <- list(x, ...)

    # first grid for all objects
    if(length(elements) > 1){
      for(i in 2:length(elements)){
        if(!any(elements[[1]]@ppm.grid == elements[[i]]@ppm.grid)){
          elements[[i]]@spectra <- apply(elements[[i]]@spectra, 2, .changeGrid,
                                         elements[[i]]@ppm.grid,
                                         elements[[1]]@ppm.grid)
        }
      }
    }

    # merge quantification
    extract_quantif_list <- lapply(elements, getQuantification)
    extract_quantif_list <- lapply(extract_quantif_list,
                                   function(x) cbind(id = rownames(x), x))

    all_quantification <- join_all(extract_quantif_list,
                                   by = "id", type = "full")
    rownames(all_quantification) <- all_quantification$id
    all_quantification$id <- NULL

    all_quantification[is.na(all_quantification)] <- 0

    return(new("ASICSResults",
               sample.name = do.call("c", lapply(elements, getSampleName)),
               ppm.grid = x@ppm.grid,
               spectra = do.call("cbind", lapply(elements, getSpectra)),
               reconstructed.spectra =
                 do.call("cbind", lapply(elements, getReconstructedSpectra)),
               quantification = all_quantification,
               deformed.library =
                 do.call("rbind", lapply(elements, getDeformedLibrary))))
  }
)


#' @aliases plot.ASICSResults
#' @param idx Index of the spectrum to plot. Default to 1.
#' @param pure.library Pure library used for the quantification. Default to
#' \code{NULL} (in which case, the library included in the package is used).
#' @param add.metab Name of one metabolite to add to the plot. Default to
#' \code{NULL} (in which case, no pure spectrum added to the plot).
#' @export
#' @rdname visualisation-methods-spectra
setMethod(
  f = "plot",
  signature = "ASICSResults",
  definition = function(x, y, idx = 1, xlim = c(0.5, 10), ylim = NULL,
                        pure.library = NULL, add.metab = NULL, ...) {
    return(.plotSpectrum(x, idx = idx, xlim = xlim, ylim = ylim,
                         pure.library = pure.library, add.metab = add.metab))
  }
)



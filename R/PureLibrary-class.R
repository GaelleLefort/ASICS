#' Class \code{PureLibrary}
#'
#' Objects of class \code{PureLibrary} contain a set of pure metabolite NMR
#' spectra, used as a reference for the quantification. This class is an 
#' extension of the class \linkS4class{Spectra}, with an additional slot
#' (number of protons for each metabolite) needed for spectrum quantification.
#'
#'
#' @slot nb.protons Numeric vector of the number of protons of each pure
#' metabolite spectra.
#'
#' @section Methods:
#'   Multiple methods can be applied on \linkS4class{PureLibrary} objects.
#'   \itemize{
#'     \item As usual for S4 object, show and summary methods are available, see
#'     \link[=summary-methods]{Object summary}
#'     \item All slots have an accessor \code{get_slot name}, see
#'     \link[=accessors-methods]{Accessors}
#'     \item Two objects can be combined or a subset can be extracted, see
#'     \link[=combineAndSubset-methods]{Combine and subset methods}
#'     \item All spectra contained in an object can be represented in a plot, 
#'     see \link[=visualization-methods-spectra]{Visualization methods}
#'   }
#'
#' @seealso \linkS4class{Spectra}
#'
#' @name PureLibrary-class
#' @exportClass PureLibrary
#' @include Spectra-class.R
#'

setClass("PureLibrary",
  slots = list(
    nb.protons = "numeric"
  ),
  contains = "Spectra"
)

setValidity(
  Class = "PureLibrary",
  function(object) {
    errors <- character()

    if (nrow(object@spectra) != length(object@ppm.grid)) {
      msg <- paste("Length of ppm grid and number of rows of spectra matrix",
                   "must be identical.")
      errors <- c(errors, msg)
    }
    if (ncol(object@spectra) != length(object@sample.name)) {
      msg <- paste("Numbers of sample names and rows of spectra matrix must be",
                   "identical.")
      errors <- c(errors, msg)
    }
    if (length(object@nb.protons) != length(object@sample.name)) {
      msg <- paste("Number of sample names and length of number of protons",
                   "vector of spectra matrix must be identical.")
      errors <- c(errors, msg)
    }
    if(any(duplicated(object@sample.name))){
      msg <- paste("Sample names need to be unique.")
      errors <- c(errors, msg)
    }

    if (length(errors) == 0) TRUE else errors
  }
)


#### Accessors

setGeneric("getNbProtons",
           function(object) standardGeneric("getNbProtons")
)


#' @rdname accessors-methods
#' @aliases getNbProtons
#' @export
setMethod("getNbProtons", "PureLibrary",
          function(object) return(object@nb.protons)
)


#### Basic methods

#' @rdname combineAndSubset-methods
#' @aliases [.PureLibrary
#' @export
setMethod(
  f = "[",
  signature(x = "PureLibrary", i = "ANY"),
  function (x, i) {
    return(new("PureLibrary",
               sample.name = x@sample.name[i],
               ppm.grid = x@ppm.grid,
               spectra = as.matrix(x@spectra[, i]),
               nb.protons = x@nb.protons[i]))
  }
)


#' @rdname combineAndSubset-methods
#' @aliases c.PureLibrary
#' @export
setMethod(
  "c",
  signature(x = "PureLibrary"),
  function(x, ...) {
    elements <- list(x, ...)

    if(any(duplicated(do.call("c", lapply(elements, getSampleName))))){
      stop("Sample names need to be unique.")
    }

    # first grid for all objects
    for (i in 2:length(elements)) {
      if (!any(elements[[1]]@ppm.grid == elements[[i]]@ppm.grid)) {
        elements[[i]]@spectra <- apply(elements[[i]]@spectra, 2, .changeGrid,
                                       elements[[i]]@ppm.grid,
                                       elements[[1]]@ppm.grid)
      }
    }

    return(new("PureLibrary",
               sample.name = do.call("c", lapply(elements, getSampleName)),
               ppm.grid = x@ppm.grid,
               spectra = do.call("cbind", lapply(elements, getSpectra)),
               nb.protons = do.call("c", lapply(elements, getNbProtons))))
  }
)

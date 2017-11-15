#' Class \code{PureLibrary}
#'
#' Objects of class \code{PureLibrary} represent a set of pure metabolite NMR
#' spectra.
#' It an extension of class \linkS4class{Spectra} with an additional slot needed for
#' spectrum quantification.
#'
#'
#' @slot nb.protons Numeric vector of the number of protons of each pure
#' metabolite spectra.
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
      msg <- paste("Number of sample names and rows of spectra matrix must be",
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

setGeneric("get_nb_protons",
           function(object) standardGeneric("get_nb_protons")
)


#' @export
#' @param object an object of class \code{PureLibrary}
#' @describeIn PureLibrary extract number of protons of the \code{PureLibrary}
#' object.
setMethod("get_nb_protons", "PureLibrary",
          function(object) return(object@nb.protons)
)


#### Basic methods

#' @describeIn PureLibrary extract some metabolites from a \code{PureLibrary}
#' object.
#' @param x an object of class \code{PureLibrary}
#' @param i indices specifying metabolites to extract
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


#' @describeIn PureLibrary combine \code{PureLibrary} objects.
#' @aliases c.PureLibrary
#' @export
setMethod(
  "c",
  signature(x = "PureLibrary"),
  function(x, ...) {
    elements <- list(x, ...)

    if(any(duplicated(do.call("c", lapply(elements, get_sample_name))))){
      stop("Sample names need to be unique.")
    }

    # first grid for all objects
    for (i in 2:length(elements)) {
      if (!any(elements[[1]]@ppm.grid == elements[[i]]@ppm.grid)) {
        elements[[i]]@spectra <- apply(elements[[i]]@spectra, 2, change_grid,
                                       elements[[i]]@ppm.grid,
                                       elements[[1]]@ppm.grid)
      }
    }

    return(new("PureLibrary",
               sample.name = do.call("c", lapply(elements, get_sample_name)),
               ppm.grid = x@ppm.grid,
               spectra = do.call("cbind", lapply(elements, get_spectra)),
               nb.protons = do.call("c", lapply(elements, get_nb_protons))))
  }
)

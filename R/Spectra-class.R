#' Class \linkS4class{Spectra}
#'
#' Objects of class \linkS4class{Spectra} contain a set of NMR spectra.
#' It includes preprocessed spectra and can be created with the function
#' \code{\link{createSpectra}}.
#'
#' @name Spectra-class
#' @import Matrix
#' @exportClass Spectra
#'
#'
#' @slot sample.name Character vector of sample names.
#' @slot ppm.grid Numeric vector of a unique grid (definition domain) for all
#' spectra (in ppm).
#' @slot spectra Numeric matrix with all spectra in columns. Columns must be in
#' the same order as for \code{sample.name} and rows correspond to points of
#' \code{ppm.grid}.
#' @slot norm.method Character specifying the normalisation method to use on
#' spectra
#' @slot norm.params List containing normalisation parameteres (see
#' \code{\link{normaliseSpectra}} for details).
#'
#' @section Methods:
#'   Multiple methods can be applied on \linkS4class{Spectra} objects.
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
setClass(
  Class = "Spectra",
  slots = list(
    sample.name = "character",
    ppm.grid = "numeric",
    spectra = "generalMatrix",
    norm.method = "character",
    norm.params = "list"
  )
)

setValidity(
  Class = "Spectra",
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
    if(any(duplicated(object@sample.name))){
      msg <- paste("Sample names need to be unique.")
      errors <- c(errors, msg)
    }

    if (length(errors) == 0) TRUE else errors
  }
)


#### Accessors

setGeneric("getSampleName",
           function(object) standardGeneric("getSampleName")
)
setGeneric("getPpmGrid",
           function(object) standardGeneric("getPpmGrid")
)
setGeneric("getSpectra",
           function(object) standardGeneric("getSpectra")
)
setGeneric("getNormMethod",
           function(object) standardGeneric("getNormMethod")
)
setGeneric("getNormParams",
           function(object) standardGeneric("getNormParams")
)

#' Accessors
#'
#' List of available accessors for each slot of all S4 classes present in the
#' package.
#'
#' @name accessors-methods
#' @param object An object of class \linkS4class{Spectra},
#' \linkS4class{PureLibrary}, \linkS4class{ASICSResults} or
#' \linkS4class{AnalysisResults}.
#'
#' @return The wanted accessor
#'
#' @examples
#' # Import data and create object
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_obj <- createSpectra(spectra_data)
#'
#' # Sample names
#' getSampleName(spectra_obj)
#' # Spectra
#' getSpectra(spectra_obj)
NULL

#' @export
#' @aliases getSampleName
#' @rdname accessors-methods
setMethod("getSampleName", "Spectra",
          function(object) return(object@sample.name)
)

#' @export
#' @aliases getPpmGrid
#' @rdname accessors-methods
setMethod("getPpmGrid", "Spectra",
          function(object) return(object@ppm.grid)
)

#' @rdname accessors-methods
#' @aliases getSpectra
#' @export
setMethod("getSpectra", "Spectra",
          function(object) {
            spectra <- object@spectra
            rownames(spectra) <- object@ppm.grid
            colnames(spectra) <- object@sample.name
            return(spectra)
          }
)

#' @export
#' @aliases getNormMethod
#' @rdname accessors-methods
setMethod("getNormMethod", "Spectra",
          function(object) return(object@norm.method)
)

#' @export
#' @aliases getNormParams
#' @rdname accessors-methods
setMethod("getNormParams", "Spectra",
          function(object) return(object@norm.params)
)

#### Basic methods

#' Summary methods
#'
#' Methods available to summarize the various S4 objects of ASICS package.
#'
#' @name summary-methods
#' @param object An object of class \linkS4class{Spectra},
#' \linkS4class{PureLibrary}, \linkS4class{ASICSResults} or
#' \linkS4class{AnalysisResults}.
#'
#' @return A summary of the object, its length or its dimensions.
#'
#' @examples
#' # Import data and create object
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_obj <- createSpectra(spectra_data)
#'
#' # Summary
#' summary(spectra_obj)
#' # Length
#' length(spectra_obj)
#' # Dimensions
#' dim(spectra_obj)
NULL

#' @aliases show.Spectra
#' @importFrom methods show
#' @export
#' @rdname summary-methods
setMethod(
  f = "show",
  signature = "Spectra",
  definition = function(object){
    cat("An object of class", class(object), "\n")
    cat("It contains", length(object@sample.name), "spectra of",
        length(object@ppm.grid), "points. \n")
  }
)

#' @rdname summary-methods
#' @export
#' @aliases summary.Spectra
setMethod(
  f = "summary",
  signature = "Spectra",
  definition = function(object) object
)

#' @rdname summary-methods
#' @param x An object of class \linkS4class{Spectra},
#' \linkS4class{PureLibrary} or \linkS4class{ASICSResults}.
#' @export
#' @aliases length.Spectra
setMethod(f = "length", signature(x = "Spectra"),
          function(x) return(length(x@sample.name))
)

#' @rdname summary-methods
#' @export
#' @aliases dim.Spectra
setMethod(f = "dim", signature(x = "Spectra"),
          function(x) return(c(nrow(x@spectra), length(x@sample.name)))
)


#' Combine or subset functions
#'
#' Methods available to combine multiple objects or to extract a subset of one
#' object in ASICS package.
#'
#' @name combineAndSubset-methods
#' @param x An object of class \linkS4class{Spectra},
#' \linkS4class{PureLibrary} or \linkS4class{ASICSResults}.
#' @param i vector of indices specifying which elements to extract
#' @param ... objects to be concatenated
#'
#' @return A \linkS4class{Spectra} object containing a part of the original
#' object or combining other \linkS4class{Spectra} objects
#'
#' @examples
#' # Import data and create object
#' current_path <- file.path(system.file("extdata", package = "ASICS"),
#'                           "spectra_example.txt")
#' spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
#' spectra_obj <- createSpectra(spectra_data)
#'
#' # Extract the first sample
#' spectra_obj[1]
NULL

#' @rdname combineAndSubset-methods
#' @aliases [.Spectra
#' @export
setMethod(
  f = "[",
  signature(x = "Spectra", i = "ANY"),
  function(x, i) {
    return(new("Spectra",
               sample.name = x@sample.name[i],
               ppm.grid = x@ppm.grid,
               spectra = Matrix(x@spectra[ ,i]),
               norm.method = x@norm.method,
               norm.params = x@norm.params))
  }
)


#' @rdname combineAndSubset-methods
#' @export
#' @aliases c.Spectra
setMethod(
  "c",
  signature(x = "Spectra"),
  function(x, ...) {
    elements <- list(x, ...)

    if(any(duplicated(do.call("c", lapply(elements, getSampleName))))){
      stop("Sample names need to be unique.")
    }

    if (!all(do.call("c", lapply(elements, getNormMethod)) ==
             elements[[1]]@norm.method)) {
      warning(paste("Elements not have the same normalisation method, the",
                    "first one is used"))
    }

    # first grid for all objects
    for(i in 2:length(elements)){
      if(!any(elements[[1]]@ppm.grid == elements[[i]]@ppm.grid)){
        elements[[i]]@spectra <- apply(elements[[i]]@spectra, 2, .changeGrid,
                                       elements[[i]]@ppm.grid,
                                       elements[[1]]@ppm.grid)
      }
    }

    return(new("Spectra",
               sample.name = do.call("c", lapply(elements, getSampleName)),
               ppm.grid = x@ppm.grid,
               spectra = do.call("cbind", lapply(elements, getSpectra)),
               norm.method = x@norm.method,
               norm.params = x@norm.params))
  }
)



#' Visualisation methods
#'
#' Methods available to plot one object in ASICS package.
#'
#' @name visualisation-methods-spectra
#' @param x An object of class \linkS4class{Spectra},
#' \linkS4class{PureLibrary} or \linkS4class{ASICSResults}.
#' @param xlim,ylim Boundaries for x and y, respectively.
#' @param y Currently not used.
#' @param ... Currently not used.
#'
#' @return
#' \itemize{
#' \item A \code{\link[ggplot2]{ggplot}} plot of all spectra (or of a subset) on
#' the same figure for \linkS4class{Spectra} and \linkS4class{PureLibrary}
#' objects.
#' \item A \code{\link[ggplot2]{ggplot}} plot of original and reconstructed
#' spectra of one sample in the same figure for \linkS4class{ASICSResults}
#' object. In addition, one pure metabolite spectrum (as provided in the
#' reference library) and the deformed one can be superimposed to the plot.
#' }
#'
#' @examples
#' # Import data and create object
#' current_path <- system.file("extdata", "example_spectra", package = "ASICS")
#' spectra_data <- importSpectraBruker(current_path)
#' spectra_obj <- createSpectra(spectra_data)
#' spectra_obj <- createSpectra(spectra_data)
#'
#' # Plot the spectra
#' plot(spectra_obj)
NULL


#' @aliases plot.Spectra
#'
#' @importFrom stats reshape
#' @import ggplot2
#' @export
#' @rdname visualisation-methods-spectra
setMethod(
  f = "plot",
  signature = "Spectra",
  definition = function(x, y, xlim = c(0.5, 10), ylim = NULL, ...) {

    # reshape data
    data_to_plot <- reshape(as.data.frame(as.matrix(x@spectra)),
                            varying = seq_len(ncol(x@spectra)),
                            times = x@sample.name,
                            ids = as.numeric(x@ppm.grid),
                            v.names = "intensity", timevar = "sample_name",
                            idvar = "ppm_grid", direction = "long",
                            new.row.names = seq_len(ncol(x@spectra) *
                                                 nrow(x@spectra)))

    # set y boundaries if ylim is NULL
    if (is.null(ylim)) ylim <- c(0, max(data_to_plot$intensity))

    # plot
    graph <- ggplot(data_to_plot, aes_string(x = "ppm_grid", y = "intensity",
                                             color = "sample_name")) +
      geom_line(na.rm = TRUE) + theme_bw() +
      labs(x = "Chemical shift (ppm)", y = "Intensity",
           color = "Spectrum names") +
      scale_x_reverse(limits = c(xlim[2], xlim[1])) +
      ylim(c(ylim[1], ylim[2]))

    # remove legend if more than 10 spectra
    if (length(x@sample.name) > 10) {
      graph <- graph + guides(color = FALSE)
    }

    return(graph)
  }
)






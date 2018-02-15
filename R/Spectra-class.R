#' Class \code{Spectra}
#'
#' Objects of class \code{Spectra} represent a set of NMR spectra of one study.
#' It contains preprocessed spectra and can be created with the function
#' \code{\link{create_spectra}}.
#'
#' @name Spectra-class
#' @exportClass Spectra
#'
#' @slot sample.name Character vector of sample names.
#' @slot ppm.grid Numeric vector of a unique grid (definition domain) for all
#' spectra (in ppm).
#' @slot spectra Numeric matrix with all spectra in columns. Columns must be in
#' the same order as for \code{sample.name} and rows correspond to points of
#' \code{ppm.grid}.
#'
setClass(
  Class = "Spectra",
  slots = list(
    sample.name = "character",
    ppm.grid = "numeric",
    spectra = "matrix"
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
      msg <- paste("Number of sample names and rows of spectra matrix must be",
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

setGeneric("get_sample_name",
           function(object) standardGeneric("get_sample_name")
)
setGeneric("get_ppm_grid",
           function(object) standardGeneric("get_ppm_grid")
)
setGeneric("get_spectra",
           function(object) standardGeneric("get_spectra")
)


#' @export
#' @describeIn Spectra extract sample names of the \code{Spectra} object.
setMethod("get_sample_name", "Spectra",
          function(object) return(object@sample.name)
)

#' @export
#' @describeIn Spectra extract p.p.m grid of the \code{Spectra} object.
setMethod("get_ppm_grid", "Spectra",
          function(object) return(object@ppm.grid)
)

#' @describeIn Spectra extract spectra matrix of the \code{Spectra} object.
#' @export
setMethod("get_spectra", "Spectra",
          function(object) {
            spectra <- object@spectra
            rownames(spectra) <- object@ppm.grid
            colnames(spectra) <- object@sample.name
            return(object@spectra)
          }
)



#### Basic methods

#' @describeIn Spectra show a summary of the \code{Spectra} object.
#' @param object an object of class \code{Spectra}
#' @aliases show.Spectra
#' @export
setMethod(
  f = "show",
  signature = "Spectra",
  definition = function(object){
    cat("An object of class", class(object), "\n")
    cat("It contains", length(object@sample.name), "spectra of",
        length(object@ppm.grid), "points. \n")
  }
)

#' @describeIn Spectra show a summary of the \code{Spectra} object.
#' @aliases summary.Spectra
#' @export
setMethod(
  f = "summary",
  signature = "Spectra",
  definition = function(object) object
)


#' @describeIn Spectra extract some samples from a \code{Spectra} object.
#' @param i indices specifying elements to extract
#' @aliases [.Spectra
#' @export
setMethod(
  f = "[",
  signature(x = "Spectra", i = "ANY"),
  function(x, i) {
    return(new("Spectra",
               sample.name = x@sample.name[i],
               ppm.grid = x@ppm.grid,
               spectra = as.matrix(x@spectra[, i])))
  }
)



#' @describeIn Spectra number of samples in a \code{Spectra} object.
#' @aliases length.Spectra
#' @export
setMethod(f = "length", signature(x = "Spectra"),
  function(x) return(length(x@sample.name))
)


#' @describeIn Spectra combine \code{Spectra} objects.
#' @aliases c.Spectra
#' @export
setMethod(
  "c",
  signature(x = "Spectra"),
  function(x, ...) {
    elements <- list(x, ...)

    if(any(duplicated(do.call("c", lapply(elements, get_sample_name))))){
      stop("Sample names need to be unique.")
    }

    # first grid for all objects
    for(i in 2:length(elements)){
      if(!any(elements[[1]]@ppm.grid == elements[[i]]@ppm.grid)){
        elements[[i]]@spectra <- apply(elements[[i]]@spectra, 2, change_grid,
                                       elements[[i]]@ppm.grid,
                                       elements[[1]]@ppm.grid)
      }
    }

    return(new("Spectra",
               sample.name = do.call("c", lapply(elements, get_sample_name)),
               ppm.grid = x@ppm.grid,
               spectra = do.call("cbind", lapply(elements, get_spectra))))
  }
)



#' @describeIn Spectra plot all spectra (or a subset) on the same figure.
#' Legend only appears if there is less than ten spectra to plot.
#' @aliases plot.Spectra
#'
#' @param x an object of class \code{Spectra}
#' @param xlim,ylim boundaries for x and y, respectively
#' @param y currently not used
#' @param ... currently not used
#'
#' @importFrom stats reshape
#' @importFrom ggplot2 ggplot aes_string geom_line theme_bw labs
#' @importFrom ggplot2 scale_x_reverse ylim guides
#' @export
setMethod(
  f = "plot",
  signature = "Spectra",
  definition = function(x, y, xlim = c(0.5, 10), ylim = NULL, ...) {

    # reshape data
    data_to_plot <- reshape(as.data.frame(x@spectra),
                            varying = 1:ncol(x@spectra),
                            times = x@sample.name,
                            ids = as.numeric(x@ppm.grid),
                            v.names = "intensity", timevar = "sample_name",
                            idvar = "ppm_grid", direction = "long",
                            new.row.names = 1:(ncol(x@spectra) *
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






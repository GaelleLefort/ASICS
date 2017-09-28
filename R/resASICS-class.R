#' An S4 class to represent results of ASICS.
#' @name resASICS-class
#' @exportClass resASICS
#'
#' @slot original_mixture original spectrum
#' @slot reconstituted_mixture reconstituted spectrum with estimated
#'   concentrations
#' @slot ppm_grid grid of the spectrum (in p.p.m.)
#' @slot present_metabolites a data frame with identified metabolites and their
#'   relative concentrations
#' @slot non_identified_metabolites a data frame with non-identified metabolites
#'   and their identification thresholds
#'
#' @note Slots can be accessed by accessor functions with the same name
#' (see \link{resASICS-getters})

setClass(
  Class = "resASICS",
  slots = list(
    original_mixture = "numeric",
    reconstituted_mixture = "numeric",
    ppm_grid = "numeric",
    present_metabolites = "data.frame",
    non_identified_metabolites = "data.frame"
  )
)


#' S4 methods to represent results of ASICS.
#' @name resASICS-methods
NULL


#' @param object an object of class resASICS
#' @param ... not used
#' @aliases summary.resASICS
#' @rdname resASICS-methods
#' @export
setMethod(
  f = "summary",
  signature = "resASICS",
  definition = function(object, ...){
    print(object)
  }
)

#' @rdname resASICS-methods
#' @aliases print.resASICS
#' @importFrom utils head
#' @export
setMethod(
  f = "print",
  signature = "resASICS",
  definition = function(x){
    cat("Number of identified metabolites:", nrow(x@present_metabolites), "\n\n")

    cat("Most concentrated metabolites: \n")
    print(head(x@present_metabolites))
  }
)


#' @param x an object of class resASICS
#' @param y not used
#' @param xmin,xmax x minimum and maximum limits
#' @param ymin,ymax y minimum and maximum limits
#' @param add_metab name of one metabolite to add to the plot
#' @aliases plot.resASICS
#' @export
#' @rdname resASICS-methods
setMethod(
  f = "plot",
  signature = "resASICS",
  definition = function(x, y, ..., xmin = 0, xmax = 10, ymin = 0, ymax = NULL,
                        add_metab = NULL){
    plot_spectrum(x, xmin, xmax, ymin, ymax, add_metab)
  }
)






setGeneric("present_metabolites",
           function(object){standardGeneric("present_metabolites")}
)
setGeneric("non_identified_metabolites",
           function(object){standardGeneric("non_identified_metabolites")}
)
setGeneric("original_mixture",
           function(object){standardGeneric("original_mixture")}
)
setGeneric("reconstituted_mixture",
           function(object){standardGeneric("reconstituted_mixture")}
)
setGeneric("ppm_grid",
           function(object){standardGeneric("ppm_grid")}
)


#' S4 methods to represent results of ASICS.
#' @name resASICS-getters
NULL

#' @param object an object of class resASICS
#' @return The respective slot from resASICS object.
#' @export
#' @aliases present_metabolites
#' @rdname resASICS-getters
setMethod("present_metabolites", "resASICS",
          function(object){
            return(object@present_metabolites)
          }
)


#' @export
#' @aliases non_identified_metabolites
#' @rdname resASICS-getters
setMethod("non_identified_metabolites", "resASICS",
          function(object){
            return(object@non_identified_metabolites)
          }
)


#' @export
#' @aliases original_mixture
#' @rdname resASICS-getters
setMethod("original_mixture", "resASICS",
          function(object){
            return(object@original_mixture)
          }
)


#' @export
#' @aliases reconstituted_mixture
#' @rdname resASICS-getters
setMethod("reconstituted_mixture", "resASICS",
          function(object){
            return(object@reconstituted_mixture)
          }
)


#' @export
#' @aliases ppm_grid
#' @rdname resASICS-getters
setMethod("ppm_grid", "resASICS",
          function(object){
            return(object@ppm_grid)
          }
)


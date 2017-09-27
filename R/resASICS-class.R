#' An S4 class to represent results of ASICS.
#' @name resASICS-class
#' @rdname resASICS-class
#' @exportClass resASICS
#'
#' @slot original_mixture original spectrum
#' @slot reconstituted_mixture reconstituted spectrum with estimated
#'   concentrations
#' @slot grid grid of the spectrum (in p.p.m.)
#' @slot present_metabolites a data frame with identified metabolites and their
#'   relative concentrations
#' @slot non_identified_metabolites a data frame with non-identified metabolites
#'   and their identification thresholds
setClass(
  Class = "resASICS",
  slots = list(
    original_mixture = "numeric",
    reconstituted_mixture = "numeric",
    grid = "numeric",
    present_metabolites = "data.frame",
    non_identified_metabolites = "data.frame"
  )
)


#' @param x an object of class resASICS
#' @param y not used
#' @param xmin,xmax x minimum and maximum limits
#' @param ymin,ymax y minimum and maximum limits
#' @param add_metab name of one metabolite to add to the plot
#' @rdname resASICS-class
#' @aliases plot.resASICS
#' @export
setMethod(
  f = "plot",
  signature = "resASICS",
  definition = function(x, y, ..., xmin = 0, xmax = 10, ymin = 0, ymax = NULL,
                        add_metab = NULL){
    plot_spectrum(x, xmin, xmax, ymin, ymax, add_metab)
  }
)

#' @param object an object of class resASICS
#' @param ... not used
#' @rdname resASICS-class
#' @aliases summary.resASICS
#' @importFrom utils head
#' @export
setMethod(
  f = "summary",
  signature = "resASICS",
  definition = function(object, ...){
    print(object)
  }
)

#' @rdname resASICS-class
#' @aliases print.resASICS
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

#' An S4 class to represent results of ASICS.
#' @name resASICS-class
#' @exportClass resASICS
#'
#' @slot original_mixture original spectrum
#' @slot reconstituted_mixture reconstituted spectrum with estimated
#' concentrations
#' @slot ppm_grid grid (definition domain) of the spectrum (in ppm)
#' @slot present_metabolites a data frame with identified metabolites and their
#' relative concentrations
#'
#' @note Slots can be accessed by accessor functions with the same name
#' (see \link{resASICS-getters}).
#'
#' @seealso \code{\link{ASICS}} \code{\link{resASICS-methods}}

setClass(
  Class = "resASICS",
  slots = list(
    original_mixture = "numeric",
    reconstituted_mixture = "numeric",
    ppm_grid = "numeric",
    present_metabolites = "data.frame"
  )
)


#' S4 methods to represent results of ASICS.
#' @name resASICS-methods
#' @examples
#' \dontshow{
#' lib_file <- system.file("extdata", "library_for_examples.rda",
#'                         package = "ASICS")
#' cur_path <- system.file("extdata", "example_spectra", "AG_faq_Beck01",
#'                         package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' result <- ASICS(path = cur_path, exclusion.areas = to_exclude,
#'                 nb.iter.signif = 10, library.metabolites = lib_file)
#'
#' #Results of ASICS
#' result
#' summary(result)
#' plot(result)
#'
#' }
#' \dontrun{
#' cur_path <- system.file("extdata", "example_spectra", "AG_faq_Beck01",
#'                         package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' result <- ASICS(path = cur_path, exclusion.areas = to_exclude)
#'
#' result
#' summary(result)
#' plot(result)
#'
#' }
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

#' @aliases show.resASICS
#' @rdname resASICS-methods
#' @importFrom methods show
#' @export
setMethod(
  f = "show",
  signature = "resASICS",
  definition = function(object){
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
#' @param xmin,xmax,ymin,ymax lower and upper bounds for x and y, respectively
#' @param add_metab name of one metabolite to add to the plot. Default to
#' \code{NULL} (no pure spectrum added to the plot)
#'
#' @aliases plot.resASICS
#'
#' @return plot the true and recomposed (as estimated by \code{\link{ASICS}})
#' spectra on one figure. In addition, one pure metabolite spectrum (as
#' provided in the reference library) can be superimposed to the plot.
#'
#' @seealso \code{\link{ASICS}} \code{\link{resASICS-class}}
#'
#' @export
#'
#' @rdname resASICS-methods
setMethod(
  f = "plot",
  signature = "resASICS",
  definition = function(x, y, xmin = 0.5, xmax = 10, ymin = 0, ymax = NULL,
                        add_metab = NULL){
    plot_spectrum(x, xmin, xmax, ymin, ymax, add_metab)
  }
)


setGeneric("present_metabolites",
           function(object){standardGeneric("present_metabolites")}
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
#' @examples
#' \dontshow{
#' lib_file <- system.file("extdata", "library_for_examples.rda",
#'                         package = "ASICS")
#' cur_path <- system.file("extdata", "example_spectra", "AG_faq_Beck01",
#'                         package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' result <- ASICS(path = cur_path, exclusion.areas = to_exclude,
#'                 nb.iter.signif = 10, library.metabolites = lib_file)
#'
#' #Identified metabolites
#' present_metabolites(result)
#'
#' }
#' \dontrun{
#' cur_path <- system.file("extdata", "example_spectra", "AG_faq_Beck01",
#'                         package = "ASICS")
#' to_exclude <- matrix(c(4.5,5.1,5.5,6.5), ncol = 2, byrow = TRUE)
#' result <- ASICS(path = cur_path, exclusion.areas = to_exclude)
#'
#' #Identified metabolites
#' present_metabolites(result)
#'
#' }
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


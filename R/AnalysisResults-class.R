#' Class \linkS4class{AnalysisResults}
#'
#' Objects of class \linkS4class{AnalysisResults} contains results of analyses
#' performed with the functions \code{\link{pca}}, \code{\link{oplsda}} and
#' \code{\link{kruskalWallis}}.
#'
#' @name AnalysisResults-class
#' @exportClass AnalysisResults
#'
#'
#' @slot type.analysis Name of the analysis (\emph{e.g.,} \code{"PCA"},
#' \code{"OPLS-DA"}, ...).
#' @slot type.data Type of data used for the analyses (\emph{e.g.,}
#' \code{"quantification"}, \code{"buckets"}...).
#' @slot dataset The object of type \code{\link{SummarizedExperiment}} used for
#' the analysis.
#' @slot results Results of the analysis. Can be a data frame for test results
#' or an object of class \code{\link{opls}} from \code{\link{ropls}} for PCA and
#' OPLS-DA.
#' @slot best.model Best model (only for OPLS-DA analyses).
#' @slot cv.error Cross validation error (only for OPLS-DA analyses).
#' @slot mean.by.group Data frame with means by group and a variable indicating
#' if there is a significant difference between groups for tests and if the VIP
#' associated to the variable is superior to the given threshold for OPLS-DA.
#'
#' @section Methods:
#'   Multiple methods can be applied on \linkS4class{AnalysisResults} objects.
#'   \itemize{
#'     \item As usual for S4 object, show and summary methods are available, see
#'     \link[=summary-methods]{Object summary}
#'     \item All slots have an accessor \code{get_slot name}, see
#'     \link[=accessors-methods]{Accessors}
#'     \item All results contained in an object can be represent in a plot, see
#'     \link[=visualisation-methods-analyses]{Visualisation methods}
#'   }
#'
setClass(
  Class = "AnalysisResults",
  slots = list(
    type.analysis = "character",
    type.data = "character",
    dataset = "SummarizedExperiment",
    results = "ANY",
    best.model = "list",
    cv.error = "numeric",
    mean.by.group = "data.frame"
  )
)

#### Accessors

setGeneric("getTypeAnalysis",
           function(object) standardGeneric("getTypeAnalysis")
)
setGeneric("getTypeData",
           function(object) standardGeneric("getTypeData")
)
setGeneric("getDataset",
           function(object) standardGeneric("getDataset")
)
setGeneric("getResults",
           function(object) standardGeneric("getResults")
)
setGeneric("getBestModel",
           function(object) standardGeneric("getBestModel")
)
setGeneric("getCVError",
           function(object) standardGeneric("getCVError")
)
setGeneric("getMeanByGroup",
           function(object) standardGeneric("getMeanByGroup")
)


#' @export
#' @aliases getTypeAnalysis
#' @rdname accessors-methods
setMethod("getTypeAnalysis", "AnalysisResults",
          function(object) return(object@type.analysis)
)

#' @export
#' @aliases getTypeData
#' @rdname accessors-methods
setMethod("getTypeData", "AnalysisResults",
          function(object) return(object@type.data)
)

#' @export
#' @aliases getDataset
#' @rdname accessors-methods
setMethod("getDataset", "AnalysisResults",
          function(object) return(object@dataset)
)


#' @rdname accessors-methods
#' @aliases getResults
#' @export
setMethod("getResults", "AnalysisResults",
          function(object) return(object@results)
)

#' @rdname accessors-methods
#' @aliases getBestModel
#' @export
setMethod("getBestModel", "AnalysisResults",
          function(object) return(object@best.model)
)

#' @export
#' @aliases getCVError
#' @rdname accessors-methods
setMethod("getCVError", "AnalysisResults",
          function(object) return(object@cv.error)
)

#' @export
#' @aliases getMeanByGroup
#' @rdname accessors-methods
setMethod("getMeanByGroup", "AnalysisResults",
          function(object) return(object@mean.by.group)
)



#### Basic methods


#' @aliases show.AnalysisResults
#' @export
#' @rdname summary-methods
setMethod(
  f = "show",
  signature = "AnalysisResults",
  definition = function(object){
    cat(paste(object@type.analysis, "performed on", object@type.data, "\n"))

    if (object@type.analysis == "PCA") {
      eigen_value <-
        data.frame(Dimension = factor(seq_len(10)),
                   "Explained variance" = object@results@modelDF$R2X * 100)

      cat("\n")
      print(eigen_value[seq_len(5), ])
      cat("[...] \n")
    } else if (object@type.analysis == "OPLS-DA") {
      cat(paste("Cross validation error:", object@cv.error))
      cat("\n\n")
      cat("Variable with the higher VIP: \n")
      print(object@mean.by.group[seq_len(10), ])
      cat("[...] \n")
    } else if (object@type.analysis == "Kruskal-Wallis tests") {
      cat("Variable with the lower adjusted p-value: \n")
      cat("\n")
      print(object@results[seq_len(10), ])
      cat("[...] \n")
    }
  }
)

#' @rdname summary-methods
#' @export
#' @aliases summary.AnalysisResults
setMethod(
  f = "summary",
  signature = "AnalysisResults",
  definition = function(object) object
)





#' Visualisation methods
#'
#' Method available to plot results of analyses in ASICS package.
#'
#' @name visualisation-methods-analyses
#' @param x An object of class \linkS4class{AnalysisResults}.
#' @param y Currently not used.
#' @param ... Currently not used.
#' @param graph A vector specifying what to plot. Allowed values are
#' \code{"eig"} for the screegraph (PCA), \code{"ind"} for plot of individuals
#' (PCA and OPLS-DA), \code{"var"} for plot of variables (PCA and OPLS-DA),
#' \code{"boxplot"} for boxplots of test results and \code{"buckets"} to show
#' significant or influential buckets on the mean spectrum.
#' Default value is \code{"default"} (\emph{i.e.,} \code{c("ind", "var")} for
#' PCA and OPLS-DA and \code{c("boxplot")} for tests).
#' @param add.label If \code{TRUE}, labels are added on individual plot.
#' @param n.label.var An integer indicating the number of label to add on
#' variable plot.
#' @param axes A numeric vector of length 2 specifying the dimensions to be
#' plotted for individual and variable plots.
#' @param col.ind A character specifying the name of the design variable used
#' to color the observations by groups for PCA individual plot.
#' @param xlim,ylim Boundaries for x and y, respectively.
#'
#'
#' @return
#' \itemize{
#' \item PCA: a \code{\link[ggplot2]{ggplot}} plot that allows for the
#' visualisation of PCA results (eigen values, individuals and variables)
#' \item OPLS-DA: a \code{\link[ggplot2]{ggplot}} plot that allows for the
#' visualisation of OPLS-DA results (individuals and variables). If
#' \code{cross.val > 1} in \code{\link{oplsda}}, the best model is plotted.
#' }
#'
#' @examples
#' # Import quantification results
#' if (require("ASICSdata", quietly = TRUE)) {
#'   quantif_path <- system.file("extdata", "results_ASICS.txt",
#'                               package = "ASICSdata")
#'   quantification <- read.table(quantif_path, header = TRUE, row.names = 1)
#'
#'   # Import design
#'   design <- read.table(system.file("extdata", "design_diabete_example.txt",
#'                                    package = "ASICSdata"), header = TRUE)
#'
#'   # Create object for analysis and remove metabolites with more than 25% of
#'   # zeros
#'   analysis_obj <- formatForAnalysis(quantification,
#'                                     zero.threshold = 25, design = design)
#'
#'   # Perform a PCA and plot results
#'   res_pca <- pca(analysis_obj)
#'   plot(res_pca)
#'
#'   # Perform an OPLS-DA and plot results
#'   res_oplsda <- oplsda(analysis_obj, "condition", orthoI = 1)
#'   plot(res_oplsda)
#' }
NULL


#' @aliases plot.AnalysisResults
#' @rdname visualisation-methods-analyses
#' @export
setMethod(
  f = "plot",
  signature = "AnalysisResults",
  definition = function(x, y, ..., graph = c("default", "ind", "var", "eig",
                                             "boxplot", "buckets"),
                        add.label = TRUE, n.label.var = 10, axes = c(1, 2),
                        col.ind = NULL, xlim = c(0.5, 10), ylim = NULL) {

      graph <- match.arg(graph)

      if (x@type.analysis == "PCA") {
        if (graph == "default") graph <- c("ind", "var")
        if(any(c("buckets", "boxplot") %in% graph))
          stop(paste("Type of plot 'buckets' or 'boxplot' are not possible",
                     "with PCA analysis"))

        .plotPCA(x@results, graph = graph, add.label = add.label,
                 axes = axes, col.ind = col.ind, nb_label = n.label.var)

      } else if (x@type.analysis == "OPLS-DA") {
        if (graph == "default") graph <- c("ind", "var")
        if(any(c("eig", "boxplot") %in% graph))
          stop(paste("Type of plot 'eig' or 'boxplot' are not possible with",
                     "OPLS-DA analysis"))
        if("buckets" %in% graph & x@type.data != "buckets")
          stop(paste("Type of plot 'buckets' is possible only with",
                     "bucket data type"))

        .plotOPLSDA(x, graph = graph, add.label = add.label,
                                xlim = xlim, ylim = ylim, nb_label = n.label.var)

      } else if (x@type.analysis == "Kruskal-Wallis tests") {
        if (graph == "default") graph <- "boxplot"
        if(any(c("eig", "var", "ind") %in% graph))
          stop(paste("Type of plot 'eig', 'var' or 'ind' are not possible with",
                     "test analysis"))
        if("buckets" %in% graph & x@type.data != "buckets")
          stop(paste("Type of plot 'buckets' is possible only with",
                     "bucket data type"))

        .plotTests(x, graph = graph, xlim = xlim, ylim = ylim)

      }
  }
)






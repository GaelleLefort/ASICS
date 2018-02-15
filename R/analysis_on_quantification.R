#' Format data for analysis
#'
#' Create an object useful for analysis with omics data and informations
#' about the omics (metabolites/genes...) and study design.
#'
#' @param data a data frame contaning omics data with samples in columns and
#' features of interest in rows (metabolites/genes...).
#' @param design a data frame describing the colums of \code{data} with at
#' least two columns and the first one corresponding to colnames of \code{data}.
#' Default to NULL (colnames of \code{data} is used for study design).
#' @param feature_info a data frame describing the rows of \code{data} with
#' at least two columns and the first one corresponding to rownames of
#' \code{data}. Default to NULL (rownames of \code{data} is used for feature
#' informations).
#'
#' @return An object of type \code{\link{SummarizedExperiment}} useful for
#' further analysis.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
format_for_analysis <- function(data, design = NULL, feature_info = NULL,
                                zero.threshold = 1, outliers = NULL){

  # if design and feature_info are NULL use rownames and colnames
  #else same order than for data
  if (is.null(design)) {
    design <- data.frame(sample.name = colnames(data))
  } else {
    if (!all(colnames(data) %in% design[, 1])){
      stop(paste("First column of 'design' matrix must contain all",
                 "colnames of 'data' matrix."))
    }
    design <- design[design[, 1] %in% colnames(data), ]
    design <- design[match(colnames(data), design[, 1]), ]
  }
  if (is.null(feature_info)) {
    feature_info <- data.frame(metabolite.name = rownames(data))
  } else {
    if (!all(rownames(data) %in% feature_info[, 1])){
      stop(paste("First column of 'feature_info' matrix must contain all",
                 "rownames of 'data' matrix."))
    }
    feature_info <- feature_info[feature_info[, 1] %in% rownames(data), ]
    feature_info <- feature_info[match(rownames(data), feature_info[, 1]), ]
  }

  raw_data <- SummarizedExperiment(assays = as.matrix(data),
                                   rowData = feature_info,
                                   colData = design)

  ## clean data
  clean_data <- raw_data
  if (zero.threshold != 1) {
    clean_data <- clean_data[apply(assay(clean_data, 1), 1,
                                   function(x) sum(x == 0)/length(x)) <
                               zero.threshold, ]
  }
  if (!is.null(outliers)) {
    clean_data <- clean_data[, !(colData(clean_data)[, 1] %in% outliers)]
  }

  return(clean_data)
}



#' @importFrom ropls opls
#' @importFrom SummarizedExperiment assay colData
#' @importFrom gridExtra grid.arrange
#' @export
#' @import ggplot2
pca <- function(analysis_data, scale.unit = TRUE){

  if (scale.unit) {
    scale_unit <- "standard"
  } else {
    scale_unit <- "none"
  }

  # PCA
  resPCA <- opls(x = t(assay(analysis_data, 1)), predI = 10,
                 scaleC = scale_unit, printL = FALSE, plotL = FALSE)

  resPCA@suppLs$yModelMN <- colData(analysis_data)

  return(resPCA)
}

#' @export
#' @importFrom gridExtra grid.arrange
plot_PCA <- function(res.pca, graph = c("ind", "var", "eig"), add.label = TRUE,
                     axes = c(1, 2), col.ind = NULL){

  eigen_value <- data.frame(dim = factor(1:10),
                            eigen_var = res.pca@pcaVarVn,
                            eigen_perc = res.pca@modelDF$R2X * 100)

  to_plot <- list()
  pos <- 1

  # graph
  if ("eig" %in% graph) {
    to_plot[[pos]] <- plot_eigen(eigen_value)
    pos <- pos + 1
  }

  if ("ind" %in% graph) {
    indiv_coord <- data.frame(res.pca@scoreMN[, axes])
    colnames(indiv_coord) <- c("x", "y")

    conditions <- NULL
    if (!is.null(col.ind)) {
      conditions <- res.pca@suppLs$yModelMN[, col.ind]
    }

    to_plot[[pos]] <- plot_indiv(indiv_coord, eigen_value, add.label, axes,
                                     conditions)
    pos <- pos + 1
  }

  if ("var" %in% graph) {
    var_coord <- data.frame(res.pca@loadingMN[, axes])
    colnames(var_coord) <- c("x", "y")

    # variable contributions on each axes
    var_contrib <- (var_coord ^ 2 * 100) / colSums(var_coord ^ 2)

    # contributions of variables on selected axes
    var_color <- (var_contrib[, 1] * eigen_value$eigen_var[axes[1]] +
                    var_contrib[, 2] * eigen_value$eigen_var[axes[2]]) /
      (eigen_value$eigen_var[axes[1]] + eigen_value$eigen_var[axes[2]])

    to_plot[[pos]] <- plot_var(var_coord, var_color, eigen_value, axes)
  }

  if (length(graph) == 1) {
    to_return <- to_plot[[1]]
    return(to_return)
  } else {
    grid.arrange(grobs = to_plot, ncol = length(graph))
  }
}


plot_eigen <- function(eigen_value){
  # label for percentage of explained variances
  text_labels <- paste0(round(eigen_value$eigen_perc, 1), "%")

  # graph
  p_eigen <- ggplot(eigen_value, aes(dim, eigen_perc, group = 1)) +
    geom_bar(stat = "identity", fill = "#377DB8", color = "#0B65B2") +
    geom_text(label = text_labels, vjust = -0.4) +
    labs(title = "PCA - Scree plot", x = "Dimensions",
         y = "Percentage of explained variances") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

  return(p_eigen)
}


plot_var <- function(var_coord, var_color, eigen_value, axes){

  # add label for the 10 most contributive variables
  coord_lab <- var_coord[var_color >= sort(var_color, decreasing = TRUE)[10], ]

  # graph limits
  limit <- max(abs(var_coord$x), abs(var_coord$y)) + 0.05


  # graph
  p_var <-
    ggplot(var_coord) +
    geom_segment(aes(x = 0, y = 0, xend = x, yend = y, color = var_color),
                 arrow = arrow(length = unit(0.10, "inches"), type = "closed"),
                 alpha = ifelse(rownames(var_coord) %in% rownames(coord_lab),
                                1, 0.3)) +
    geom_text(data = coord_lab, aes(x = x, y = y, label = rownames(coord_lab)),
              hjust = 0.5, vjust = ifelse(coord_lab$y > 0, -0.5, 1.2),
              size = 3, col = "gray30") +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2)  +
    labs(title = "PCA - Plot of variables",
         x = paste0("Dimension ", axes[1], " (",
                    round(eigen_value$eigen_perc[axes[1]], 1), "%)"),
         y = paste0("Dimension ", axes[2], " (",
                    round(eigen_value$eigen_perc[axes[2]], 1), "%)"),
         color = "Contributions") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradientn(colours = c("#1A9641", "#F7F76C", "#D7191C")) +
    coord_fixed() + scale_x_continuous(limits = c(- limit, limit)) +
    scale_y_continuous(limits = c(- limit, limit))

  return(p_var)
}


plot_indiv <- function(indiv_coord, eigen_value, add.label = TRUE,
                       axes = c(1, 2), condition = NULL){

  # labels to add
  if (add.label) {
    label <- rownames(indiv_coord)
    # label colors
    col.label <- geom_text(aes(label = label), hjust = 0.5, vjust = -0.5,
                           size = 3, col = "gray30", show.legend = F)
  } else {
    label <- ""
  }


  # if color on individuals
  ellipse <- NULL
  plot_centroids <- NULL

  if (!is.null(condition)) {
    # label colors (same color as points)
    col.label <- geom_text(aes(label = label, color = condition),
                           hjust = 0.5, vjust = -0.5, size = 3, show.legend = F)
    # add ellipses and centroids
    ellipse <-
      stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.1)
    centroids <- aggregate(cbind(x, y)~condition, data = indiv_coord, mean)
    plot_centroids <-
      geom_point(data = centroids, aes_string(x = "x", y = "y",
                                              color = "condition"), size = 4)
  }

  # graph
  p_indiv <-
    ggplot(indiv_coord, aes_string(x = "x", y = "y", color = condition)) +
    geom_point() + col.label + ellipse + plot_centroids +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(title = "PCA - Plot of individuals",
         x = paste0("Dimension ", axes[1], " (",
                    round(eigen_value$eigen_perc[axes[1]], 1), "%)"),
         y = paste0("Dimension ", axes[2], " (",
                    round(eigen_value$eigen_perc[axes[2]], 1), "%)"),
         color = "Conditions", fill = "Conditions") +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") + theme_minimal()  +
    theme(plot.title = element_text(hjust = 0.5))

  return(p_indiv)
}


#' @export
#' @importFrom ropls opls getSubsetVi predict getVipVn
#' @importFrom stats aggregate
OPLSDA <- function(analysis_data, condition, scale.unit = TRUE,
                   cross.val = 1, thres.VIP = 1){

  ## scale.unit
  if (scale.unit) {
    scale_unit <- "standard"
  } else {
    scale_unit <- "none"
  }

  ## opls-da and cross-validation
  if (cross.val > 1) {
    folds <- sample(cut(1:ncol(analysis_data), breaks = cross.val,
                        labels = FALSE), ncol(analysis_data))
    cv_oplsda_ASICS <-
      sapply(1:cross.val,
             function(i) opls(t(assay(analysis_data, 1)),
                              colData(analysis_data)[, condition], permI = 0,
                              predI = 1, orthoI = NA, scaleC = scale_unit,
                              subset = which(folds != i),
                              plotL = FALSE, printL = FALSE))
  } else {
    cv_oplsda_ASICS <-
      sapply(1:1,
             function(i) opls(t(assay(analysis_data, 1)),
                              colData(analysis_data)[, condition], permI = 0,
                              predI = 1, orthoI = NA, scaleC = scale_unit,
                              subset = "odd",
                              plotL = FALSE, printL = FALSE))
  }



  ## prediction error
  cv_error_ASICS <- round(mean(sapply(cv_oplsda_ASICS, error_pred_ASICS,
                                      analysis_data, condition)), 2)

  if (cross.val == 1) {
    cat(paste("Prediction error:", cv_error_ASICS, "\n"))
    trainVi <- getSubsetVi(cv_oplsda_ASICS[[1]])
    comparison <- table(colData(analysis_data)[, condition][-trainVi],
                        predict(cv_oplsda_ASICS[[1]],
                                t(assay(analysis_data, 1))[-trainVi, ]))
    cat("Confusion matrix:")
    print(comparison)
    cat("\n\n")
  } else {
    cat(paste("Cross-validation prediction error:", cv_error_ASICS, "\n"))
    cat("\n\n")
  }


  ## Influencial metabolites
  VIP_all <- sapply(cv_oplsda_ASICS, getVipVn)
  mean_VIP <- apply(VIP_all, 1, mean)
  VIP_all <- as.data.frame(VIP_all[order(mean_VIP, decreasing = TRUE), ])

  stab_metab <- rownames(VIP_all)[apply(VIP_all, 1, function(x) sum(x >= thres.VIP)) >= cross.val - 1]

  mean_by_group <- aggregate(t(assay(analysis_data, 1)),
                             list(colData(analysis_data)[, condition]), mean)
  rownames(mean_by_group) <- mean_by_group[, 1]
  mean_by_group[, 1] <- NULL
  mean_by_group <- t(mean_by_group)
  mean_by_group <- mean_by_group[order(mean_VIP, decreasing = TRUE),]
  mean_by_group_influencial <- mean_by_group[rownames(mean_by_group) %in%
                                               stab_metab, ]


  cat("Most influential metabolites (mean of concentration in each group): \n")
  print(head(mean_by_group_influencial))

  return(list(cv_oplsda = cv_oplsda_ASICS, cv_error = cv_error_ASICS,
              mean_influencial = mean_by_group_influencial))

}

error_pred_ASICS <- function(x, analysis_data, condition){
  trainVi <- getSubsetVi(x)
  comparison <- table(colData(analysis_data)[, condition][-trainVi],
                      predict(x, t(assay(analysis_data, 1))[-trainVi, ]))
  return(sum(comparison[row(comparison) != col(comparison)]) / sum(comparison))
}


#' @export
plot_OPLSDA <- function(res.oplsda, graph = c("ind", "var"), add.label = TRUE){

  to_plot <- list()
  pos <- 1

  # eigen values
  eigen_value <- data.frame(dim = factor(1:2),
                            eigen_var = c(1, 1),
                            eigen_perc = res.oplsda$cv_oplsda[[1]]@modelDF$R2X[1:2] * 100)


  if ("ind" %in% graph) {
    indiv_coord <- data.frame(x = res.oplsda$cv_oplsda[[1]]@scoreMN,
                              y = res.oplsda$cv_oplsda[[1]]@orthoScoreMN,
                              row.names = rownames(res.oplsda$cv_oplsda[[1]]@scoreMN))
    colnames(indiv_coord) <- c("x", "y")



    to_plot[[pos]] <-
      plot_indiv(indiv_coord, eigen_value, add.label = add.label, axes = c(1, 2),
                 condition = res.oplsda$cv_oplsda[[1]]@suppLs$y[getSubsetVi(res.oplsda$cv_oplsda[[1]])]) +
      coord_cartesian() + ylab("Orthogonal component 1") +
      xlab("Dimension 1") +
      ggtitle("OPLS-DA - Plot of individuals")

    pos <- pos + 1
  }

  if ("var" %in% graph) {
    var_coord <- data.frame(x = res.oplsda$cv_oplsda[[1]]@loadingMN,
                              y = res.oplsda$cv_oplsda[[1]]@orthoLoadingMN,
                              row.names = rownames(res.oplsda$cv_oplsda[[1]]@loadingMN))
    colnames(var_coord) <- c("x", "y")

    # contributions of variables on selected axes
    var_color <- apply(sapply(res.oplsda$cv_oplsda, getVipVn), 1, mean)

    to_plot[[pos]] <- plot_var(var_coord, var_color, eigen_value, axes = c(1, 2)) +
      labs(title = "OPLS-DA - Plot of variables",
           y = "Orthogonal component 1",
           x = "Dimension 1",
           color = "VIP")
  }

  if (length(graph) == 1) {
    to_return <- to_plot[[1]]
    return(to_return)
  } else {
    grid.arrange(grobs = to_plot, ncol = length(graph))
  }
}




#' @importFrom stats wilcox.test p.adjust
#' @export
kruskal_on_quantification <- function(analysis_data, condition,
                                       correction = "BH", alpha = 0.05){

  kruskal_pval <-
    apply(t(assay(analysis_data, 1)), 2,
          function(x) kruskal.test(x~colData(analysis_data)[, condition])$p.value)

  p_adj <- p.adjust(kruskal_pval, method = correction)

  mean_by_group <- aggregate(t(assay(analysis_data, 1)),
                             list(colData(analysis_data)[, condition]), mean)
  rownames(mean_by_group) <- mean_by_group[, 1]
  mean_by_group[, 1] <- NULL
  mean_by_group <- t(mean_by_group)
  mean_by_group_sig <- mean_by_group[rownames(mean_by_group) %in%
                                               names(p_adj[p_adj < alpha]), ]
  colnames(mean_by_group_sig) <- paste("mean_", colnames(mean_by_group_sig))

  res_tests <- cbind(p_value_adj = p_adj[p_adj < alpha], mean_by_group_sig)

  return(res_tests)
}




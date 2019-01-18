library(ASICS)

context("Analysis on ASICS relative quantifications")

check_ASICSdata <- function() {
  if (!require("ASICSdata")) {
    skip("ASICSdata package not available")
  }
}


test_that("PCA on quantification goes well", {
  check_ASICSdata()

  # import quantification results
  quantif_path <- system.file("extdata", "results_ASICS.txt",
                              package = "ASICSdata")
  quantification <- read.table(quantif_path, header = TRUE, row.names = 1)

  analysis_obj <- formatForAnalysis(quantification,
                                    zero.threshold = 25)

  expect_error(pca(analysis_obj), NA)
  expect_error(plot(pca(analysis_obj)), NA)
  expect_error(plot(pca(analysis_obj), graph = "buckets"),
               paste("Type of plot 'buckets' or 'boxplot' are not possible",
                     "with PCA analysis"))
})


test_that("OPLS-DA on quantifications works well", {
  check_ASICSdata()

  # import quantification results
  quantif_path <- system.file("extdata", "results_ASICS.txt",
                              package = "ASICSdata")
  quantification <- read.table(quantif_path, header = TRUE, row.names = 1)

  # import design
  design <- read.table(system.file("extdata", "design_diabete_example.txt",
                                   package = "ASICSdata"), header = TRUE)

  analysis_obj <- formatForAnalysis(quantification,
                                    zero.threshold = 25, design = design)

  res_oplsda <- oplsda(analysis_obj, "condition", orthoI = 1)
  expect_error(plot(res_oplsda), NA)

  spectra_path <- system.file("extdata", "spectra_diabetes_example.txt",
                              package = "ASICSdata")
  spectra <- read.table(spectra_path, header = TRUE, row.names = 1)
  analysis_obj_buck <- formatForAnalysis(binning(spectra),
                                         zero.threshold = 25, design = design)
  expect_error(plot(oplsda(analysis_obj_buck, "condition", orthoI = 1,
                           type.data = "buckets"),
                    graph = "buckets"), NA)

})


test_that("Kruskall-Wallis on quantifications works well", {
  check_ASICSdata()

  # import quantification results
  quantif_path <- system.file("extdata", "results_ASICS.txt",
                              package = "ASICSdata")
  quantification <- read.table(quantif_path, header = TRUE, row.names = 1)

  # import design
  design <- read.table(system.file("extdata", "design_diabete_example.txt",
                                   package = "ASICSdata"), header = TRUE)

  # create object for analysis and remove metabolites with more than 25% of
  #zeros
  analysis_obj <- formatForAnalysis(quantification,
                                    zero.threshold = 25, design = design)
  expect_error(kruskalWallis(analysis_obj, "condition"), NA)
  expect_error(plot(kruskalWallis(analysis_obj, "condition")), NA)
})

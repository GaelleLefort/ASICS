library(ASICS)

context("Quantification of NMR spectra")

test_that("Quantification with ASICS works well and gives a well formatted
S4 object", {
  current_path <- file.path(system.file("extdata", package = "ASICS"),
                            "spectra_example.txt")
  spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
  spectra_obj <- createSpectra(spectra_data)

  to_exlude <- matrix(c(4.5, 5.1, 5.5, 6.5), ncol = 2)

  expect_s4_class(ASICS(spectra_obj, exclusion.areas = to_exlude),
                  "ASICSResults")
})

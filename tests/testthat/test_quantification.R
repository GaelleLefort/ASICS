library(ASICS)

context("Quantification of NMR spectra")

test_that("Quantification with ASICS works well and gives a well formatted
S4 object (first version of algorithm)", {
  current_path <- system.file("extdata", package = "ASICS")
  spectra_data <- importSpectra(name.dir = current_path,
                                name.file = "spectra_example.txt",
                                type.import = "txt")
  spectra_obj <- createSpectra(spectra_data)

  to_exclude <- matrix(c(4.5, 10), ncol = 2)

  expect_s4_class(ASICS(spectra_obj, exclusion.areas = to_exclude, combine = FALSE),
                  "ASICSResults")
})


test_that("Quantification with ASICS works well and gives a well formatted
S4 object (second version of algorithm)", {
  current_path <- system.file("extdata", "example_spectra", package = "ASICS")
  spectra_data <- importSpectra(name.dir = current_path, type.import = "1r")
  spectra_obj <- createSpectra(spectra_data)
  expect_s4_class(ASICS(spectra_obj, combine = TRUE),
                  "ASICSResults")
})

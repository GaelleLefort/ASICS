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
  pure_lib <- pure_library[getSampleName(pure_library) %in%
                             c("Lactate", "L-Alanine")]

  expect_s4_class(ASICS(spectra_obj[1], exclusion.areas = to_exclude,
                        pure.library = pure_lib),
                  "ASICSResults")
})


test_that("Quantification with ASICS works well and gives a well formatted
S4 object (second version of algorithm)", {
  current_path <- system.file("extdata", package = "ASICS")
  spectra_data <- importSpectra(name.dir = current_path,
                                name.file = "spectra_example.txt",
                                type.import = "txt")
  spectra_obj <- createSpectra(spectra_data)

  to_exclude <- matrix(c(4.5, 10), ncol = 2)
  pure_lib <- pure_library[getSampleName(pure_library) %in%
                             c("Lactate", "L-Alanine")]

  expect_s4_class(ASICS(spectra_obj[1], exclusion.areas = to_exclude,
                        pure.library = pure_lib, combine = TRUE),
                  "ASICSResults")
})

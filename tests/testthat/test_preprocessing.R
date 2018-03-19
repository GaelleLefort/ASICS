library(ASICS)

context("Preprocessing the spectra before quantification")

test_that("Importation from Bruker files works well and gives a well formatted
data frame", {
  lib_path <- system.file("extdata", "example_library",
                          package = "ASICS")
  expect_is(importSpectraBruker(lib_path), "data.frame")

  example_path <- system.file("extdata", "example_spectra",
                              package = "ASICS")
  expect_warning(importSpectraBruker(example_path),
                 regexp = "There is a problem in files for spectra:",
                 fixed = TRUE)

  expect_is(as.numeric(rownames(importSpectraBruker(lib_path))),
            "numeric")
})


test_that("Preprocessing functions work well and give a data frame", {
  current_path <- system.file("extdata", "spectra_example.txt",
                              package = "ASICS")
  spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
  expect_is(normalisation(spectra_data), "data.frame")
  expect_is(baselineCorrection(spectra_data), "data.frame")
  expect_is(alignment(spectra_data), "data.frame")
})

test_that("Bucketing function works well and gives a data frame", {
  current_path <- system.file("extdata", "spectra_example.txt",
                              package = "ASICS")
  spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
  expect_is(binning(spectra_data), "data.frame")
})


test_that("Functions of PureLibrary class work without error", {
  expect_equal(length(pure_library), length(pure_library@sample.name))
  expect_equal(length(pure_library[1]), 1)
  expect_equal(length(c(pure_library[1], pure_library[2])), 2)
})


test_that("Creation of a 'Spectra' objet and functions of Spectra class work
without error", {
  current_path <- system.file("extdata", "spectra_example.txt",
                              package = "ASICS")
  spectra_data <- read.table(current_path, header = TRUE, row.names = 1)
  spectra_data <- baselineCorrection(spectra_data)
  spectra_obj <- createSpectra(spectra_data)

  expect_s4_class(spectra_obj, "Spectra")
  expect_equal(length(spectra_obj[1]), 1)
  expect_equal(length(c(spectra_obj[1], spectra_obj[2])), 2)

})


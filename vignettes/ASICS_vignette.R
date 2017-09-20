## ----lib-----------------------------------------------------------------
data(pure_library)
pure_library$name

## ----ASICS_onefile-------------------------------------------------------
result <- ASICS(path = system.file("extdata", "example_spectra",
                                   "AG_faq_Beck01", package = "ASICS"),
                exclusion.areas = matrix(c(4.5,5.1,5.5,6.5),
                                         ncol = 2, byrow = TRUE),
                max.shift = 0.02, which.spectra = "last",
                library.metabolites = NULL, threshold.noise = 0.02)


## ----plot_spectrum, warning=FALSE, fig.width=10, fig.height=6------------
plot_spectrum(result, xmin = 1, xmax = 1.5, ymax = 10, add_metab = "Lactate")

## ----rel_conc------------------------------------------------------------
result$Present

## ----ASICS_multifile, results='hide'-------------------------------------
res_multi <- ASICS_multiFiles(name.dir = system.file("extdata",
                                                     "example_spectra",
                                                     package = "ASICS"),
                              exclusion.areas = matrix(c(4.5,5.1,5.5,6.5),
                                                       ncol = 2, byrow = TRUE),
                              max.shift = 0.02, which.spectra = "last",
                              library.metabolites = NULL,
                              threshold.noise = 0.02, ncores = 4)

## ----extract_concentrations----------------------------------------------
quantification <- extract_concentrations(res_multi)
head(quantification)


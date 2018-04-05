## ----ASICSload-----------------------------------------------------------
library(ASICS)
library(ASICSdata)

## ----lib-----------------------------------------------------------------
head(getSampleName(pure_library), n = 8)

## ----import_create_lib, results='hide'-----------------------------------
pure_spectra <- importSpectraBruker(system.file("extdata", "example_library", 
                                                package = "ASICS"))
new_pure_library <- createPureLibrary(pure_spectra, 
                                        nb.protons = c(5, 4))

## ----select_merge_lib----------------------------------------------------
merged_pure_library <- c(pure_library[1:10], new_pure_library)

## ----import_from_Bruker, results='hide'----------------------------------
spectra_data <- importSpectraBruker(system.file("extdata", 
                                                "Human_diabetes_example", 
                                                package = "ASICSdata"))

## ----import_from_txt, results='hide'-------------------------------------
diabetes <- system.file("extdata", "spectra_diabetes_example.txt", 
                        package = "ASICSdata")
spectra_data_txt <- read.table(diabetes, header = TRUE, row.names = 1)

## ----preprocessing, results='hide'---------------------------------------
spectra_base_cor <- baselineCorrection(spectra_data_txt)

## ----create_spectra, results='hide'--------------------------------------
spectra_obj <- createSpectra(spectra_base_cor)

## ----ASICS, results='hide'-----------------------------------------------
# part of the spectrum to exclude (water and urea)
to_exclude <- matrix(c(4.5, 5.1, 5.5, 6.5), ncol = 2, byrow = TRUE)
ASICS_results <- ASICS(spectra_obj, exclusion.areas = to_exclude)

## ----summary_res---------------------------------------------------------
ASICS_results

## ----plot_spectrum, warning=FALSE, fig.width=12, fig.height=8------------
plot(ASICS_results, idx = 1, xlim = c(2.8, 3.3), add.metab = "Creatinine")

## ----rel_conc------------------------------------------------------------
head(getQuantification(ASICS_results), 10)[, 1:2]

## ----design--------------------------------------------------------------
design <- read.table(system.file("extdata", "design_diabete_example.txt", 
                                 package = "ASICSdata"), header = TRUE)

## ----analyses_obj--------------------------------------------------------
analysis_data <- formatForAnalysis(getQuantification(ASICS_results),
                                   design = design, zero.threshold = 75,
                                   zero.group = "condition", 
                                   outliers = c("ADG10003u_007", 
                                                "ADG19007u_163"))

## ----pca, fig.width=10---------------------------------------------------
resPCA <- pca(analysis_data)
plot(resPCA, graph = "ind", col.ind = "condition")
plot(resPCA, graph = "var")

## ----oplsda--------------------------------------------------------------
resOPLSDA <- oplsda(analysis_data, condition = "condition", orthoI = 1)
resOPLSDA

## ----oplsda_plot, fig.width=10, results='hide'---------------------------
plot(resOPLSDA)

## ----perform_tests-------------------------------------------------------
resTests <- kruskalWallis(analysis_data, "condition")
resTests

## ----test_plot, fig.width=10, fig.height=10, results='hide'--------------
plot(resTests)

## ----align_and_binning, results='hide'-----------------------------------
spectra_align <- alignment(spectra_data_txt)
spectra_bucket <- binning(spectra_align)

## ----create_ana_obj------------------------------------------------------
analysis_data_bucket <- formatForAnalysis(spectra_bucket, design = design,
                                          zero.threshold = 75)

## ----oplsda_buckets------------------------------------------------------
resOPLSDA_buckets <- oplsda(analysis_data_bucket, condition = "condition",
                            type.data = "buckets")
resOPLSDA_buckets

## ----oplsda_buckets_plot, fig.width=10-----------------------------------
plot(resOPLSDA_buckets, graph = "buckets")

## ----sysinfo-------------------------------------------------------------
sessionInfo()


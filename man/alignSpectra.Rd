% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data.R
\name{alignSpectra}
\alias{alignSpectra}
\title{Alignment}
\usage{
alignSpectra(
  spectra,
  reference = NULL,
  max.shift = 0.02,
  ncores = 1,
  verbose = TRUE
)
}
\arguments{
\item{spectra}{Data frame with spectra in columns and chemical shift in rows.
Colnames of this data frame correspond to pure metabolite names and rownames
to chemical shift grid (in ppm).}

\item{reference}{Index of the reference spectrum used for the alignment.
Default to \code{NULL}, \emph{i.e.} the reference spectrum is automatically
detected.}

\item{max.shift}{Maximum shift allowed for the alignment. Default to 0.002.}

\item{ncores}{Number of cores used in parallel evaluation. Default to
\code{1}.}

\item{verbose}{A boolean value to allow print out process information.}
}
\value{
A data frame with aligned spectra in columns and chemical shifts
(in ppm) in rows.
}
\description{
Align spectra of a data frame by a method based on the CluPA algorithm
(Vu et al., (2011))
}
\examples{
current_path <- system.file("extdata", package = "ASICS")
spectra_data <- importSpectra(name.dir = current_path,
                     name.file = "spectra_example.txt", type.import = "txt")
spectra_align <- alignSpectra(spectra_data)
}
\references{
Vu, T. N., Valkenborg, D., Smets, K., Verwaest, K. A., Dommisse,
R., Lemiere, F., ... & Laukens, K. (2011). An integrated workflow for robust
alignment and simplified quantitative analysis of NMR spectrometry data.
\emph{BMC Bioinformatics}, \strong{12}(1), 405.
}

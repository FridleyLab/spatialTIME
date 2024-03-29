% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_immunoflo.R
\name{plot_immunoflo}
\alias{plot_immunoflo}
\title{Generate plot of TMA point process}
\usage{
plot_immunoflo(
  mif,
  plot_title,
  mnames,
  mcolors = NULL,
  cell_type = NULL,
  filename = NULL,
  path = NULL,
  xloc = NULL,
  yloc = NULL
)
}
\arguments{
\item{mif}{MIF object created using create_MIF().}

\item{plot_title}{Character string or vector of character strings of variable name(s) to serve as plot title(s).}

\item{mnames}{Character vector containing marker names.}

\item{mcolors}{Character vector of color names to display markers in the plot.}

\item{cell_type}{Character vector of cell type}

\item{filename}{Character string of file name to store plots. Plots are generated as single .pdf file.}

\item{path}{Different path than file name or to use in conjunction with filename ???}

\item{xloc, yloc}{columns in the spatial files containing the x and y locations of cells. Default is `NULL` which will result in `xloc` and `yloc` being calculated from `XMin`/`YMin` and `XMax`/`YMax`}
}
\value{
mif object and the ggplot objects can be viewed form the derived slot of the mif object
}
\description{
This function generates plot of point process in rectangular or circular window.
}
\examples{
#Create mif object
library(dplyr)
x <- create_mif(clinical_data = example_clinical \%>\% 
mutate(deidentified_id = as.character(deidentified_id)),
sample_data = example_summary \%>\% 
mutate(deidentified_id = as.character(deidentified_id)),
spatial_list = example_spatial,
patient_id = "deidentified_id", 
sample_id = "deidentified_sample")

mnames_good <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
"FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
"PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")

x <- plot_immunoflo(x, plot_title = "deidentified_sample", mnames = mnames_good, 
cell_type = "Classifier.Label")

x[["derived"]][["spatial_plots"]][[4]]
}

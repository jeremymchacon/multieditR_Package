% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect_edits_batch.R
\name{detect_edits_batch}
\alias{detect_edits_batch}
\title{run edit detection over multiple samples}
\usage{
detect_edits_batch(params = NULL)
}
\arguments{
\item{params}{A data.frame with the following columns: sample_file, control_file,
motif, motif_fwd, wt, edit, phred_cutoff, p_value. see multiEditR::detect_edits for details.}
}
\value{
A list of multiEditR objects
}
\description{
run edit detection over multiple samples
}

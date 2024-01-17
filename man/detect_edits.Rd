\name{detect_edits}
\alias{detect_edits}
\title{detect_edits}
\usage{
detect_edits("path/to/sample.ab1","path/to/control.ab1", motif = "AGTAGCTGGGATTACAGATG", wt = "A", edit = "G")
}
\description{
Compares a sanger sequence to a control sanger sequence or ab1 to predict edits
}
\examples{
sample_file = system.file("extdata", "RP272_cdna_wt.ab1", package="multieditR")
ctrl_file = system.file("extdata", "RP272_cdna_ko.ab1", package="multieditR")
motif = "AGTAGCTGGGATTACAGATG"
fit = detect_edits(sample_file, ctrl_file,
            motif = motif, wt = "A", edit = "G")
# chromatogram of predicted edits in sample
geom_chromatogram(fit$sample_sanger,fit$sample_locs[1], fit$sample_locs[2])

# chromatogram of predicts "edits" in control
geom_chromatogram(fit$ctrl_sanger,fit$ctrl_locs[1], fit$ctrl_locs[2])

# percent signal of basecalls in sample
plot_raw_sample(fit)

# percent noise of basecalls 
plot_trimmed_sample(fit)
}
\arguments{
  \item{sample_file}{path to the sample .ab1 file}
  \item{ctrl_file}{path to the control .ab1 file or .fasta file}
  \item{motif}{target sequence to examine for edits}
  \item{wt}{base targeted for editing}
  \item{edit}{base the wt base should be altered to}

}

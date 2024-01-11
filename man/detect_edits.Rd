\name{detect_edits}
\alias{detect_edits}
\title{detect_edits}
\usage{
detect_edits(sample_name, "path/to/sample.ab1","path/to/control.ab1")
}
\description{
Compares a sanger sequence to a control sanger sequence or ab1 to predict edits
}
\examples{
sample_file = system.file("extdata", "RP272_cdna_wt.ab1", package="multieditR")
ctrl_file = system.file("extdata", "RP272_cdna_ko.ab1", package="multieditR")
motif = "AGTAGCTGGGATTACAGATG"
detect_edits("example", sample_file, ctrl_file,
            motif = motif, wt = "A", edit = "G")
}
\arguments{
  \item{sample_name}{string describing the sample name}
  \item{sample_file}{path to the sample .ab1 file}
  \item{ctrl_file}{path to the control .ab1 file or .fasta file}

}

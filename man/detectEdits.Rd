\name{detectEdits}
\alias{detectEdits}
\title{detectEdits}
\usage{
detectEdits(sample_name, "path/to/sample.ab1","path/to/control.ab1")
}
\description{
Compares a sanger sequence to a control sanger sequence or ab1 to predict edits
}
\examples{
detectEdits(sample_name, "path/to/sample.ab1","path/to/control.ab1")
}
\arguments{
  \item{sample_name}{string describing the sample name}
  \item{sample_file}{path to the sample .ab1 file}
  \item{ctrl_file}{path to the control .ab1 file or .fasta file}

}

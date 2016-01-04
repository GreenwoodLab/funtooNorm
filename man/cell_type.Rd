\name{cell_type}
\alias{cell_type}
\docType{data}
\title{
Cell types
}
\description{
Vector indicating the cell type information for the samples to be normalized by funtooNorm.  The length of the vector should correspond to the number of samples to be normalized, and the order should correspond to the order of the columns in the signal A and signal B matrices.
}
\usage{cell_type}

\format{
  A vector with 93 rows, one per sample, indicating the cell type
  information for each sample. Coded as factor. The sample names are
  used to name the vector elements.
}

\name{funtoonorm}
\alias{funtoonorm}

\title{
A function to normalize Illumina Infinium Human Methylation 450 BeadChip (Illumina 450K) with multiple tissues or cell types.
}
\description{
This function performs normalization of Illumina Infinium Human Methylation 450 BeadChip  data using the information contained in the control probes.  It implements different corrections at different quantiles, and allows for the normalization corrections to vary across tissues/cell types. 
}
\usage{
funtoonorm(sigA, sigB, Annot = NULL, controlred, controlgrn,
cp.types = NULL, cell_type, ncmp = 4, save.quant = TRUE, save.loess = TRUE,
apply.loess = TRUE,  validate = FALSE)
}

\arguments{
\item{sigA, sigB}{
Matrices containing the signal A and signal B results extracted from the IDAT files.
}
\item{controlred, controlgrn}{
Matrices containing control probe data. 
}
\item{Annot}{
Annotation matrix.  Supplied by default.
}
 
\item{cp.types}{
Vector of types of control probes. 
}
\item{cell_type}{
Vector of cell type (or tissue type) information.
}
\item{ncmp}{
Number of partial least squares components used in the model fitting.
}
\item{save.quant}{
Logical, whether to save calculated quantiles.
}
\item{save.loess}{
Logical, whether to save calculated results of loess regression.
}
\item{apply.loess}{
Logical, whether to apply results of loess regression. If TRUE, two matrices are returned, one the data before normalization and one after normalization.  
normalised beta values is returned.
}
\item{validate}{
Either FALSE, or the maximum number of PLS components to be explored in cross-validation.
If FALSE, the normalization corrections are calculated using \verb{ncmp} partial least squares (PLS) components.  if not FALSE, then a number must be supplied.  This number will be the maximum number of PLS components used when exploring model fit performance across a range of \verb{ncmp} values ranging from 1 to the supplied number.
}
}
\details{
The funtooNorm function operates in one of two modes. If validate=FALSE, then the normalization corrections are calculated using the supplied value of \verb{ncmp} to fix the number of partial least squares (PLS) components.  If validate is an integer, K>1, (e.g. K=5), then cross-validation is performed exploring performance across a range of values for \verb{ncmp} ranging from 1 to K.
}
\value{The values returned depend on the parameters chosen.
\itemize{
\item If validate is not FALSE (i.e. validate=K), the function creates a pdf file containing a series of plots showing residual error curves across percentiles of the signal distributions, to facilitate the choice of an appropriate value for \verb{ncmp}.   No object is returned by the function.

\item If validate = FALSE, then funtoonorm has the following behaviour:
  \itemize{
   \item If apply.loess = FALSE the function will not return any object.  
         However, if save.loess=TRUE or if save.quant=TRUE then RData objects will be saved to disk for future use.
   \item If apply.less= TRUE, then the function returns a list of 2 objects.  The first, 'origBeta', is the matrix of Beta avalues before normalization, and the second, 'newBeta' is the Beta values after normalization.
  }
}
}

\examples{
%% to normalize methylation data:
ncmp <- 4
funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample, Annot=Annotsample, 
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type, 
                       ncmp=ncmp, save.quant=TRUE, save.loess=TRUE, apply.loess=TRUE, 
                       validate=FALSE)
%%to choose the number of components: 
funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample,
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type,
                      ncmp=4, save.quant=TRUE, save.loess=TRUE, 
                      apply.loess=FALSE, validate=5)
}

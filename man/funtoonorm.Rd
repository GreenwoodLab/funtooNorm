\name{funtoonorm}
\alias{funtoonorm}

\title{
A function to normalize Illumina Infinium Human Methylation 450 BeadChip (Illumina 450K) with multiple tissues or cell types.
}
\description{
Function performs normalization Illumina Infinium Human Methylation 450 BeadChip  data, 
correcting for quantile and tissue/cell type. 
}
\usage{
funtoonorm(sigA, sigB, Annot = NULL, controlred, controlgrn,
cp.types = NULL, cell_type, ncmp = 4, save.quant = TRUE, save.loess = TRUE,
apply.loess = TRUE, logit.quant = TRUE, validate = FALSE)
}

\arguments{
\item{sigA, sigB}{
Matrices containing the signal A and signal B results extracted from the IDAT files.
}
\item{controlred, controlgrn}{
Matrices containing control probe data. 
}
\item{Annot}{
Annotation matrix.
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
Logical, whether to apply results of loess regression. If TRUE, matrix of 
normalised beta values is returned.
}
\item{logit.quant}{
Whether to apply logit transformation to quantiles.
}
\item{validate}{
Either FALSE, or the maximum number of PLS components to be explored in cross-validation.
If FALSE, the normalization corrections are calculated using \emph{ncmp} partial least squares (PLS) components.  if not FALSE, then a number must be supplied.  This number will be the maximum number of PLS components used when exploring model fit performance across a range of \emph{ncmp} values ranging from 1 to the supplied number.
}
}
\details{
The funtooNorm function operates in one of two modes. If validate=FALSE, then the normalization corrections are calculated using the supplied value of \emph{ncmp} to fix the number of partial least squares (PLS) components.  If validate is an integer, $K>1$, (e.g.\ $K=5$), then cross-validation is performed exploring performance across a range of values for \emph{ncmp} ranging from 1 to $K$.
}
\value{
If validate is not FALSE (i.e. validate=$K$), the function creates a pdf file containing a series of plots showing residual error curves across percentiles of the signal distributions, to facilitate the choice of an appropriate value for \emph{ncmp}. 

If validate = FALSE, then the function may return either completely normalized data (if apply.loess=TRUE), or save some intermediate calulations (if apply.loess=FALSE but save.loess=TRUE or save.quant=TRUE).
}

\examples{
%% to normalize methylation data:
ncmp <- 4
funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample, Annot=Annotsample, 
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type, 
                       ncmp=ncmp, save.quant=TRUE, save.loess=TRUE, apply.loess=TRUE, 
                       logit.quant=TRUE, validate=FALSE)
%%to validate: 
funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample,
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type,
                      ncmp=4, save.quant=TRUE, save.loess=TRUE, 
                      apply.loess=FALSE, logit.quant=TRUE, validate=5)
}

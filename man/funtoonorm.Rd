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
cp.types = NULL, cell_type, ncmp = 4, 
ncv.fold = 10, logged2.data=FALSE, save.quant=TRUE,
type.fits="PCR", apply.fits = TRUE, validate = FALSE)
}

\arguments{
\item{sigA, sigB}{
Matrices containing the signal A and signal B results extracted from the IDAT files.  (The user should know whether or not these data have undergone a log transformation prior to being sent into funtooNorm.)
}
\item{Annot}{
Annotation matrix.  Supplied by default.
}
\item{controlred, controlgrn}{
Matrices containing control probe data. (Again, the user should know whether the data have been log transformed prior to calling funtooNorm).
}
\item{cp.types}{
Vector of types of control probes. (Often available as the row names from the control probe information).
}
\item{cell_type}{
Vector of cell type (or tissue type) information.
}
\item{ncmp}{
Number of components, from either principal component regression or partial least squares regression, used in the model fitting when \verb{validate=FALSE}.
}
\item{ncv.fold}{
Number of cross-validation partitions.
}
\item{logged2.data}{
Logical, \verb{TRUE} if data have been previously log2 transformed (sigA, sigB, controlred, controlgrn), and \verb{FALSE} if not.
}

\item{save.quant}{
Logical, whether to save calculated quantiles.
}
\item{type.fits}{
Values are either \verb{"PLS"}, or \verb{"PCR"} for partial least squares or principal component regression, respectively. The default is \verb{"PCR"}.
}
\item{apply.fits}{
Logical, whether to apply results of the PCR/PLS fits to normalize the data. If \verb{TRUE}, two matrices are returned, one the data before normalization and one after normalization.  
This would normally be set to \verb{FALSE} when \verb{validate=TRUE}. 
}
\item{validate}{
Logical.  If \verb{FALSE}, the normalization corrections are calculated using \verb{ncmp} components.  
If \verb{TRUE}, then a range of numbers of components are explored graphically. 
}
}
\details{
The funtooNorm function operates in one of two modes. If \verb{validate=FALSE}, then the normalization corrections are calculated using the supplied value of \verb{ncmp} to specify the number of components.  
If \verb{validate=TRUE}, then cross-validation is performed exploring performance across a range of values.  The maximum value is set to 8 to facilite the visualization of results.
}
\value{The values returned depend on the parameters chosen.
\itemize{
\item If validate is TRUE the function creates two pdf files, each containing a series of plots showing root mean squared error curves across percentiles of the signal distributions, obtained from cross-validation, to facilitate the choice of an appropriate value for \verb{ncmp}.   No object is returned by the function.

\item If validate is FALSE, then \verb{funtoonorm} has the following behaviour:
  \itemize{
   \item If apply.fits = FALSE the function will not return any object.  
         However, if save.quant=TRUE then RData objects will be saved to disk for future use.
   \item If apply.fits= TRUE, then the function returns a list of 2 objects.  The first, 'origBeta', is the matrix of Beta values before normalization, and the second, 'newBeta' is the Beta values after normalization.
  }
}
}

\seealso{\link[minfi]{preprocessFunnorm} in the \pkg{minfi} package}

\examples{
## to choose the number of components:
funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample,
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=NULL, cell_type = cell_type,
                      logged2.data=FALSE, save.quant=TRUE, apply.fits=FALSE,  
                      validate=TRUE)
## to normalize methylation data, assuming save.quant was set to TRUE in the previous call:
ncmp <- 4
funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample, Annot=Annot,
                      controlred=matred, controlgrn=matgrn,
                      cp.types=NULL, cell_type = cell_type,
                      logged2.data=FALSE, save.quant=FALSE, ncmp=ncmp, apply.fits=TRUE,
                      validate=FALSE)

}

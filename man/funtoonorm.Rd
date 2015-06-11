\name{funtoonorm}
\alias{funtoonorm}

\title{
Function to normalize Illumina Infinium Human Methylation 450 BeadChip (Illumina 450K)
}
\description{
Function performs normalization Illumina Infinium Human Methylation 450 BeadChip  data, 
correcting for quantile and tissue/cell type. 
}
\usage{
funtoonorm(sigA, sigB, Annot = NULL, quantiledat = NULL, controlred, controlgrn,
cp.types = NULL, cell_type, ncmp = 4, save.quant = TRUE, save.loess = TRUE,
apply.loess = TRUE, logit.quant = TRUE, validate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{sigA, sigB}{
Matrices with signals A and B.
}
\item{controlred, controlgrn}{
Matrices with control probes signal.
}

\item{Annot}{
Annotation matrix.
}
 
\item{cp.types}{
Vector of types of control probes.
}
\item{cell_type}{
Cell types vector.
}
\item{ncmp}{
Number of PLS components.
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
If not FALSE, takes the maximum number of PLS components for validation.
}
}
\details{
Function works in two regimes, validation and  normalization. Validation regime 
is triggered by setting validate variable to the integer value bigger that 1.      
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If the validation regime is triggered (see details) the function creates a pdf 
file with validation curves. In normalization regime output can be either results 
of loess fitting (Celia, could you add better definitions here) or normalised 
beta values, depending on the value of apply.loess.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
%% to normalize methylation data:
data("data")
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
                      apply.loess=FALSE, logit.quant=FALSE, validate=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

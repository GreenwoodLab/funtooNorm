\name{agreement}
 \alias{agreement}
 
 \title{
 Function to measure intra-replicate agreement in methylation data.
 }
 \description{
 Function calculates the sum of squared differences between sets of replicates, using
 methylation data, in order to evaluate the performance of the funtooNorm function.
 }
 \usage{
 agreement(Beta, individualID)
 }
 %- maybe also 'usage' for other objects documented here.
 \arguments{
   \item{Beta}{
 Matrix with beta-values, rows corresponding to probes, columns corresponding to samples. 
 }
  \item{individualID}{
 Vector of individual ID's matching the order of the samples (columns) in Beta. 
 Only individualID's that occur more than once (i.e. replicates) are used in calculating agreement.
 }
 }
 \details{
 We expect that the values returned by the agreement function after normalization by funtooNorm to be smaller than before.
 }
 \value{
 Function returns a measure of agreement between replicates in
 methylation data, for details, see our paper.
 }
 \references{
 Klein K.O. et al., (2015). \bold{funtooNorm:} An improvement of the funNorm
 normalization method for methylation data from multiple cell or tissue types.
 Manuscript in preparation.
 }
 
 \examples{

     #This call  will perform normalization of the data:
     ncmp <- 4
     funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample, Annot=Annot,
                      controlred=matred, controlgrn=matgrn,
                      cp.types=NULL, cell_type = cell_type,
                      logged.data=FALSE, save.quant=TRUE, ncmp=ncmp, apply.fits=TRUE,
                      validate=FALSE)
 
     #To calculate measures of agreement before and  after normalization
     agreement(funtoonormout[[1]], individualID)
     agreement(funtoonormout[[2]], individualID)

 }
 
 
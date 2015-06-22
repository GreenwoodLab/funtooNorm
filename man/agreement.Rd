\name{agreement}
 \alias{agreement}
 
 \title{
 Function to access intra-replicate agreement in methylation data.
 }
 \description{
 Function calculates the sum of squared differences between sets of replicates, using
 methylation data, in order to evaluate the performance of the funtooNorm function.
 }
 \usage{
 agreement(Beta, patientid)
 }
 %- maybe also 'usage' for other objects documented here.
 \arguments{
   \item{Beta}{
 Matrix with beta-values, rows corresponding to probes, columns corresponding to samples. 
 }
  \item{patientid}{
 Vector of patient ID's matching the order of the samples (columns) in Beta. 
 Only patientid's that occur more than once (i.e. replicates) are used in calculating agreement.
 }
 }
 \details{
 We expect that the values returned by the agreement function after normalization by funtooNorm to be smaller than before.
 }
 \value{
 Function returns a measure of agreement between biological replicates in
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
     funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample,
 Annot=NULL,
                       controlred=matred, controlgrn=matgrn,
                       cp.types=NULL, cell_type = cell_type,
                       ncmp=4, save.quant=TRUE, save.loess=TRUE, apply.loess=TRUE,
 validate=FALSE)
 
     #To calculate measures of agreement before and  after normalization
     agreement(funtoonormout[[1]], patientid)
     agreement(funtoonormout[[2]], patientid)

 }
 
 
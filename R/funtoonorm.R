################################################################################
## Package to normalize Illumina Infinium Human Methylation 
## 450 BeadChip (Illumina 450K) with multiple tissues or cell types.
## Signal tables are organized with columns for samples and row for positions
## Inversely quantiles columns have quantiles as columns and rows as samples


################################################################################
## construct control probe summaries, averages by type of control probe
## then specifically create columns by cell type
constructProbCovMat <- function(controlred, controlgrn, cp.types,cell_type){
    ## For one control type, return the mean signal intensity 
    ## per color per sample
    getMeanProbesIntensity <- function(x){
    cbind(colMeans(array(controlred[cp.types==x,])),
          colMeans(array(controlgrn[cp.types==x,])))}
    
    ## means for each of the 2 signals and the 15 type of control probe
    ## mat.by.ct is 30 columns
    mat.by.ct <- do.call(cbind,lapply(unique(cp.types),getMeanProbesIntensity))
    
    ## now add additional sets of 30 columns, each specific to one cell type
    ctl.covmat <- mat.by.ct
    for(ct in unique(cell_type)){
    copy=mat.by.ct
    copy[cell_type!=ct,]=0
    ctl.covmat <- cbind(ctl.covmat, copy)
    }
    return(ctl.covmat)
}


################################################################################
## This function  allow to get a list of quantile between  0 and 1 that are more
## densely distributed at the extremes
buildQuantileList <- function(nCpG){
    qntllist <- c(seq((0.5)/nCpG, 0.0009, 0.0001), seq(0.001, 0.009, 0.001),
                seq(0.01,0.99,0.002), seq(0.991,0.999,0.001),
                seq(0.9991 ,(nCpG-0.5)/nCpG,0.0001 ))
    qntllist <- c(qntllist, (nCpG - 0.5)/nCpG)  ## add one more at the top end!
    return(qntllist[order(qntllist)])
}


################################################################################
## This function  allow to calculate Beta for two consistent tables 
## to avoid dividing by zero Illumina uses 100!
## But for instance minfi package use 0 by default
calcBeta <- function(A,B,offset=100) {
    return((2^B-1)/(2^A + 2^B-2 + offset))
}



################################################################################
## Validation graphs allowing to choose the number of components
## look better with no more than 8 components
## ncv.fold correspond to 
plotValidate <- function(quantiles, qntllist, ctl.covmat, numcompmax,
                         signaltype, ncv.fold=10, type.fits="PCR"){

    if (type.fits=="PCR") method.mvr <- "svdpc"
    else if (type.fits=="PLS") method.mvr <- "kernelpls"
    else stop("type.fits must be PCR or PLS")
    
    nqnt=ncol(quantiles)
    ns=nrow(quantiles)
    # TODO determine a minimum number of sample
    ncv.fold=min(as.integer(ns/2),ncv.fold) 
    mat <- array(NA, dim=c(ns, ncol=nqnt, numcompmax))
    ## potential for speedup by using an apply function here somehow?   
    for (j in (1:ncol(quantiles)))  {
      cvindex <- sample(ceiling((1:ns)/(ns/ncv.fold)))   
      for (k in (1:ncv.fold)) {
        wh.cv <- which(cvindex==k) ## Here is a list of unique samples
        tempfit <- pls::mvr(quantiles[-wh.cv,j] ~ ctl.covmat[-wh.cv,],
                       ncomp=numcompmax,
                       method = method.mvr)
        mat[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
      }
    }
    rmse <- t(apply(mat, 3, function(x) sqrt(colSums((x-quantiles)^2))))
    
    sq <- 2; while (max(rmse[,sq])>1.3*max(rmse[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse)[sq:nqnt,]),
            type='l', lty=1, col=rainbow(numcompmax),
            xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse[,sq:nqnt])), cex=1.0,
         paste(type.fits,", ",signaltype,sep=" "))
}


################################################################################
## on a numeric matrix return the quantile normalized matrix
## http://davetang.org/muse/2014/07/07/quantile-normalisation-in-r 
quantileNormalization <- function(mat){
    averagePerQuantiles=rowMeans(data.frame(apply(mat, 2, sort)))
    mat_ranked=apply(mat,2,rank,ties.method="min")
    return(apply(mat_ranked,2,function(x) averagePerQuantiles[x]))
}



################################################################################
## Main function of the package, apply a fit to a signal according to a number 
## of component 
##  !! TO DO  !! Add special treatment of Y chromosome.  See Fortin et al.
funtooNormApply <- function(signal, quantiles, qntllist,
                            ctl.covmat, ncmp=4, type.fits="PCR"){
  
    if (type.fits=="PCR") method.mvr <- "svdpc"
    else if (type.fits=="PLS") method.mvr <- "kernelpls"
    else stop("type.fits must be PCR or PLS")
    
    ## Fitting the model
    regression <- function(x) pls::mvr(x ~ ctl.covmat,
                                ncomp=ncmp,
                                method=method.mvr)$fitted.values[,1,ncmp]
    prediction <- apply(quantiles,2,regression)
    rankmat <- (apply(signal,2,  rank) - 0.5)/nrow(signal)
    predmat <- sapply(1:ncol(signal),function(x){
    stats::approx(qntllist,prediction[x,],xout=rankmat[,x])$y
    })
    return(predmat)
}

################################################################################
#' Function to measure intra-replicate agreement for methylation data.
#'
#' @param Beta : Matrix with beta-values, rows corresponding to probes, columns
#' corresponding to samples.  
#' @param individualID : a vector where 2 replicates have the exact same value
#' for two technical replicates. Order of samples should nmatch the samples
#' (columns) in Beta
#'
#' @details We expect that the values returned by the agreement function after
#'  normalization by funtooNorm to be smaller than before.
#' @return The average value of the square distance between replicates: 
#' a measure of agreement between replicates in methylation data.
#' @export
#'
#' @examples 
#' agreement(cbind(rnorm(n = 10),rnorm(n = 10),rnorm(n = 10)),c(1,1,1))
#' 
agreement <- function(Beta, individualID) {
    # List of duplicated values
    listOfReplicates=unique(individualID[duplicated(individualID)])
    listOfReplicates=listOfReplicates[!is.na(listOfReplicates)]
    sum.over.pt <- 0
    nbpairs=0
    for (id.rep in listOfReplicates)  {
    # Subtable of replicate
    mat.temp <- Beta[, which(individualID == id.rep)]
    # tall pair of columns for general case
    combs=utils::combn(1:ncol(mat.temp),2)
    # square differece beween each pair of value
    mat.diffs=apply(combs,2,function(x){(mat.temp[,x[1]]-mat.temp[,x[2]])^2})
    
    #Adding a score for each pair
    nbpairs=nbpairs+ncol(combs)
    sum.over.pt <- sum.over.pt + sum(colMeans(mat.diffs))
    }
    return(sum.over.pt)
}


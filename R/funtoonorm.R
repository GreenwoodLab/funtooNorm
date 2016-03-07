# Function to normalize Illumina Infinium Human Methylation 
# 450 BeadChip (Illumina 450K) with multiple tissues or cell types.




## construct control probe summaries, averages by type of control probe
## then specifically create columns by cell type
constructProbCovMat <- function(controlred, controlgrn, cp.types,cell_type){

  ## mat.by.ct is 30 columns - means for each type of control probe                          
  mat.by.ct <- do.call(cbind,lapply(unique(cp.types),function(x) cbind(colMeans(controlred[cp.types==x,]),colMeans(controlgrn[cp.types==x,]))))    
  
  ## now add additional sets of 30 columns, each specific to one cell type
  ctl.covmat <- mat.by.ct
  for(ct in rev(unique(cell_type))){
   copy=mat.by.ct
   copy[cell_type!=ct,]=0
   ctl.covmat <- cbind(ctl.covmat, copy)
  }
  return(ctl.covmat)
}


applyFits <- function(signal,qntlist,quantiles,regressionFunction){
  pred <- apply(quantiles,2,regressionFunction)
  rankmat<- (apply(signal,2,  rank) - 0.5)/nrow(signal)
  return(sapply(1:ncol(signal),function(x) approx(qntlist,pred[x,],xout=rankmat[,x])$y))
}


# Main function of the package

funtoonorm <- function(sigA, sigB, Annot, 
                       controlred, controlgrn, cp.types, cell_type, ncmp=4,  ncv.fold=10, 
                       logged2.data = FALSE, save.quant=TRUE, type.fits="PCR", validate=FALSE)
{
  
  ####################################################################################
  ## function to test if a matrix have non numeric values
  notNumeric<-function(a){
    return(any(!is.finite(as.matrix(a)))) 
  }
  ## to avoid dividing by zero Illumina uses 100!
  ## minfi package use 0 by default
  denom.offset <- 100   
  ##
  calcbeta <- function(A,B,offset) {
    return((2^B-1)/(2^A + 2^B-2 + offset))
  }
  
  
  if (is.null(Annot))  {
    message("Since Annot is NULL using default annotation.", '\n')
    data('Annot', package='funtooNorm', envir = environment())
  }

  
  #checking sanity of the data
  message("Checking sanity of the data...", '\n')
  if (notNumeric(sigA)){stop("There are non-numeric values in the matrix", '\n')}
  if (notNumeric(sigB)){stop("There are non-numeric values in the matrix", '\n')}
  if (notNumeric(controlred)){stop("There are non-numeric values in the matrix", '\n')}
  if (notNumeric(controlgrn)){stop("There are non-numeric values in the matrix", '\n')}
  if (any(cell_type == '' | typeof(levels(cell_type)) != 'character' | is.na(cell_type))) {{stop("There are non-character values in cell_type", '\n')}}
  if (length(unique(cell_type))<2) {{stop("There should be more that one tissue or cell type in cell_type variable", '\n')}}
  if (any(cp.types == '' | typeof(cp.types) != 'character' | is.na(cp.types))) {{stop("There are non-character values in cp.types", '\n')}}
  if (type.fits !="PCR" & type.fits !="PLS")  {{stop("type.fits must be either PCR or PLS","\n")}}
  
  #check that ID's match between sigA, sigB and control probe data
  if (identical(colnames(controlgrn), colnames(sigA)) &
      identical(colnames(controlred), colnames(sigA)) &
      identical(colnames(sigB), colnames(sigA))){
        orderednames <- order(colnames(controlred))
        controlgrn<-controlgrn[,orderednames]
        controlred<-controlred[,orderednames]
        sigA<-sigA[,orderednames]
        sigB<-sigB[,orderednames]
      }else{
        stop("Sample names do not match.", '\n')
      }
      
  if (type.fits!="PCR" & type.fits!="PLS") {stop("type.fits must be either PCR or PLS", "\n")}
  if (type.fits=="PCR") method.mvr <- "svdpc"
  if (type.fits=="PLS") method.mvr <- "kernelpls"
  
  
  
  if ((max(controlred)<20 & max(sigA)>25) | (max(controlred)>25 & max(sigA)<20)) {
    stop("apparent inconsistency w.r.t. log transformation of sigA/B data and control data \n")
  }
  if (max(sigA, na.rm=TRUE)>25) {
    message("Assuming data have not been previously log2 transformed, and applying a log2 transformation, \n")
    logged2.data <- FALSE
  }
  
  #check if annotation is ok, and subset annotation on the probes from Sig files
  if(any(!(rownames(sigA) %in% rownames(Annot)))){stop("Probe names do not match annotation entries", '\n')}
  cpg<-intersect(rownames(sigA), rownames(Annot))
  Annot=Annot[cpg,]
  
  if (!logged2.data) {        ## log transformation if data are not already logged at input
    sigA <- log2(1 + sigA)
    sigB <- log2(1 + sigB)
    controlgrn <- log2(1 + controlgrn)
    controlred <- log2(1 + controlred)
  }
  
  message("Data is ok.", '\n')
  
  wh.red <- which(Annot$Type=='I' & Annot$Color=="Red")
  wh.grn <- which(Annot$Type=='I' & Annot$Color=="Grn")
  wh.II <- which(Annot$Type=='II')
  
  ####! TO DO good to add special treatment of Y chromosome.  See Fortin et al. particularly for Y
  ## columns of quantilesA,B are the desired quantiles.  Rows are samples
  if (save.quant)  {
    
    
    qntllist <- c(seq((0.5)/nr, 0.0009, 0.0001), seq(0.001, 0.009, 0.001), seq(0.01,0.99,0.002), 
                  seq(0.991,0.999,0.001), seq(0.999 ,(nr-0.5)/nr,0.0001 ))
    qntllist <- c(qntllist, (nr - 0.5)/nr)  ## add one more at the top end!
    nqnt <- length(qntllist)  ## number of desired quantiles
    
    
    
    ## *temp* seem to take 20% less of time on the small data set
    quantilesA.red <- colQuantiles(sigA[wh.red,],prob=qntllist) 
    quantilesA.grn <- colQuantiles(sigA[wh.grn,],prob=qntllist)
    quantilesA.II <-  colQuantiles(sigA[wh.II,],prob=qntllist)
    quantilesB.red <- colQuantiles(sigB[wh.red,],prob=qntllist)
    quantilesB.grn <- colQuantiles(sigB[wh.grn,],prob=qntllist)
    quantilesB.II <-  colQuantiles(sigB[wh.II,],prob=qntllist)
    
    save(list=c("quantilesA.red","quantilesB.red",
                "quantilesA.grn","quantilesB.grn",
                "quantilesA.II","quantilesB.II"),file="quantiles.RData")

    ctl.covmat <- constructProbCovMat(controlred,controlgrn,cp.types,cell_type)

  } else  {
    load("quantiles.RData")
    #load("svd.ctlcovmat.RData")
  }
  
  
  
  
  
  
  
  
  ######################################################################################
  #validation chunk
  if (validate){
    if (notNumeric(numcomp) | numcomp>20)  {
      stop('There may be a problem with the SVD decomposition of the control matrix: number of components needed is either>20 or missing', "\n")
    }
    
    message('\n', 'Starting validation with ', numcomp , ' components...',  '\n')
    
    svd.ctlcovmat=svd(ctl.covmat)
    numcomp <- min(c(which(cumsum(svd.ctlcovmat$d)/sum(svd.ctlcovmat$d)>=0.98),8))  ## set max to 8
    ## validation graphs look better with no more than 8
    pcs.ctlcovmat <- svd.ctlcovmat$u[,1:numcomp]
    

    plotRMSE <- function(quant,ncomp,nfold,covmat,meth,label){
      mat <- array(NA, dim=c(nrow(quant), ncol=ncol(quant), ncomp))
      ## potential for speedup
      ## by using an apply function here somehow?   
      for (j in (1:ncol(quant)))  {
        cvindex <- sample(ceiling((1:nc)/(nrow(quant)/nfold)))   ## cross validation groups
        for (k in (1:nfold)) {
          wh.cv <- which(cvindex==k) ## Here is a list of unique samples
          tempfit <- mvr(quant[-wh.cv,j] ~ covmat[-wh.cv,],ncomp=ncomp, method = meth)
          mat[wh.cv, j, ] <- predict(tempfit, newdat = covmat[wh.cv,])
        }
      }
      rmse <- t(apply(mat, 3, function(x) sqrt(colSums((x-quant)^2))))
      
      sq <- 2; while (max(rmse[,sq])>1.3*max(rmse[,0.4*nqnt])) sq<-sq+1
      matplot(qntllist[sq:nqnt], sqrt(t(rmse)[sq:nqnt,]), type='l', lty=1, col=rainbow(numcomp), xlab = 'Percentiles', ylab='RMSE')
      text(0.6, max(sqrt(rmse[,sq:nqnt])), label, cex=1.0)
      return(mat)
    }

    ############# plot
    pdfFile=paste("validationcurves",type.fits,"pdf",sep='.')
    pdf(file = pdfFile, height=7, width=7)
    
    par(omi=c(0,0,0,0)) 
    m <- matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
    layout(mat = m, heights = c(0.3,0.3,0.15))
    par(omi=c(0,0,0,0),mar = c(2,2,0.2, 0.2), mgp = c(1,0,0))
    
    plotRMSE(quantilesA.red,numcomp,ncv.fold,ctl.covmat,method.mvr,paste(type.fits,", A I red",sep=""))
    plotRMSE(quantilesA.grn,numcomp,ncv.fold,ctl.covmat,method.mvr,paste(type.fits,", A I grn",sep=""))
    plotRMSE(quantilesA.II,numcomp,ncv.fold,ctl.covmat,method.mvr,paste(type.fits,", A II",sep=""))
    plotRMSE(quantilesB.red,numcomp,ncv.fold,ctl.covmat,method.mvr,paste(type.fits,", B I red",sep=""))
    plotRMSE(quantilesB.grn,numcomp,ncv.fold,ctl.covmat,method.mvr,paste(type.fits,", B I grn",sep=""))
    plotRMSE(quantilesB.II,numcomp,ncv.fold,ctl.covmat,method.mvr,paste(type.fits,", B II",sep=""))
    
    
    par(mar=c(0, 0, 1, 0)) 
    plot.new()
    ltylist <- rep(1, min(numcomp,8))
    if (numcomp>8) ltylist <- c(ltylist, rep(2, numcomp-8))
    legend(title='Number of components:', x = "top",inset = 0,
           legend = 1:numcomp, col=rainbow(numcomp), lty=ltylist, cex=1.1, horiz = TRUE)
    
    dev.off()
    ###############
    message('Done. Check your working directory for the file ',pdfFile,'\n')
    return()
    }else{
  
      ## could add in here different values of ncmp for different probe types  
      
      regression <- function(x) return(mvr(x ~ ctl.covmat, ncomp=ncmp, method=method.mvr)$fitted.values[,1,ncmp])
      
      ## applying the normalisation
    
      predmatA <- sigA
      predmatB <- sigB
      predmatA[wh.red,]=applyFits(sigA[wh.red,],qntllist,quantilesA.red,regression)
      predmatA[wh.grn,]=applyFits(sigA[wh.grn,],qntllist,quantilesA.grn,regression)
      predmatA[wh.II,]=applyFits(sigA[wh.II,],qntllist,quantilesA.II,regression)
      predmatB[wh.red,]=applyFits(sigB[wh.red,],qntllist,quantilesB.red,regression)
      predmatB[wh.grn,]=applyFits(sigB[wh.grn,],qntllist,quantilesB.grn,regression)
      predmatB[wh.II,]=applyFits(sigB[wh.II,],qntllist,quantilesB.II,regression)
      
      origBeta <- calcbeta(sigA, sigB, denom.offset)
      newBeta <- calcbeta(predmatA, predmatB, denom.offset)
      
      return(list(origBeta=origBeta, newBeta=newBeta))
  }
}  ## end function




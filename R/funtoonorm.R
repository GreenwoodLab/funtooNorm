# Function to normalize Illumina Infinium Human Methylation 450 BeadChip (Illumina 450K) with multiple tissues or cell types.
#Main function of the package.s
funtoonorm <- function(sigA, sigB, Annot=NULL, 
                       controlred, controlgrn, cp.types=NULL, cell_type, ncmp=4,  ncv.fold=10, 
                       logged.data = FALSE, save.quant=TRUE, type.fits="PCR", apply.fits=TRUE, validate=FALSE)
{
  ####################################################################################
  # functions
  logitfn <- function(x) { log(x/(1-x)) }
  ##
  mult1 <- function(x, vec1) return(x*vec1)
  ##
  sum1 <- function(x, v2) { return(x + v2) }
  ##
  calcbeta <- function(A,B,offset) {
    return((exp(B)-1)/(exp(A) + exp(B)-2 + offset))
  }
  extractqnt <- function(x, i, AB, ncmp)  { 
    return(x$fitted.values[i,AB,ncmp])  }
  ##
  NotNumeric<-function(a){
    return(any(!is.finite(as.matrix(a)))) 
  }
  ##
  is.wholenumber <-  function(x){abs(x - round(x)) == 0}
  ##
  denom.offset <- 100   # to avoid dividing by zero  Illumina uses 100!
  
  #Function adjusts beta values with results of loess regression
  applyfuntoonorm <- function(lfits, sigA, sigB, Annot, qntllist)  {
    
    wh.red <- which(Annot$Type=='I' & Annot$Color=="Red")
    wh.grn <- which(Annot$Type=='I' & Annot$Color=="Grn")
    wh.II <- which(Annot$Type=='II')
    
    rankmatA.red <- (apply(sigA[wh.red,],2,  rank) - 0.5)/length(wh.red)
    rankmatA.grn <- (apply(sigA[wh.grn,],2,  rank) - 0.5)/length(wh.grn)
    rankmatA.II <- (apply(sigA[wh.II,],2,  rank) - 0.5)/length(wh.II)
    
    rankmatB.red <- (apply(sigB[wh.red,],2,  rank) - 0.5)/length(wh.red)
    rankmatB.grn <- (apply(sigB[wh.grn,],2,  rank) - 0.5)/length(wh.grn)
    rankmatB.II <- (apply(sigB[wh.II,],2,  rank) - 0.5)/length(wh.II)
    
    predmatA.red <- matrix(NA, nrow=nrow(rankmatA.red), ncol = ncol(rankmatA.red))
    predmatA.grn <- matrix(NA, nrow=nrow(rankmatA.grn), ncol = ncol(rankmatA.grn))
    predmatA.II <- matrix(NA, nrow=nrow(rankmatA.II), ncol = ncol(rankmatA.II))
    predmatB.red <- matrix(NA, nrow=nrow(rankmatB.red), ncol = ncol(rankmatB.red))
    predmatB.grn <- matrix(NA, nrow=nrow(rankmatB.grn), ncol = ncol(rankmatB.grn))
    predmatB.II <- matrix(NA, nrow=nrow(rankmatB.II), ncol = ncol(rankmatB.II))
    for (i in (1:ncol(predmatA.red)))  {
      predmatA.red[,i]<-approx(qntllist,lfits[[1]][i,],xout=rankmatA.red[,i])$y
      predmatA.grn[,i]<-approx(qntllist,lfits[[2]][i,],xout=rankmatA.grn[,i])$y
      predmatA.II[,i]<-approx(qntllist,lfits[[3]][i,],xout=rankmatA.II[,i])$y
      predmatB.red[,i]<-approx(qntllist,lfits[[4]][i,],xout=rankmatB.red[,i])$y
      predmatB.grn[,i]<-approx(qntllist,lfits[[5]][i,],xout=rankmatB.grn[,i])$y
      predmatB.II[,i]<-approx(qntllist,lfits[[6]][i,],xout=rankmatB.II[,i])$y
    }
    predmatA <- matrix(NA, nrow(sigA), ncol(predmatA.red))
    predmatB <- predmatA
    predmatA[wh.red,] <- predmatA.red
    predmatA[wh.grn,] <- predmatA.grn
    predmatA[wh.II,] <- predmatA.II
    predmatB[wh.red,] <- predmatB.red
    predmatB[wh.grn,] <- predmatB.grn
    predmatB[wh.II,] <- predmatB.II
    
    newBeta <- calcbeta(predmatA, predmatB, denom.offset)        
    origBeta <- as.matrix(calcbeta(sigA, sigB, denom.offset))
    rownames(newBeta) <- rownames(origBeta)
    colnames(newBeta) <- colnames(origBeta)
    return(list(origBeta=origBeta, newBeta=newBeta))
  }
  ###############################################################################
  
  
  if (is.null(Annot))  {
    cat("Since Annot is NULL using default annotation.", '\n')
    data('Annot', package='funtooNorm', envir = environment())
  }
  
  # obtain names of types of control probes if not provided
  if (!is.null(cp.types)) { cp.types.tab <- table(cp.types) }
  if (is.null(cp.types))  {
    cp.types.tab <- table(rownames(controlgrn))
    if (any(names(cp.types.tab)=="NEGATIVE")) {
      cp.types <- rownames(controlgrn)
    }
    if (!(any(names(cp.types.tab)=="NEGATIVE")))  {
      cat("Since cp.types is NULL and rownames of control probe matrices do not contain,", "\n", 
          "this information, analysis will use default cp.types.", '\n')
      data('cp.types', package='funtooNorm', envir = environment())
      cp.types.tab <- table(cp.types)
    } 
  }
  
  #checking sanity of the data
  cat("Checking sanity of the data...", '\n')
  if (NotNumeric(sigA)){stop("There are non-numeric values in the matrix", '\n')}
  if (NotNumeric(sigB)){stop("There are non-numeric values in the matrix", '\n')}
  if (NotNumeric(controlred)){stop("There are non-numeric values in the matrix", '\n')}
  if (NotNumeric(controlgrn)){stop("There are non-numeric values in the matrix", '\n')}
  if (any(cell_type == '' | typeof(levels(cell_type)) != 'character' | is.na(cell_type))) {{stop("There are non-character values in cell_type", '\n')}}
  if (length(unique(cell_type))<2) {{stop("There should be more that one tissue or cell type in cell_type variable", '\n')}}
  if (any(cp.types == '' | typeof(cp.types) != 'character' | is.na(cp.types))) {{stop("There are non-character values in cp.types", '\n')}}
  if (type.fits !="PCR" & type.fits !="PLS")  {{stop("type.fits must be either PCR or PLS","\n")}}
  
  #check that ID's match between sigA, sigB and control probe data
  if (!(identical(colnames(controlgrn), colnames(sigA)) &
        identical(colnames(controlred), colnames(sigA)) &
        identical(colnames(sigB), colnames(sigA)))) {
    if ((ncol(controlgrn) == ncol(sigA)) & (ncol(controlred) == ncol(sigA)) 
        & (ncol(controlred) == ncol(sigB)) & (length(cp.types) == ncol(sigB)) & (length(cell_type) == ncol(sigB))){
      ix <- sort(colnames(controlred), index.return = TRUE)$ix
      controlgrn<-controlgrn[,ix]; sigA<-sigA[,ix]; sigB<-sigB[,ix];
      if (!(identical(colnames(controlgrn), colnames(sigA)) &
            identical(colnames(controlred), colnames(sigA)) &
            identical(colnames(sigB), colnames(sigA)))){
        stop("Sample names do not match.", '\n')}
      
    } else{stop("Data dimensions or samples names do not match.", '\n')}
  }
  if (type.fits!="PCR" & type.fits!="PLS") {stop("type.fits must be either PCR or PLS", "\n")}
  
  if ((max(controlred)<20 & max(sigA)>25) | (max(controlred)>25 & max(sigA)<20)) {
    stop("apparent inconsistency w.r.t. log transformation of sigA/B data and control data \n")
  }
  if (max(sigA, na.rm=TRUE)>25) {
    cat("Assuming data have not been previously log transformed, and applying a log transformation, \n")
    logged.data <- FALSE
  }
  
  
  #check if annotation is ok, and subset annotation on the probes from Sig files
  if(any(!(rownames(sigA) %in% rownames(Annot)))){stop("Probe names do not match annotation entries", '\n')}
  cpg<-intersect(rownames(sigA), rownames(Annot))
  Annot=Annot[cpg,]
  
  cat("Data is ok.", '\n')
  
  if (!logged.data) {        # log transformation if data are not already logged at input
    sigA <- log2(1 + sigA)
    sigB <- log2(1 + sigB)
  }
  
  wh.red <- which(Annot$Type=='I' & Annot$Color=="Red")
  wh.grn <- which(Annot$Type=='I' & Annot$Color=="Grn")
  wh.II <- which(Annot$Type=='II')
  
  nr <- nrow(sigA)  
  qntllist <- c(seq((0.5)/nr, 0.0009, 0.0001), seq(0.001, 0.009, 0.001), seq(0.01,0.99,0.002), 
                seq(0.991,0.999,0.001), seq(0.999 ,(nr-0.5)/nr,0.0001 ))
  qntllist <- c(qntllist, (nr - 0.5)/nr)  # add one more at the top end!
  nqnt <- length(qntllist)  # number of desired quantiles
  
  if (type.fits=="PCR") method.mvr <- "svdpc"
  if (type.fits=="PLS") method.mvr <- "kernelpls"
  
  ####! TO DO good to add special treatment of Y chromosome.  See Fortin et al. particularly for Y
  ## columns of quantilesA,B are the desired quantiles.  Rows are samples
  if (save.quant)  {
    quantilesA.red <- matrix(NA, ncol(sigA), nqnt)   
    quantilesA.grn <- quantilesA.red;  quantilesA.II <- quantilesA.red
    quantilesB.red <- matrix(NA, ncol(sigA), nqnt)  
    quantilesB.grn <- quantilesA.red; quantilesB.II <- quantilesA.red      
    for (i in (1:ncol(sigA)))  {
      quantilesA.red[i,] <- quantile(sigA[wh.red,i],qntllist) 
      quantilesA.grn[i,] <- quantile(sigA[wh.grn,i],qntllist)
      quantilesA.II[i,] <-  quantile(sigA[wh.II,i],qntllist)
      quantilesB.red[i,] <- quantile(sigB[wh.red,i],qntllist)
      quantilesB.grn[i,] <- quantile(sigB[wh.grn,i],qntllist)
      quantilesB.II[i,] <-  quantile(sigB[wh.II,i],qntllist)
    }
    save(quantilesA.red, file="quantilesA.red.RData")
    save(quantilesA.grn, file="quantilesA.grn.RData")
    save(quantilesA.II, file="quantilesA.II.RData")
    save(quantilesB.red, file="quantilesB.red.RData")
    save(quantilesB.grn, file="quantilesB.grn.RData")
    save(quantilesB.II, file="quantilesB.II.RData")
  }  # if save.quant
  if (!save.quant)  {
    load("quantilesA.red.RData")
    load("quantilesA.grn.RData")
    load("quantilesA.II.RData")
    load("quantilesB.red.RData")
    load("quantilesB.grn.RData")
    load("quantilesB.II.RData")
  }
  
  if (!logged.data) {   # same assumption as for signalA, signalB w.r.t. log transformation
    controlgrn <- log2(1 + controlgrn)
    controlred <- log2(1 + controlred)
  }
  
  # construct control probe summaries, averages by type of control probe
  # then specifically create columns by cell type
  for (k in (1:length(cp.types.tab))) {
    temp1 <- apply(controlred[cp.types==names(cp.types.tab[k]),],2,mean) 
    temp2 <- apply(controlgrn[cp.types==names(cp.types.tab[k]),],2,mean)
    if (k==1) {
      mat.by.ct <- cbind(temp1,temp2)  }
    if (k>1) {  mat.by.ct <- cbind(mat.by.ct, cbind(temp1,temp2))  }
  } # end for   mat.by.ct is 15 columns - means for each type of control probe
  
  ind=matrix(FALSE,ncol=length(cell_type),nrow=length(unique(cell_type)))
  for(j in 1:length(unique(cell_type))){
    ind[j,]<- (cell_type==names(table(cell_type))[j])
  }
  # now add additional sets of 15 columns, each specific to one cell type
  ctl.covmat <- cbind(mat.by.ct, apply(mat.by.ct, 2, mult1, ind[1,]))
  for(j in 2:length(unique(cell_type))){
    ctl.covmat <- cbind(ctl.covmat, apply(mat.by.ct, 2, mult1, ind[j,]))
  }
  
  svd.ctlcovmat <- svd(ctl.covmat)
  numcomp <- min(c(which(cumsum(svd.ctlcovmat$d)/sum(svd.ctlcovmat$d)>=0.98),8))  # set max to 8
  # validation graphs look better with no more than 8
  pcs.ctlcovmat <- svd.ctlcovmat$u[,1:numcomp]
  
  if (save.quant)  {
    save(svd.ctlcovmat, file="svd.ctlcovmat.RData")
  }
  if (!save.quant)  {
    load("svd.ctlcovmat.RData")
  }
  
  ######################################################################################
  #validation chunk
  if (validate){
    
    if (NotNumeric(numcomp) | numcomp>20)  {
      stop('There may be a problem with the SVD decomposition of the control matrix: number of components needed is either>20 or missing', "\n")
    }
    
    cat('\n', 'Starting validation with ', numcomp , ' components...',  '\n')
    
    cores=1 #probably does not make sense to use parallelisation
    pls.options(parallel = cores)
    
    ns <- nrow(quantilesA.red)
    mat1P.Ared <- array(NA, dim=c(ns, ncol=nqnt, numcomp))
    mat1P.Agrn <- mat1P.Ared;  mat1P.AII <- mat1P.Ared;     mat1P.Bred <- mat1P.Ared
    mat1P.Bgrn <- mat1P.Ared;  mat1P.BII <- mat1P.Ared;     mat2P.Ared <- mat1P.Ared; 
    mat2P.Agrn <- mat1P.Ared;  mat2P.AII <- mat1P.Ared;     mat2P.Bred <- mat1P.Ared; 
    mat2P.Bgrn <- mat1P.Ared;  mat2P.BII <- mat1P.Ared; 
    
    # potential for speedup
    # by using an apply function here somehow?   
    for (j in (1:nqnt))  {
      cvindex <- sample(ceiling((1:ns)/(ns/ncv.fold)))   # cross validation groups
      for (k in (1:10)) {
        wh.cv <- which(cvindex==k)
        tempfit <- mvr(quantilesA.red[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "svdpc")
        mat1P.Ared[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesA.grn[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "svdpc")
        mat1P.Agrn[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesA.II[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "svdpc")
        mat1P.AII[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesB.red[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "svdpc")
        mat1P.Bred[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesB.grn[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "svdpc")
        mat1P.Bgrn[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesB.II[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "svdpc")
        mat1P.BII[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        
        tempfit <- mvr(quantilesA.red[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "kernelpls")
        mat2P.Ared[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesA.grn[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "kernelpls")
        mat2P.Agrn[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesA.II[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "kernelpls")
        mat2P.AII[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesB.red[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "kernelpls")
        mat2P.Bred[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesB.grn[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "kernelpls")
        mat2P.Bgrn[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
        tempfit <- mvr(quantilesB.II[-wh.cv,j] ~ ctl.covmat[-wh.cv,],ncomp=numcomp, method = "kernelpls")
        mat2P.BII[wh.cv, j, ] <- predict(tempfit, newdat = ctl.covmat[wh.cv,])
      }  # for k
    }  # for j
    
    rmse1.Ared <- matrix(NA, nrow = numcomp, ncol=nqnt)
    rmse1.Agrn <- rmse1.Ared; rmse1.AII <- rmse1.Ared
    rmse2.Ared <- rmse1.Ared;  
    rmse2.Agrn <- rmse1.Ared; rmse2.AII <- rmse1.Ared
    rmse1.Bred <- rmse1.Ared
    rmse1.Bgrn <- rmse1.Ared; rmse1.BII <- rmse1.Ared
    rmse2.Bred <- rmse1.Ared;  
    rmse2.Bgrn <- rmse1.Ared; rmse2.BII <- rmse1.Ared
    
    for (i.n in (1:numcomp)) {
      rmse1.Ared[i.n,] <- sqrt(apply((mat1P.Ared[,,i.n]-quantilesA.red)^2, 2, sum))
      rmse1.Agrn[i.n,] <- sqrt(apply((mat1P.Agrn[,,i.n]-quantilesA.grn)^2, 2, sum))
      rmse1.AII[i.n,] <- sqrt(apply((mat1P.AII[,,i.n]-quantilesA.II)^2, 2, sum))
      rmse1.Bred[i.n,] <- sqrt(apply((mat1P.Bred[,,i.n]-quantilesB.red)^2, 2, sum))
      rmse1.Bgrn[i.n,] <- sqrt(apply((mat1P.Bgrn[,,i.n]-quantilesB.grn)^2, 2, sum))
      rmse1.BII[i.n,] <- sqrt(apply((mat1P.BII[,,i.n]-quantilesB.II)^2, 2, sum))
      
      rmse2.Ared[i.n,] <- sqrt(apply((mat2P.Ared[,,i.n]-quantilesA.red)^2, 2, sum))
      rmse2.Agrn[i.n,] <- sqrt(apply((mat2P.Agrn[,,i.n]-quantilesA.grn)^2, 2, sum))
      rmse2.AII[i.n,] <- sqrt(apply((mat2P.AII[,,i.n]-quantilesA.II)^2, 2, sum))
      rmse2.Bred[i.n,] <- sqrt(apply((mat2P.Bred[,,i.n]-quantilesB.red)^2, 2, sum))
      rmse2.Bgrn[i.n,] <- sqrt(apply((mat2P.Bgrn[,,i.n]-quantilesB.grn)^2, 2, sum))
      rmse2.BII[i.n,] <- sqrt(apply((mat2P.BII[,,i.n]-quantilesB.II)^2, 2, sum))
    }
    
    pdf(file = 'validationcurves.PCR.pdf', height=7, width=7)
    
    par(omi=c(0,0,0,0)) 
    m <- matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
    layout(mat = m, heights = c(0.3,0.3,0.15))
    
    par(mar = c(2,2,0.2, 0.2), mgp = c(1,0,0))
    sq <- 2; while (max(rmse1.Ared[,sq])>1.3*max(rmse1.Ared[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse1.Ared)[sq:nqnt,]), type='l', lty=1, col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse1.Ared[,sq:nqnt])), "PCR, A I red", cex=1.0)
    
    sq <- 2; while (max(rmse1.Agrn[,sq])>1.3*max(rmse1.Agrn[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse1.Agrn)[sq:nqnt,]), type='l', lty=1,col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse1.Agrn[,sq:nqnt])), "PCR, A I green", cex=1.0)
    
    sq <- 2; while (max(rmse1.AII[,sq])>1.3*max(rmse1.AII[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse1.AII)[sq:nqnt,]), type='l', lty=1, col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse1.AII[,sq:nqnt])), "PCR, A II", cex=1.0)
    
    sq <- 2; while (max(rmse1.Bred[,sq])>1.3*max(rmse1.Bred[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse1.Bred)[sq:nqnt,]), type='l', lty=1, col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse1.Bred[,sq:nqnt])), "PCR, B I red", cex=1.0)
    
    sq <- 2; while (max(rmse1.Bgrn[,sq])>1.3*max(rmse1.Bgrn[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse1.Bgrn)[sq:nqnt,]), type='l', lty=1,col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse1.Bgrn[,sq:nqnt])), "PCR, B I green", cex=1.0)
    
    sq <- 2; while (max(rmse1.BII[,sq])>1.3*max(rmse1.BII[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse1.BII)[sq:nqnt,]), type='l', lty=1, col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse1.BII[,sq:nqnt])), "PCR, B II", cex=1.0)
    
    par(mar=c(0, 0, 1, 0)) 
    plot.new()
    ltylist <- rep(1, min(numcomp,8))
    if (numcomp>8) ltylist <- c(ltylist, rep(2, numcomp-8))
    legend(title='Number of components:', x = "top",inset = 0,
           legend = 1:numcomp, col=1:numcomp, lty=ltylist, cex=1.1, horiz = TRUE)
    
    dev.off()
    
    pdf(file = 'validationcurves.PLS.pdf', height=7, width=7)
    par(omi=c(0,0,0,0)) 
    m <- matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
    layout(mat = m, heights = c(0.3,0.3,0.15))
    
    par(mar = c(2,2,0.2, 0.2), mgp = c(1,0,0))
    sq <- 2; while (max(rmse2.Ared[,sq])>1.3*max(rmse2.Ared[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse2.Ared)[sq:nqnt,]), type='l', lty=1, col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse2.Ared[,sq:nqnt])), "PLS, A I red", cex=1.0)
    
    sq <- 2; while (max(rmse2.Agrn[,sq])>1.3*max(rmse2.Agrn[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse2.Agrn)[sq:nqnt,]), type='l', lty=1,col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse2.Agrn[,sq:nqnt])), "PLS, A I green", cex=1.0)
    
    sq <- 2; while (max(rmse2.AII[,sq])>1.3*max(rmse2.AII[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse2.AII)[sq:nqnt,]), type='l', lty=1, col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse2.AII[,sq:nqnt])), "PLS, A II", cex=1.0)
    
    sq <- 2; while (max(rmse2.Bred[,sq])>1.3*max(rmse2.Bred[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse2.Bred)[sq:nqnt,]), type='l', lty=1, col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse2.Bred[,sq:nqnt])), "PLS, B I red", cex=1.0)
    
    sq <- 2; while (max(rmse2.Bgrn[,sq])>1.3*max(rmse2.Bgrn[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse2.Bgrn)[sq:nqnt,]), type='l', lty=1,col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse2.Bgrn[,sq:nqnt])), "PLS, B I green", cex=1.0)
    
    sq <- 2; while (max(rmse2.BII[,sq])>1.3*max(rmse2.BII[,0.4*nqnt])) sq<-sq+1
    matplot(qntllist[sq:nqnt], sqrt(t(rmse2.BII)[sq:nqnt,]), type='l', lty=1, col=1:numcomp, xlab = 'Percentiles', ylab='RMSE')
    text(0.6, max(sqrt(rmse2.BII[,sq:nqnt])), "PLS, B II", cex=1.0)
    
    par(mar=c(0, 0, 1, 0)) 
    plot.new()
    ltylist <- rep(1, min(numcomp,8))
    if (numcomp>8) ltylist <- c(ltylist, rep(2, numcomp-8))
    legend(title='Number of components:', x = "top",inset = 0,
           legend = 1:numcomp, col=1:numcomp, lty=ltylist, cex=1.1, horiz = TRUE)
    dev.off()
    
    selquant <- sort(sample((1:nqnt), 10))
    rmse.red <- array(NA, dim = c(length(selquant), length(selquant), numcomp))
    betaP.red <- rmse.red
    rmse.grn <- rmse.red;   betaP.grn <- betaP.red
    rmse.II <- rmse.red;  betaP.II <- betaP.red
    for (i.n in (1:numcomp)) {
      for (i.s1 in (1:length(selquant))) {
        for (i.s2 in (1:length(selquant))) {
          betaP <- calcbeta(mat1P.Ared[,selquant[i.s1], i.n], mat1P.Bred[,selquant[i.s2],i.n],offset=0.0001)
          betaO <- calcbeta(quantilesA.red[,selquant[i.s1]], quantilesB.red[,selquant[i.s2]], offset=0.0001)
          rmse.red[i.s1,i.s2, i.n] <- sqrt(sum((betaP-betaO)^2))
          betaP.red[i.s1,i.s2, i.n] <- mean(betaP)
          betaP <- calcbeta(mat1P.Agrn[,selquant[i.s1], i.n], mat1P.Bgrn[,selquant[i.s2],i.n],offset=0.0001)
          betaO <- calcbeta(quantilesA.grn[,selquant[i.s1]], quantilesB.grn[,selquant[i.s2]], offset=0.0001)
          rmse.grn[i.s1,i.s2, i.n] <- sqrt(sum((betaP-betaO)^2))
          betaP.grn[i.s1,i.s2, i.n] <- mean(betaP)
          betaP <- calcbeta(mat1P.AII[,selquant[i.s1], i.n], mat1P.BII[,selquant[i.s2],i.n],offset=0.0001)
          betaO <- calcbeta(quantilesA.II[,selquant[i.s1]], quantilesB.II[,selquant[i.s2]], offset=0.0001)
          rmse.II[i.s1,i.s2, i.n] <- sqrt(sum((betaP-betaO)^2))
          betaP.II[i.s1,i.s2, i.n] <- mean(betaP)
        }}    # i.s1, i.s2
    }  # i.n
    
    #pdf(file = 'validationcurves.beta.PCR.pdf', height=7, width=7)
    #par(mfrow = c(2,2), omi=c(0,0,0,0), mar = rep(0,4), mgp = c(-1,0,0)) 
    
    #plot(c(0,1), c(min(rmse.red), max(rmse.red)), type='n', xlab='', ylab='RMSE')
    #for (i.n in (1:8)) {
    #  points(betaP.red[,,i.n], rmse.red[,,i.n], col=i.n, pch=16, cex=0.7)  }
    
    #plot(c(0,1), c(min(rmse.grn), max(rmse.grn)), type='n', xlab='', ylab='RMSE')
    #for (i.n in (1:8)) {
    #  points(betaP.grn[,,i.n], rmse.grn[,,i.n], col=i.n, pch=16, cex=0.7)  }
    
    #plot(c(0,1), c(min(rmse.II), max(rmse.II)), type='n', xlab='', ylab='RMSE')
    #for (i.n in (1:8)) {
    #  points(betaP.II[,,i.n], rmse.II[,,i.n], col=i.n, pch=16, cex=0.7)  }
    
    #plot.new()
    #legend(title='Number of components:', x = "top",
    #       legend = 1:numcomp, col=1:numcomp, cex=1.1, lty=1, horiz = FALSE)
    
    #dev.off()
    
    cat('Done. Check your working directory for the file "validationcurves.PCR.pdf" and validationcurves.PLS.pdf"', '\n')
    return()
  }    # end if validate
  #########################
  
  # could add in here different values of ncmp for different probe types  
  predA.red <- matrix(NA, nrow = ncol(sigA), ncol = nqnt)
  predA.grn <- matrix(NA, nrow = ncol(sigA), ncol = nqnt)
  predA.II <- matrix(NA, nrow = ncol(sigA), ncol = nqnt)
  predB.red <- matrix(NA, nrow = ncol(sigA), ncol = nqnt)
  predB.grn <- matrix(NA, nrow = ncol(sigA), ncol = nqnt)
  predB.II <- matrix(NA, nrow = ncol(sigA), ncol = nqnt)
  for (j in (1:nqnt)) {
    predA.red[,j] <- mvr(quantilesA.red[,j] ~ ctl.covmat, ncomp=ncmp, method=method.mvr)$fitted.values[,1,ncmp]
    predB.red[,j] <- mvr(quantilesB.red[,j] ~ ctl.covmat, ncomp=ncmp, method=method.mvr)$fitted.values[,1,ncmp]
    predA.grn[,j] <- mvr(quantilesA.grn[,j] ~ ctl.covmat, ncomp=ncmp, method=method.mvr)$fitted.values[,1,ncmp]
    predB.grn[,j] <- mvr(quantilesB.grn[,j] ~ ctl.covmat, ncomp=ncmp, method=method.mvr)$fitted.values[,1,ncmp]
    predA.II[,j] <- mvr(quantilesA.II[,j] ~ ctl.covmat, ncomp=ncmp, method=method.mvr)$fitted.values[,1,ncmp]
    predB.II[,j] <- mvr(quantilesB.II[,j] ~ ctl.covmat, ncomp=ncmp, method=method.mvr)$fitted.values[,1,ncmp]
  } # for j
  
  lfits <- list(predA.red, predA.grn, predA.II, 
                predB.red, predB.grn, predB.II)
  if (apply.fits)  {
    betamatrices <- applyfuntoonorm(lfits, sigA, sigB, Annot, qntllist)
    return(betamatrices)
  } else {   return(lfits)  }
}  # end function




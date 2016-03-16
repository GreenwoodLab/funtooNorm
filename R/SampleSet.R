################################################################################
## This file allow to create some 'SampleSet' objects, they are list containing
## signal data and different variables usefull for funTooNorm
## The data is separated into the 3 kind positions are quantified
## Each way have then 2 channel (methylated and unmethylated ie : A and B
## We define then the 6 (2*3) labels: AIGrn BIGrn AIRed BIRed AII BII
##
## Each SampleSet object are list containing :
## - type: is 'minfi' or 'GenomeStudio'
##     -'minfi' imply a field rgset containing the object
##     -'GenomeStudio' imply 2 fields containing the file names
## - sampleNames: conatain the list of sample names in the order used in all tables
## - sampleSize: is the number of samples
## - nPos is the number of positions in the chip
## - cell_type is a list matching each sample to define the different batches
## - signal: is the list of 6 tables containing the 6 type of signals
## - predmat: is the list of 6 table containing the normalize signals
## - Annotation: is a list of name in the same order as the postions in the tables
## - qntllist, quantiles, ctl.covmat are internal variable to avoid recomputing

################################################################################
## To create a sample set from RGSet from minfi package
## RGset should contain a cell_type vector in it s phenotypes data (pdata)
fromRGChannelSet <- function(argRGset){
  o <- list(type="minfi",
            rgset=argRGset)
  
  if(! "cell_type" %in% colnames(pData(o$rgset))){
    stop("Your object should contain a field \"cell_type\" in it's phenotype
         data pData(rgset) in order to use funTooNorm")
  }
  
  return(Initialize.SampleSet(o,pData(o$rgset)$cell_type))

}

################################################################################
## To create a sample set from GenomeStudio files
## control file and signal file should have the exact same samples
## cell_type vector should have names matching the 2 files from genome studios
fromGenStudFiles <- function(controlProbeFile,signalFile,cell_type,path=""){
  o <- list(type="genomeStudio",
            controlProbeFile=paste(path,controlProbeFile,sep=""),
            signalFile=paste(path,signalFile,sep=""))
  return(Initialize.SampleSet(o,cell_type))
}


################################################################################
## Print information about the SampleSet
print.SampleSet <- function(object){
  cat("SampleSet object built from ",object$type,'\n')
  cat("Signal data: ",object$nPos,"positions and ",object$sampleSize, "samples",'\n')
  cat("   cell type:",levels(object$cell_type),'\n')
  cat("  ",length(object$qntllist),"quantiles",'\n')
  if(is.null(object$predmat)){
    cat("funTooNorm Normalization was not applied",'\n')
  }else{
    cat("funTooNorm Normalization was applied",'\n')
  }
}



################################################################################
getAnnotation <- function(object){
  return(c(object$Annotation$IGrn,
         object$Annotation$IRed,
         object$Annotation$II))
}

################################################################################
getLogRawSigA <- function(object){
  return(rbind(object$signal$AIGrn,
               object$signal$AIRed,
               object$signal$AII))
}
getLogRawSigB <- function(object){
  return(rbind(object$signal$BIGrn,
               object$signal$BIRed,
               object$signal$BII))
}
getLogNormSigA <- function(object){
  if(any(is.na(c(object$predmat$AIGrn,object$predmat$AIRed,object$predmat$AII)))){
    message("WARNING: you cant t access signal A normalization, please call SampleSet=funTooNorm(SampleSet)")
  }else{
    return(rbind(object$predmat$AIGrn,
                 object$predmat$AIRed,
                 object$predmat$AII))
  }
}
getLogNormSigB <- function(object){
    if(any(is.na(c(object$predmat$BIGrn,object$predmat$BIRed,object$predmat$BII)))){
      message("WARNING: you cant t access signal B normalization, please call SampleSet=funTooNorm(SampleSet)")
    }else{
      return(rbind(object$predmat$BIGrn,
                   object$predmat$BIRed,
                   object$predmat$BII))    
    }
}

getRawBeta <- function(object,offset=100){
  return(calcBeta(getLogRawSigA(object),getLogRawSigB(object),offset))
}

getNormBeta <- function(object,offset=100){
  return(calcBeta(getLogNormSigA(object),getLogNormSigB(object),offset))
}


funTooNorm <- function(object, type.fits="PCR",ncmp=4,force=FALSE){
  
  if(force | is.null(object$predmat)){
    object$predmat=list()
  }else{
    message("Normalization already exist")
  }
  
  for(signal in names(object$quantiles)){
    if(force | is.null(object$predmat[[signal]])){
      message("Normalization of signal : ",signal)
      object$predmat[[signal]]=funtoonormApply(object$signal[[signal]],
                                               object$quantiles[[signal]],
                                               object$qntllist,
                                               object$ctl.covmat,
                                               ncmp,type.fits)
    }else{
      message("Already done : ",signal)
    }
  }
  return(object)
}


plotValidationGraph <- function(object, type.fits="PCR",file=""){
  
  
  svd.ctlcovmat=svd(object$ctl.covmat)
  numcomp <- min(c(which(cumsum(svd.ctlcovmat$d)/sum(svd.ctlcovmat$d)>=0.98),8))  ## set max to 8
  if (!is.finite(numcomp) | numcomp>20)  {
    stop('There may be a problem with the SVD decomposition of the control matrix: number of components needed is either>20 or missing', "\n")
  }
  message('\n', 'Starting validation with a max of ', numcomp , ' components...',  '\n')
  
  if(file!=""){
    message("plotting in the pdf file:",file)
    pdf(file = file, height=7, width=7)
  }
  
  layout(mat =  matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE),
         heights = c(0.3,0.3,0.15))
  par(omi=c(0,0,0,0),mar = c(2,2,0.2, 0.2), mgp = c(1,0,0))
  
  
  for(i in names(object$quantiles)){
    plotValidate(object$quantiles[[i]],object$qntllist,object$ctl.covmat,numcomp,i)
  }
  
  
  par(mar=c(0, 0, 1, 0)) 
  plot.new()
  ltylist <- rep(1, min(numcomp,8))
  if (numcomp) ltylist <- c(ltylist, rep(2, max(8,numcomp)-8))
  legend(title='Number of components:', x = "top",inset = 0,
         legend = 1:numcomp, col=rainbow(numcomp), lty=ltylist, cex=1.1, horiz = TRUE)
  if(file!=""){
    dev.off()
  }
  
}


################################################################################
## Common long function to load and format the data.
Initialize.SampleSet <- function(object,cell_type){
  if (any(cell_type == '' | is.na(cell_type))) {
    stop("There are NA values in cell_type", '\n')
  }
  if (length(unique(cell_type))<2) {
    stop("There should be more that one tissue or cell type in cell_type variable", '\n')
  }
  class(object) <- "SampleSet"
  object$sampleSize=length(cell_type)
  object$cell_type=cell_type
  
  
  ## Genome Studio Data Files
  if(object$type=="genomeStudio"){
    
    ## Construct the control Table
    controlTable=read.table(object$controlProbeFile,sep='\t',header=TRUE)
    object$sampleNames=unlist(strsplit(colnames(controlTable)[1:object$sampleSize*3+1],".Signal_Grn"))
    
    controlred=controlTable[,paste(object$sampleName,".Signal_Red",sep="")]
    controlgrn=controlTable[,paste(object$sampleName,".Signal_Grn",sep="")]
    cp.types=controlTable$TargetID
    
    if(any(! object$sampleNames %in% names(object$cell_type))){
      i=which(! object$sampleNames%in%names(object$cell_type))[1]
      stop("STOP: ",object$sampleNames[i],"from file ",object$controlProbeFile,
           " should be present in the cell_types names:\n")
    }
    object$cell_type=object$cell_type[object$sampleNames]
    
  }else if (object$type=="minfi") {
    object$sampleNames=rownames(pData(object$rgset))
    
    ## Formating the control probes data for the covariance matrix
    controlTable <- getProbeInfo(getManifest(object$rgset), type="Control")
    
    #Here we have to remove the 2 probes that are absent from the Chip
    controlTable=controlTable[! controlTable$Address %in% c("21630339","24669308"),]
    #Here we remove the 15 probes that cannot be in the GenomeStudio output
    controlTable=controlTable[controlTable$Color!="-99",]
    controlred=as.data.frame(getRed(object$rgset)[controlTable$Address,])
    
    controlgrn=as.data.frame(getGreen(object$rgset)[controlTable$Address,])
    
    cp.types=controlTable$Type
  }
  
  
  
  controlgrn <- log2(1 + controlgrn)
  controlred <- log2(1 + controlred)
  object$ctl.covmat=constructProbCovMat(controlred,controlgrn,cp.types,object$cell_type)
  message("A covariance Matrix was build")
  
  
  ## Genome Studio Data Files
  if(object$type=="genomeStudio"){
    
    ## Building the colClasses to help reading the file
    colClasses <- rep("double", object$sampleSize*2+4)
    colClasses[1:4]=list(NULL)
    signalTable=read.table(object$signalFile,sep='\t',header=TRUE,colClasses=colClasses)
    sigA=signalTable[,1:object$sampleSize*2-1]
    sigB=signalTable[,1:object$sampleSize*2]
    
    rm(signalTable)
    
    #checking sanity of the data
    if (any(!is.finite(as.matrix(sigA)))){stop("There are non-numeric values in the matrix", '\n')}
    if (any(!is.finite(as.matrix(sigB)))){stop("There are non-numeric values in the matrix", '\n')}
    if(any(unlist(strsplit(colnames(sigA),".Signal_A"))!=object$sampleNames)){stop("Those two files should contain the same samples in the same order:\n",controlProbeFile,signalFile)}
    
    #Naming the column consistently
    colnames(sigA)=object$sampleNames
    colnames(sigB)=object$sampleNames
    
    object$nPos=nrow(sigA)
    
    load("data/Order_450K.rda",envir=environment())
    object$signal=list(AIGrn=sigA[orderIGrn,],
                       BIGrn=sigB[orderIGrn,],
                       AIRed=sigA[orderIRed,],
                       BIRed=sigB[orderIRed,],
                       AII=sigA[orderII,],
                       BII=sigB[orderII,])
    message("Signal data loaded")
    rm(sigA);rm(sigB)
    
  } else if (object$type=="minfi") {
    
    
    SnpI <- getProbeInfo(object$rgset, type = "SnpI")
    
    TypeII <- rbind(getProbeInfo(object$rgset, type = "II"),
                    getProbeInfo(object$rgset, type = "SnpII"))
    
    TypeI.Red <- rbind(getProbeInfo(object$rgset, type = "I-Red"),
                       SnpI[SnpI$Color == "Red",])
    
    TypeI.Green <- rbind(getProbeInfo(object$rgset, type = "I-Green"),
                         SnpI[SnpI$Color == "Grn",])
    object$nPos=(nrow(TypeII)+nrow(TypeI.Green)+nrow(TypeI.Red))
    object$signal=list(AIGrn=getGreen(object$rgset)[TypeI.Green$AddressA,],
                       BIGrn=getGreen(object$rgset)[TypeI.Green$AddressB,],
                       AIRed=getRed(object$rgset)[TypeI.Red$AddressA,],
                       BIRed=getRed(object$rgset)[TypeI.Red$AddressB,],
                       AII=getRed(object$rgset)[TypeII$AddressA,],
                       BII=getGreen(object$rgset)[TypeII$AddressA,])
    
    message("Signal data loaded")
    
    object$Annotation=list(IGrn=TypeI.Green$Name,
                           IRed=TypeI.Red$Name,
                           II=TypeII$Name)
    
  }
  
  for(i in names(object$signal)){
    object$signal[[i]]=log2(1 + object$signal[[i]])
  }
  
  object$qntllist=buildQuantileList(object$nPos)
  object$quantiles=list()
  for(i in names(object$signal)){
    object$quantiles[[i]]=colQuantiles(object$signal[[i]],prob=object$qntllist)
  }
  message("Quantiles done")
  
  return(object)
  
}







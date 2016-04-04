################################################################################
## This file allow to create some 'SampleSet' objects, they are list containing
## signal data and different variables usefull for funtooNorm
## The data is separated into the 3 kind positions are quantified
## Each way have then 2 channel (methylated and unmethylated ie : A and B
## We define then the 6 (2*3) labels: AIGrn BIGrn AIRed BIRed AII BII


## This file allow to create some 'SampleSet' objects, they are list containing

#' SampleSet is an S3 class defined for the purpose of running the
#' funtooNorm algorithm. They are lists containing signal data and different
#' variables useful for funtooNorm. The data is separated into the 3 probes
#' types, each having 2 channels (methylated and unmethylated ie : A and B)
#' We then define then the 6 (2*3) labels: AIGrn BIGrn AIRed BIRed AII BII
#' 
#'
#' @slot type character: is 'minfi' or 'GenomeStudio'
#' @slot sampleNames character vector:
#' contain the list of sample names in order used 
#' @slot sampleSize numeric: the number of samples
#' @slot nPos numeric: the number of positions in the ILLUMINA chip
#' @slot annotation IlluminaMethylationAnnotation: the annotation object from
#' mnfi package
#' @slot cell_type list: list matching each sample to define the categories
#' @slot qntllist numeric: vector of ordered quantiles
#' @slot quantiles numeric: list of  6 quantiles tables for 6 type of signals
#' @slot ctl.covmat numeric: covariance matrix for the model fit
#' @slot signal numeric: list of 6 signal tables the 6 type of signals
#'
#'
#' @return a SampleSet object
#' @export
#'
#' @examples showClass("SampleSet")
#' 
setClass("SampleSet", representation(type="character",
                                     sampleNames="character",
                                     sampleSize="numeric",
                                     nPos="numeric",
                                     annotation="IlluminaMethylationAnnotation",
                                     cell_type="list",
                                     qntllist="numeric",
                                     quantiles="numeric",
                                     ctl.covmat="numeric",
                                     signal="numeric")
         )








################################################################################
#' create a SampleSet from RGChannelSet from minfi package
#'
#' @param RGChannelSet, from minfi package, should contain a cell_type vector
#' in it s phenotypes data pData
#'
#' @return a SampleSet object
#' @export
#'
#' @examples require(funtooNorm)
#' require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' 
fromRGChannelSet <- function(myRGChannelSet){
  object <- list(type="minfi")
  class(object) <- "SampleSet"
  
  if(! "cell_type" %in% colnames(minfi::pData(myRGChannelSet))){
    stop("Your object should contain a field \"cell_type\" in it's phenotype
         data minfi::pData(rgset) in order to use funtooNorm")
  }
  
  cell_type <- minfi::pData(myRGChannelSet)$cell_type
  object$sampleSize=length(cell_type)
  object$cell_type=cell_type
  object$sampleNames=rownames(minfi::pData(myRGChannelSet))
  object$annotation=myRGChannelSet@annotation
    

  if (any(cell_type == '' | is.na(cell_type))) {
    stop("There are NA values in cell_type", '\n')
  }
  if (length(unique(cell_type))<2) {
    stop("There should be AT LEAST 2 cell type in cell_type variable", '\n')
  }
  
  ## Formating the control probes data for the covariance matrix
  controlTable <- getProbeInfo(getManifest(myRGChannelSet), type="Control")
  
  #Here we have to remove the 2 probes that are absent from the Chip
  controlTable=controlTable[!controlTable$Address%in%c("21630339","24669308"),]
  #Here we remove the 15 probes that cannot be in the GenomeStudio output
  controlTable=controlTable[controlTable$Color!="-99",]
  controlred=as.data.frame(getRed(myRGChannelSet)[controlTable$Address,])
  
  controlgrn=as.data.frame(getGreen(myRGChannelSet)[controlTable$Address,])
  
  cp.types=controlTable$Type
  
  controlgrn <- log2(1 + controlgrn)
  controlred <- log2(1 + controlred)
  object$ctl.covmat=funtooNorm:::constructProbCovMat(controlred,controlgrn,
                                        cp.types,object$cell_type)
  message("A covariance Matrix was build")

  loc=getLocations(sprintf("%sanno.%s",
                           object$annotation["array"],
                           object$annotation["annotation"]),
                   orderByLocation = FALSE)
  pos=cbind(names(loc),as.character(seqnames(loc)),start(loc))
  
  chrYnames=names(loc)[as.character(seqnames(loc))=="chrY"]
  
  object$signal=list()
  object$names=list()
  
  SnpI <- getProbeInfo(object$annotation, type = "SnpI")
  
  ## Type I Green
  TypeI.Green <- rbind(getProbeInfo(object$annotation, type = "I-Green"),
                       SnpI[SnpI$Color == "Grn",])
  sub=TypeI.Green$Name %in% chrYnames
  object$names$IGrn=TypeI.Green$Name[!sub]
  object$names$chrY=TypeI.Green$Name[sub]
  sigA=getGreen(myRGChannelSet)[TypeI.Green$AddressA,]
  sigB=getGreen(myRGChannelSet)[TypeI.Green$AddressB,]
  object$signal$AIGrn=sigA[!sub,]
  object$signal$BIGrn=sigB[!sub,]
  object$signal$BchrY=sigB[sub,]
  object$signal$AchrY=sigA[sub,]
  
  ## Type I Red
  TypeI.Red <- rbind(getProbeInfo(object$annotation, type = "I-Red"),
                     SnpI[SnpI$Color == "Red",])
  sub=TypeI.Red$Name %in% chrYnames
  object$names$IRed=TypeI.Red$Name[!sub]
  object$names$chrY=c(object$names$chrY,TypeI.Red$Name[sub])
  sigA=getRed(myRGChannelSet)[TypeI.Red$AddressA,]
  sigB=getRed(myRGChannelSet)[TypeI.Red$AddressB,]
  object$signal$AIRed=sigA[!sub,]
  object$signal$BIRed=sigB[!sub,]
  object$signal$AchrY=rbind(object$signal$AchrY,sigA[sub,])
  object$signal$BchrY=rbind(object$signal$BchrY,sigB[sub,])
  
  ## Type II
  TypeII <- rbind(getProbeInfo(object$annotation, type = "II"),
                  getProbeInfo(object$annotation, type = "SnpII"))
  sub=TypeII$Name %in% chrYnames
  object$names$II=TypeII$Name[!sub]
  object$names$chrY=c(object$names$chrY,TypeII$Name[sub])
  sigA=getRed(myRGChannelSet)[TypeII$AddressA,]
  sigB=getGreen(myRGChannelSet)[TypeII$AddressA,]
  object$signal$AII=sigA[!sub,]
  object$signal$BII=sigB[!sub,]
  object$signal$AchrY=rbind(object$signal$AchrY,sigA[sub,])
  object$signal$BchrY=rbind(object$signal$BchrY,sigB[sub,])
  
  object$nPos=sum(length(object$names$chrY),
                  length(object$names$II),
                  length(object$names$IRed),
                  length(object$names$IGrn))
  
  for(i in names(object$signal)){
    object$signal[[i]]=log2(1 + object$signal[[i]])
  }
  
  message("Signal data loaded")
  
  
  object$qntllist=funtooNorm:::buildQuantileList(object$nPos)
  object$quantiles=list()
  for(i in names(object$signal)[!grepl("Y",names(object$signal))]){
    object$quantiles[[i]]=matrixStats::colQuantiles(object$signal[[i]],
                                                    prob=object$qntllist)
  }
  message("Quantiles done")

  return(object)

}


################################################################################
#' create a SampleSet from GenomeStudio files
#'
#' @param controlProbeFile file of control probe data exported from GenomeStudio
#' @param signalFile file exported from GenomeStudio with the exact same samples
#' as control probe File
#' @param cell_type this vector should have names matching all the samples in
#' the files from genome studios, and at least 2 different cell types.
#'
#' @return a SampleSet object
#' @export
#'
fromGenStudFiles <- function(controlProbeFile,signalFile,cell_type){
  object <- list(type="genomeStudio")

  if (any(cell_type == '' | is.na(cell_type))) {
    stop("There are NA values in cell_type", '\n')
  }
  if (length(unique(cell_type))<2) {
    stop("There should be at least 2 cell types\n")
  }
  class(object) <- "SampleSet"
  object$sampleSize=length(cell_type)
  object$cell_type=cell_type
  object$annotation=c(array="IlluminaHumanMethylation450k",
                      annotation="ilmn12.hg19")
  
  
  ## Construct the control Table
  controlTable=utils::read.table(controlProbeFile,sep='\t',header=TRUE)
  object$sampleNames=unlist(strsplit(colnames(controlTable)
                                     [1:object$sampleSize*3+1],".Signal_Grn"))
  
  controlred=controlTable[,paste(object$sampleName,".Signal_Red",sep="")]
  controlgrn=controlTable[,paste(object$sampleName,".Signal_Grn",sep="")]
  cp.types=controlTable$TargetID
  
  if(any(! object$sampleNames %in% names(object$cell_type))){
    i=which(! object$sampleNames%in%names(object$cell_type))[1]
    stop("STOP: ",object$sampleNames[i],"from file ",controlProbeFile,
         " should be present in the cell_types names:\n")
  }
  object$cell_type=object$cell_type[object$sampleNames]
  
  controlgrn <- log2(1 + controlgrn)
  controlred <- log2(1 + controlred)
  object$ctl.covmat=funtooNorm:::constructProbCovMat(controlred,controlgrn,
                                        cp.types,object$cell_type)
  message("A covariance Matrix was build")
  

  ## Building the colClasses to help reading the file
  colClasses <- rep("double", object$sampleSize*2+4)
  colClasses[1:4]=list(NULL)
  colClasses[2]="character"
  signalTable=utils::read.table(signalFile,sep='\t',
                         header=TRUE,colClasses=colClasses)
  sigA=signalTable[,1:object$sampleSize*2]
  sigB=signalTable[,1:object$sampleSize*2+1]
  names=signalTable[,1]
  rm(signalTable)
  
  #checking sanity of the data
  if (any(!is.finite(as.matrix(sigA)))){
    stop("There are non-numeric values in the matrix", '\n')}
  if (any(!is.finite(as.matrix(sigB)))){
    stop("There are non-numeric values in the matrix", '\n')}
  if(any(unlist(strsplit(colnames(sigA),".Signal_A"))!=object$sampleNames)){
    stop("Those two files should contain the same samples in the same order:\n",
         controlProbeFile,signalFile)}
  
  #Naming the column consistently
  colnames(sigA)=object$sampleNames
  colnames(sigB)=object$sampleNames
  
  object$nPos=nrow(sigA)
  
  data("Order_450K")
  
 object$signal=list(AIGrn=sigA[orderIGrn,],
                    BIGrn=sigB[orderIGrn,],
                    AIRed=sigA[orderIRed,],
                    BIRed=sigB[orderIRed,],
                    AII=sigA[orderII,],
                    BII=sigB[orderII,],
                    AchrY=sigA[orderchrY,],
                    BchrY=sigB[orderchrY,])
 
 object$names=list(IGrn=names[orderIGrn],
                   IRed=names[orderIRed],
                   II=names[orderII],
                   chrY=names[orderchrY])
 
 
  message("Signal data loaded")
  rm(sigA);rm(sigB)
  

  for(i in names(object$signal)){
    object$signal[[i]]=log2(1 + object$signal[[i]])
  }
  
  object$qntllist=funtooNorm:::buildQuantileList(object$nPos)
  object$quantiles=list()
  object$quantiles=list()
  for(i in names(object$signal)[!grepl("Y",names(object$signal))]){
    object$quantiles[[i]]=matrixStats::colQuantiles(object$signal[[i]],
                                                    prob=object$qntllist)
  }
  message("Quantiles done")
  
  return(object)
}


################################################################################
#' Print information about the SampleSet
#'
#' @param object  of type SampleSet
#'
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' mySampleSet
#' 
print.SampleSet <- function(object){
  cat("SampleSet object built from ",object$type,'\n')
  cat("Data: ",object$nPos,"positions ")
  cat("and ",object$sampleSize, "samples",'\n')
  cat("   cell type:",levels(object$cell_type),'\n')
  cat("  ",length(object$qntllist),"quantiles",'\n')
  if(is.null(object$predmat)){
    cat("funtooNorm Normalization was not applied",'\n')
  }else{
    cat("funtooNorm Normalization was applied",'\n')
  }
}


getLogSigA <- function(signal){
  return(rbind(signal$AIGrn,
               signal$AIRed,
               signal$AII,
               signal$AchrY))
}

getLogSigB <- function(signal){
  return(rbind(signal$BIGrn,
               signal$BIRed,
               signal$BII,
               signal$BchrY))
}

################################################################################
#' internal function to get the position names returning a vector of
#' position names the preserving the order define by this package
getPositionNames <- function(names){
  return(c(names$IGrn,
           names$IRed,
           names$II,
           names$chrY))
}

    
#' Return a list
#'
#' @param object 
#'
#' @return a GRange object of all the methylated positions
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' getGRanges(mySampleSet)
#' 
getGRanges <- function(object){
  methWithoutSNPs=getPositionNames(object$names)
  methWithoutSNPs=methWithoutSNPs[!grepl("^rs",methWithoutSNPs)]
  loc=getLocations(sprintf("%sanno.%s",
                           object$annotation["array"],
                           object$annotation["annotation"]),
                   orderByLocation = FALSE)
  return(loc[methWithoutSNPs])
}


################################################################################
#' compute the beta value of the raw signal for each position and each
#' sample
#'
#' @param object  of type SampleSet
#' @param offset default is 100 as Illumina standard
#'
#' @return a matrix containing the raw beta value for each position and each
#' samples
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' getRawBeta(mySampleSet)
#' 
getRawBeta <- function(object,offset=100){
  mat=funtooNorm:::calcBeta(getLogSigA(object$signal),
                            getLogSigB(object$signal),
                            offset)
  colnames(mat)=object$sampleNames
  rownames(mat)=getPositionNames(object$names)
  return(mat[!grepl("^rs",rownames(mat)),])
}

################################################################################
#' compute the beta value after normalization for each position and each
#' sample
#'
#' @param object  of type SampleSet
#' @param offset default is 100 as Illumina standard
#'
#' @return a matrix containing beta after normalization value for each CpG 
#' position and each samples
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' b=getNormBeta(funtooNorm(mySampleSet))
#' 
getNormBeta <- function(object,offset=100){
  if(any(is.null(object$predmat))){
    stop("WARNING: please call funtooNorm")
  }
  mat=funtooNorm:::calcBeta(getLogSigA(object$predmat),
                            getLogSigB(object$predmat),
                            offset)
  colnames(mat)=object$sampleNames
  rownames(mat)=getPositionNames(object$names)
  return(mat[!grepl("^rs",rownames(mat)),])
}

################################################################################
#' compute the M value after normalization for each position and each
#' sample
#'
#' @param object  of type SampleSet
#' @param offset default is 100 as Illumina standard
#'
#' @return a matrix containing M after normalization value for each position
#' and each samples log2(Meth/Unmeth)
#' @export
#' 
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' m=getNormM(funtooNorm(mySampleSet))
#' 
getNormM <- function(object,offset=100){
  if(any(is.null(object$predmat))){
    stop("WARNING: please call funtooNorm")
  }
  mat=getLogSigA(object$predmat)-getLogSigB(object$predmat)
  colnames(mat)=object$sampleNames
  rownames(mat)=getPositionNames(object$names)
  return(mat[!grepl("^rs",rownames(mat)),])
}

################################################################################
#' compute the M value after normalization for each SNP position and
#' each sample
#'
#' @param object  of type SampleSet
#' @param offset default is 100 as Illumina standard
#'
#' @return a matrix containing M after normalization value for each SNP of the 
#' chip and each sample log2(Meth/Unmeth)
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' snp=getSnpM(funtooNorm(mySampleSet))
#' 
getSnpM <- function(object){
  if(any(is.null(object$predmat))){
    stop("WARNING: please call funtooNorm")
  }
  mat=getLogSigA(object$predmat)-getLogSigB(object$predmat)
  colnames(mat)=object$sampleNames
  rownames(mat)=getPositionNames(object$names)
  return(mat[grepl("^rs",rownames(mat)),])
}

################################################################################
#' This function applies the normalization method central to the package
#' to each signal.
#' The chrY have a deserve a specific treatment, men are asses using the median
#' beta estimation on the raw data positions with a cutoff at 60%. We perform
#' on the mens a quantile normalization and we do not change the women values
#'
#' @param object of type SampleSet
#' @param type.fits can be "PCR" or "PLS" (default "PCR")
#' @param ncmp number of components used in the analysis (default 4)
#' @param force set it to TRUE in order to re-compute the normalisation whent 
#' it is already done
#' @param sex boolean vector: when not null force the chrY normalization to use 
#' treat the TRUE values as mens
#'
#' @return a SampleSet containing the normalised signal
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' mySampleSet=funtooNorm(mySampleSet)
#' 
funtooNorm <- function(object, type.fits="PCR",ncmp=4,force=FALSE,sex=NULL){
  
  if(force | is.null(object$predmat)){
    object$predmat=list()
  }else{
    message("Normalization already exist")
  }
  
  
  ######## this part deal with chrY
  if(is.null(sex)){
    mens=matrixStats::colMedians(funtooNorm:::calcBeta(object$signal$AchrY,
                                                       object$signal$BchrY))<0.6
    message("we found ",sum(mens)," men and ",
            sum(!mens)," women in your data set base on Y probes only")
  }else{
    mens=sex
    message("There is ",sum(mens)," men and ",
            sum(!mens)," women")
  }
  # women will not have any corrections
  object$predmat$AchrY=object$signal$AchrY 
  object$predmat$BchrY=object$signal$BchrY
  if(1<sum(mens)){
    object$predmat$AchrY[,mens]=
      funtooNorm:::quantileNormalization(object$signal$AchrY[,mens])
    object$predmat$BchrY[,mens]=
      funtooNorm:::quantileNormalization(object$signal$BchrY[,mens])
  }
  
  
  for(signal in names(object$quantiles)){
    if(force | is.null(object$predmat[[signal]])){
      message("Normalization of signal : ",signal)
      object$predmat[[signal]]=
        funtooNorm:::funtooNormApply(object$signal[[signal]],
                                    object$quantiles[[signal]],
                                    object$qntllist,
                                    object$ctl.covmat,
                                    ncmp,type.fits)
      colnames(object$predmat[[signal]])=object$sampleNames
    }else{
      message("Already done : ",signal)
    }
  }
  return(object)
}


################################################################################
#' Plot a series of graphs with different numbers of components for
#' each signal
#'
#' @param object of type SampleSet
#' @param type.fits can be "PCR" or "PLS" (default "PCR")
#' @param file if not empty will write a pdf using this name, path can be
#' included
#'
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' plotValidationGraph(mySampleSet)
#' 
plotValidationGraph <- function(object, type.fits="PCR",file=""){

  if(file=="") message("Will plot in the usual output")
  else if (!file.access(file, mode=2))
    stop("cannot write in ",file)
  else if (tools::file_ext(file)!="pdf")
    stop("your file extension should be pdf for ",file)
    
  svd.ctlcovmat=svd(object$ctl.covmat)$d
  ## set max is 8
  numcomp <- min(c(which(cumsum(svd.ctlcovmat)/sum(svd.ctlcovmat)>=0.98),8))
  if (!is.finite(numcomp) | numcomp>20)  {
    stop('There may be a problem with the SVD decomposition of the control
matrix: number of components needed is either>20 or missing', "\n")
  }
  message("\nStarting validation with a max of ", numcomp , " components...\n")
  
  if(file!=""){
    message("plotting in the pdf file:",file)
    pdf(file = file, height=7, width=7)
  }
  
  layout(mat =  matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE),
         heights = c(0.3,0.3,0.15))
  par(omi=c(0,0,0,0),mar = c(2,2,0.2, 0.2), mgp = c(1,0,0))
  
  
  for(i in names(object$quantiles)){
    funtooNorm:::plotValidate(object$quantiles[[i]],object$qntllist,
                              object$ctl.covmat,numcomp,i)
  }
  
  
  par(mar=c(0, 0, 1, 0)) 
  plot.new()
  ltylist <- rep(1, min(numcomp,8))
  if (numcomp) ltylist <- c(ltylist, rep(2, max(8,numcomp)-8))
  legend(title='Number of components:', x = "top",inset = 0,
         legend = 1:numcomp, col=rainbow(numcomp), lty=ltylist,
         cex=1.1, horiz = TRUE)
  if(file!=""){
    dev.off()
  }
  
}





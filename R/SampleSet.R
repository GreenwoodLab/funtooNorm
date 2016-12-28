################################################################################
## This function creates a S4 object of class 'SampleSet'.
## SampleSet object contain both signal data as well as
## other covariables needed for the funtooNorm function.  
## The data is divided by colour and probe type.
## Each type have 2 channels (methylated and unmethylated ie : A and B
## We define the six groups as (2*3): AIGrn BIGrn AIRed BIRed AII BII

#' @title S4 class object SampleSet
#'
#' @description SampleSet is an S4 class defined for the purpose of running the
#' funtooNorm algorithm. They are lists containing signal data and different
#' variables useful for funtooNorm. The data is separated into the 3 probes
#' types, each having 2 channels (methylated and unmethylated ie : A and B)
#' We then define then the 6 (2*3) labels: AIGrn BIGrn AIRed BIRed AII BII
#' 
#'
#' @slot type Character: is 'minfi' or 'GenomeStudio'
#' @slot sampleNames character vector:
#' contain the list of sample names in order used 
#' @slot sampleSize numeric: the number of samples
#' @slot nPos numeric: the number of positions in the ILLUMINA chip
#' @slot annotation character: the annotation object from
#' minfi package
#' @slot cell_type factor: vector of the cell type for each sample as factors
#' @slot qntllist numeric: vector of ordered quantiles
#' @slot quantiles list: list of  6 quantiles tables for the 6 signal types
#' @slot ctl.covmat matrix: covariance matrix for the model fit
#' @slot signal list: list of the values for all 6 probe types.
#' @slot names list: list of probes for each type
#' @slot predmat list: list of the normalized values for all 6 probe types.
#'
#' @return a S4 object of class SampleSet
#' @export
#'
#' @examples showClass("SampleSet")
#' @importClassesFrom minfi IlluminaMethylationAnnotation
#' 
setClass("SampleSet", representation(type = "character",
                                     sampleNames = "character",
                                     sampleSize = "numeric",
                                     nPos = "numeric",
                                     annotation = "character",
                                     cell_type = "factor",
                                     qntllist = "numeric",
                                     quantiles = "list",
                                     ctl.covmat = "matrix",
                                     signal = "list",
                                     names = "list",
                                     predmat = "list")
         )








################################################################################

#' @title  Creates an object of class SampleSet from a RGChannelSet {minfi}
#' 
#' @description Creates a object of class SampleSet from the raw unprocessed
#' data in RGChannelSet
#' 
#' @param myRGChannelSet : RGChannelSet, from minfi package, should contain a
#'   cell_type vector in pData
#'   
#' @return An object of class 'SampleSet'
#' @export
#' 
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' 
#' @importClassesFrom minfi RGChannelSet
fromRGChannelSet <- function(myRGChannelSet){
    object = new("SampleSet")
    object@type = "minfi"

    if(! "cell_type" %in% colnames(minfi::pData(myRGChannelSet))){
    stop("Your object should contain a field \"cell_type\" in it's phenotype
         data minfi::pData(rgset) in order to use funtooNorm")
    }
    
    cell_type <- minfi::pData(myRGChannelSet)$cell_type
    object@sampleSize=length(cell_type)
    object@cell_type=as.factor(cell_type)
    object@sampleNames=rownames(minfi::pData(myRGChannelSet))
    object@annotation=myRGChannelSet@annotation
    
    
    if (any(cell_type == '' | is.na(cell_type))) {
    stop("There are NA values in cell_type")
    }
    if (length(unique(cell_type))<2) {
    stop("There should be AT LEAST 2 cell types in cell_type variable")
    }
    
    ## Formating the control probes data for the covariance matrix
    controlTable <- minfi::getProbeInfo(minfi::getManifest(myRGChannelSet), 
                                        type="Control")
    
    #Here we have to remove the 2 probes that are absent from the Chip
    controlTable=controlTable[!controlTable$Address%in%c("21630339","24669308"),
                              ]
    #Here we remove the 15 probes that cannot be in the GenomeStudio output
    controlTable=controlTable[controlTable$Color!="-99",]
    controlred=as.data.frame(minfi::getRed(myRGChannelSet)[controlTable$Address,
                                                           ])
    
    controlgrn=as.data.frame(minfi::getGreen(myRGChannelSet)[controlTable$Address,
                                                             ])
    
    cp.types=controlTable$Type
    
    controlgrn <- log2(1 + controlgrn)
    controlred <- log2(1 + controlred)
    object@ctl.covmat=constructProbCovMat(controlred,controlgrn,
                                        cp.types,object@cell_type)
    message("A covariance Matrix was build")
    
    loc=minfi::getLocations(sprintf("%sanno.%s",
                           object@annotation["array"],
                           object@annotation["annotation"]),
                   orderByLocation = FALSE)
    pos=cbind(names(loc),as.character(GenomeInfoDb::seqnames(loc)),start(loc))
    
    chrYnames=names(loc)[as.character(GenomeInfoDb::seqnames(loc))=="chrY"]
    
    SnpI <- minfi::getProbeInfo(object@annotation, type = "SnpI")
    
    ## Type I Green
    TypeI.Green <- rbind(minfi::getProbeInfo(object@annotation, 
                                             type = "I-Green"),
                       SnpI[SnpI$Color == "Grn",])
    sub=TypeI.Green$Name %in% chrYnames
    object@names$IGrn=TypeI.Green$Name[!sub]
    object@names$chrY=TypeI.Green$Name[sub]
    sigA=minfi::getGreen(myRGChannelSet)[TypeI.Green$AddressA,]
    sigB=minfi::getGreen(myRGChannelSet)[TypeI.Green$AddressB,]
    object@signal$AIGrn=sigA[!sub,]
    object@signal$BIGrn=sigB[!sub,]
    object@signal$BchrY=sigB[sub,]
    object@signal$AchrY=sigA[sub,]
    
    ## Type I Red
    TypeI.Red <- rbind(minfi::getProbeInfo(object@annotation, type = "I-Red"),
                     SnpI[SnpI$Color == "Red",])
    sub=TypeI.Red$Name %in% chrYnames
    object@names$IRed=TypeI.Red$Name[!sub]
    object@names$chrY=c(object@names$chrY,TypeI.Red$Name[sub])
    sigA=minfi::getRed(myRGChannelSet)[TypeI.Red$AddressA,]
    sigB=minfi::getRed(myRGChannelSet)[TypeI.Red$AddressB,]
    object@signal$AIRed=sigA[!sub,]
    object@signal$BIRed=sigB[!sub,]
    object@signal$AchrY=rbind(object@signal$AchrY,sigA[sub,])
    object@signal$BchrY=rbind(object@signal$BchrY,sigB[sub,])
    
    ## Type II
    TypeII <- rbind(minfi::getProbeInfo(object@annotation, type = "II"),
                  minfi::getProbeInfo(object@annotation, type = "SnpII"))
    sub=TypeII$Name %in% chrYnames
    object@names$II=TypeII$Name[!sub]
    object@names$chrY=c(object@names$chrY,TypeII$Name[sub])
    sigA=minfi::getRed(myRGChannelSet)[TypeII$AddressA,]
    sigB=minfi::getGreen(myRGChannelSet)[TypeII$AddressA,]
    object@signal$AII=sigA[!sub,]
    object@signal$BII=sigB[!sub,]
    object@signal$AchrY=rbind(object@signal$AchrY,sigA[sub,])
    object@signal$BchrY=rbind(object@signal$BchrY,sigB[sub,])
    
    object@nPos=sum(length(object@names$chrY),
                  length(object@names$II),
                  length(object@names$IRed),
                  length(object@names$IGrn))
    
    for(i in names(object@signal)){
    object@signal[[i]]=log2(1 + object@signal[[i]])
    }
    
    message("Signal data loaded")
    
    
    object@qntllist=buildQuantileList(object@nPos)
    for(i in names(object@signal)[!grepl("Y",names(object@signal))]){
    object@quantiles[[i]]=matrixStats::colQuantiles(object@signal[[i]],
                                                    prob=object@qntllist)
    }
    message("Quantiles done")
    
    return(object)
}


################################################################################
#' Creates a S4 object of class 'SampleSet' from GenomeStudio files
#'
#' @param controlProbeFile The control probe file exported from GenomeStudio
#' @param signalFile The signals exported from GenomeStudio samples must be in
#' same order as the control probe File
#' @param cell_type A vector of cell types, names must match control probes and
#'  signal files.
#'
#' @return An object of class 'SampleSet'.
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
    object <- new("SampleSet")
    object@sampleSize=length(cell_type)
    object@cell_type=cell_type
    object@annotation=c(array="IlluminaHumanMethylation450k",
                      annotation="ilmn12.hg19")
    
    
    ## Construct the control Table
    controlTable=utils::read.table(controlProbeFile,sep='\t',header=TRUE)
    object@sampleNames=unlist(strsplit(colnames(controlTable)
                                     [1:object@sampleSize*3+1],".Signal_Grn"))
    
    controlred=controlTable[,paste(object@sampleName,".Signal_Red",sep="")]
    controlgrn=controlTable[,paste(object@sampleName,".Signal_Grn",sep="")]
    cp.types=controlTable$TargetID
    
    if(any(! object@sampleNames %in% names(object@cell_type))){
    i=which(! object@sampleNames%in%names(object@cell_type))[1]
    stop("STOP: ",object@sampleNames[i],"from file ",controlProbeFile,
         " should be present in the cell_types names:\n")
    }
    object@cell_type=object@cell_type[object@sampleNames]
    
    controlgrn <- log2(1 + controlgrn)
    controlred <- log2(1 + controlred)
    object@ctl.covmat=constructProbCovMat(controlred,controlgrn,
                                        cp.types,object@cell_type)
    message("A covariance Matrix was build")
    
    
    ## Building the colClasses to help reading the file
    colClasses <- rep("double", object@sampleSize*2+4)
    colClasses[1:4]=list(NULL)
    colClasses[2]="character"
    signalTable=utils::read.table(signalFile,sep='\t',
                         header=TRUE,colClasses=colClasses)
    sigA=signalTable[,1:object@sampleSize*2]
    sigB=signalTable[,1:object@sampleSize*2+1]
    names=signalTable[,1]
    rm(signalTable)
    
    #checking sanity of the data
    if (any(!is.finite(as.matrix(sigA)))){
    stop("There are non-numeric values in the matrix", '\n')}
    if (any(!is.finite(as.matrix(sigB)))){
    stop("There are non-numeric values in the matrix", '\n')}
    if(any(unlist(strsplit(colnames(sigA),".Signal_A"))!=object@sampleNames)){
    stop("Those two files should contain the same samples in the same order:\n",
         controlProbeFile,signalFile)}
    
    #Naming the column consistently
    colnames(sigA)=object@sampleNames
    colnames(sigB)=object@sampleNames
    
    object@nPos=nrow(sigA)
    
    object@signal=list(AIGrn=sigA[orderIGrn,],
                    BIGrn=sigB[orderIGrn,],
                    AIRed=sigA[orderIRed,],
                    BIRed=sigB[orderIRed,],
                    AII=sigA[orderII,],
                    BII=sigB[orderII,],
                    AchrY=sigA[orderchrY,],
                    BchrY=sigB[orderchrY,])
    
    object@names=list(IGrn=names[orderIGrn],
                   IRed=names[orderIRed],
                   II=names[orderII],
                   chrY=names[orderchrY])
    
    
    message("Signal data loaded")
    
    for(i in names(object@signal)){
    object@signal[[i]]=log2(1 + object@signal[[i]])
    }
    
    object@qntllist=buildQuantileList(object@nPos)
    object@quantiles=list()
    object@quantiles=list()
    for(i in names(object@signal)[!grepl("Y",names(object@signal))]){
    object@quantiles[[i]]=matrixStats::colQuantiles(object@signal[[i]],
                                                    prob=object@qntllist)
    }
    message("Quantiles done")
    
    return(object)
}


################################################################################
#' Show Object SampleSet
#'
#' @description Display informations about the SampleSet object
#' @param object an object of class SampleSet
#' @param ... optional arguments passed to or from other methods.
#'
#' @return No value is returned. The function prints the summary of object of 
#' class SampleSet to screen
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' mySampleSet
#' 

setMethod(f="show",
          signature = "SampleSet",
          definition = function(object){
              cat("SampleSet object built from ",object@type,'\n')
              cat("Data: ",object@nPos,"positions ")
              cat("and ",object@sampleSize, "samples",'\n')
              cat("   cell type:",levels(object@cell_type),'\n')
              cat("  ",length(object@qntllist),"quantiles",'\n')
              if(length(object@predmat)==0){
                  cat("funtooNorm Normalization was not applied",'\n')}else{
                      cat("funtooNorm Normalization was applied",'\n')
                  }
              }
)

setGeneric(name="getLogSigA",
           def=function(object) standardGeneric("getLogSigA")
)
setMethod("getLogSigA",
          signature = "SampleSet",
          definition = function(object){
              return(rbind(object@signal$AIGrn,
                           object@signal$AIRed,
                           object@signal$AII,
                           object@signal$AchrY))
          }
)

setGeneric(name="getLogSigB",
           def=function(object) standardGeneric("getLogSigB")
)
setMethod("getLogSigB",
          signature = "SampleSet",
          definition = function(object){
              return(rbind(object@signal$BIGrn,
                           object@signal$BIRed,
                           object@signal$BII,
                           object@signal$BchrY))
          }
)

################################################################################
## internal function to get the position names returning a vector of
## position names the preserving the order define by this package
setGeneric(name="getPositionNames",
           def=function(object) standardGeneric("getPositionNames")
)

setMethod("getPositionNames",
          signature =  "SampleSet",
          definition = function(object){
              return(c(object@names$IGrn,
                       object@names$IRed,
                       object@names$II,
                       object@names$chrY))
          }
)

    
#' Build GRange object of methylation probes
#'
#' @param object Object of class SampleSet.
#'
#' @return A GRange object of the positions of each cpg.
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' gr=getGRanges(mySampleSet)
#' 
setGeneric(name="getGRanges",
           def=function(object) standardGeneric("getGRanges")
)
#' @describeIn getGRanges Build GRange object of methylation probes
setMethod("getGRanges",
          signature = "SampleSet",
          definition = function(object){
              methWithoutSNPs=getPositionNames(object)
              methWithoutSNPs=methWithoutSNPs[!grepl("^rs",methWithoutSNPs)]
              loc=minfi::getLocations(sprintf("%sanno.%s",
                                              object@annotation["array"],
                                              object@annotation["annotation"]),
                                      orderByLocation = FALSE)
    return(loc[methWithoutSNPs])
          }
)


################################################################################
#' Computes Beta value from raw signals
#'
#' @param object object of class SampleSet
#' @param offset default is 100 as Illumina standard
#'
#' @return a matrix containing the raw beta value for each position and each
#' samples
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' r=getRawBeta(mySampleSet)
#' 
setGeneric(name="getRawBeta",
           def=function(object, offset=100) standardGeneric("getRawBeta")
)
#' @describeIn getRawBeta Computes Beta value from raw signals
setMethod("getRawBeta",
          signature = "SampleSet",
          definition = function(object,offset){
              mat=calcBeta(getLogSigA(object),
                           getLogSigB(object),
                           offset)
              colnames(mat)=object@sampleNames
              rownames(mat)=getPositionNames(object)
              return(mat[!grepl("^rs",rownames(mat)),])
          }
)

################################################################################
#' Computes Beta values from normalized signals
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
setGeneric(name="getNormBeta",
           def=function(object, offset=100) standardGeneric("getNormBeta")
)
#' @describeIn getNormBeta Computes Beta values from normalized signals
setMethod("getNormBeta",
          signature = "SampleSet",
          definition = function(object,offset){
              if(length(object@predmat)==0){
                  stop("WARNING: please call funtooNorm")
                  }
              mat=calcBeta(getLogSigA(object),
                           getLogSigB(object),
                           offset)
              colnames(mat)=object@sampleNames
              rownames(mat)=getPositionNames(object)
              return(mat[!grepl("^rs",rownames(mat)),])
          }
)

################################################################################
#' Computes M values,log2(Meth/Unmeth), from normalized signals
#'
#' @param object  An object of class SampleSet
#'
#' @return a matrix containing M values, log2(Meth/Unmeth), after normalization
#' @export
#' 
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' m=getNormM(funtooNorm(mySampleSet))
#'
setGeneric(name="getNormM",
           def=function(object) standardGeneric("getNormM")
)
#' @describeIn getNormM Computes M values, log2(Meth/Unmeth),
#'  from normalized signals 
setMethod("getNormM",
          signature = "SampleSet",
          definition = function(object){
              if(length(object@predmat)==0){
                  stop("WARNING: please call funtooNorm")
                  }
              mat=getLogSigB(object)-getLogSigA(object)
              colnames(mat)=object@sampleNames
              rownames(mat)=getPositionNames(object)
              return(mat[!grepl("^rs",rownames(mat)),])
          }
)

################################################################################
#' Computes M values after normalization of SNP data.
#'
#' @param object  of class SampleSet
#'
#' @return a matrix containing M values, log2(Meth/Unmeth), after normalization
#' for SNP data
#'
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' snp=getSnpM(funtooNorm(mySampleSet))
#' 
setGeneric(name="getSnpM",
           def=function(object) standardGeneric("getSnpM")
)
#' @describeIn getSnpM Computes M values, log2(Meth/Unmeth), for normalized
#' SNP data
setMethod("getSnpM",
          signature = "SampleSet",
          definition = function(object){
              if(length(object@predmat)==0){
                  stop("WARNING: please call funtooNorm")
                  }
              mat=getLogSigA(object)-getLogSigB(object)
              colnames(mat)=object@sampleNames
              rownames(mat)=getPositionNames(object)
              return(mat[grepl("^rs",rownames(mat)),])
          }
)

################################################################################
#' @title The funtooNorm normalization function
#' 
#' @description 
#' \code{funtooNorm} Returns the normalized signals to the SampleSet object
#' 
#' @details
#' This is a generic function which applies to autosomes and the X 
#' chromosome. Chromosome Y requires separate analysis as there are few probes 
#' on Y.  We use a straightforward quantile normalization applied to males only.
#'
#' @param object Object of class SampleSet
#' @param type.fits Choice between "PCR" or "PLS" (default="PCR")
#' @param ncmp Number of components included in the analysis (default=4)
#' @param force If set to TRUE, forces the normalization procedure to re-compute
#' @param sex Boolean vector if male. if NULL Beta values from ChrY are used for
#'  classification.
#'
#' @return a S4 object of class SampleSet containing the normalized signal
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' mySampleSet=funtooNorm(mySampleSet)
#'
setGeneric(name="funtooNorm",
           def=function(object, type.fits="PCR",ncmp=4,force=FALSE,sex=NULL) 
               standardGeneric("funtooNorm")
)
#' @describeIn funtooNorm The funtooNorm normalization function
setMethod("funtooNorm",
          signature = "SampleSet",
          definition = function(object,type.fits,ncmp,force,sex)
              {if(length(object@predmat)!=0)
              {message("Normalization already exist")
                  }
              ###### this part deal with chrY
              if(is.null(sex)){
                  mens=matrixStats::colMedians(calcBeta(object@signal$AchrY,
                                                    object@signal$BchrY))<0.6
                  message("we found ",sum(mens)," men and ",sum(!mens),
                          " women in your data set base on Y probes only")
                  }else{
                      mens=sex
                      message("There is ",sum(mens)," men and ",
                              sum(!mens)," women")
                      }
              # no correction for women
              object@predmat$AchrY=object@signal$AchrY
              object@predmat$BchrY=object@signal$BchrY
              if(1<sum(mens)){
                  object@predmat$AchrY[,mens]=
                      quantileNormalization(object@signal$AchrY[,mens])
                  object@predmat$BchrY[,mens]=
                      quantileNormalization(object@signal$BchrY[,mens])
                  }
              for(signal in names(object@quantiles)){
                  if(force | is.null(object@predmat[[signal]])){
                      message("Normalization of signal : ",signal)
                      object@predmat[[signal]]=
                          funtooNormApply(object@signal[[signal]],
                                          object@quantiles[[signal]],
                                          object@qntllist,
                                          object@ctl.covmat,
                                          ncmp,type.fits)
                      colnames(object@predmat[[signal]])=object@sampleNames
                      }else{
                          message("Already done : ",signal)
                      }
                  }
              return(object)
          }
)


################################################################################
#' @title plot of Validation Graph for determing number of components
#' @description Plots a series of graphs for each signal type, to determine 
#' the number of components to include in the normalization procedure.
#'
#' @param object of class SampleSet
#' @param type.fits can be "PCR" or "PLS" (default "PCR")
#' @param pdf.file if no file name is provided print pdf file 
#' plotValidationGraph.pdf in working directory.
#'
#' @return No value is returned.  The function prints the plots to a pdf file.
#' @export
#'
#' @examples require(minfiData)
#' pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
#' mySampleSet=fromRGChannelSet(RGsetEx)
#' plotValidationGraph(mySampleSet)
#' 
setGeneric(name="plotValidationGraph",
           def=function(object, type.fits="PCR",pdf.file=NULL) 
               standardGeneric("plotValidationGraph")
)
#' @describeIn plotValidationGraph Plots a series of graphs for each 
#' signal type, to determine the number of components to include 
#' in the normalization procedure.
setMethod("plotValidationGraph",
          signature = "SampleSet",
          definition = function(object,type.fits,pdf.file){
    
    if(is.null(pdf.file)){pdf.file=file.path(getwd(),"plotValidationGraph.pdf");
    message("PDF file was not provided, will print plotValidationGraph.pdf
            to working directory")}else {if(file.access(dirname(pdf.file), mode=2)==-1){
    stop("cannot write in ",dirname(pdf.file))}else {if (tools::file_ext(pdf.file)!="pdf"){
    stop("your file extension should be .pdf for ",file)} 
    } 
                }
        
    
    svd.ctlcovmat=svd(object@ctl.covmat)$d
    ## set max is 8
    numcomp <- min(c(which(cumsum(svd.ctlcovmat)/sum(svd.ctlcovmat)>=0.98),8))
    if (!is.finite(numcomp) | numcomp>20)  {
    stop('There may be a problem with the SVD decomposition of the control
    matrix: number of components needed is either>20 or missing')
    }
    message("\nStarting validation with a max of ", numcomp , 
            " components...\n")
    
    message("writing to file:",pdf.file)
    pdf(file = pdf.file, height=7, width=7)
    layout(mat =  matrix(c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, 
                         byrow = TRUE),
           heights = c(0.3,0.3,0.15))
    par(omi=c(0,0,0,0),mar = c(2,2,0.2, 0.2), mgp = c(1,0,0))
    
    
    for(i in names(object@quantiles)){
    plotValidate(object@quantiles[[i]],object@qntllist,
                 object@ctl.covmat,numcomp,i)
    }
    
    
    par(mar=c(0, 0, 1, 0)) 
    plot.new()
    ltylist <- rep(1, min(numcomp,8))
    if (numcomp) ltylist <- c(ltylist, rep(2, max(8,numcomp)-8))
    legend(title='Number of components:', x = "top",inset = 0,
         legend = 1:numcomp, col=rainbow(numcomp), lty=ltylist,
         cex=1.1, horiz = TRUE)
    dev.off()
    
          }

    )

#' @import IlluminaHumanMethylation450kmanifest 
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import methods
#' @importFrom stats predict start
#' @importFrom graphics matplot text layout legend par plot.new
#' @importFrom grDevices rainbow dev.off pdf
NULL



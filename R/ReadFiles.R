
## Here we define a way to read manifest files and datasets
## chipType can be 450K or 850K
## Most of the complexity is handdled in the way the manifest file is organised

# Function to read manifest from ILLUMINA

readManifestFile <- function(chipType = "450K",manifestFile= NULL,minimal=TRUE) {
  if(is.null(manifestFile)){
    if(chipType=="850K"){
      manifestFile="../../ILLUMINA_DOC/comparingmanifests/MethylationEPIC_v-1-0_B1.csv"
    } else if (chipType=="450K"){
      manifestFile="../../ILLUMINA_DOC/comparingmanifests/HumanMethylation450_15017482_v1-2.csv"
    } else {
        stop(paste(chipType," is not a recognize manifest type, use 450K or 850K"))
    }
  }
  
  if(!file.exists(manifestFile)){
      stop("The manifest file do not exist:", manifestFile)
  }
    
  message(paste("Reading manifest file for",chipType,"array"))
  
  # Line number of the end of the classic positions list
  control.line <- system(sprintf("grep -n \\\\[Controls\\\\] %s", manifestFile), intern = TRUE)
  control.line <- as.integer(sub(":.*", "", control.line))
  
  # Line number of the head of the classic probes table
  assay.line <- system(sprintf("grep -n \\\\[Assay\\\\] %s", manifestFile), intern = TRUE)
  assay.line <- as.integer(sub(":.*", "", assay.line))
  
  # Sanity check
  if(length(control.line) != 1 || !is.integer(control.line) || is.na(control.line)
     || length(assay.line) != 1 || !is.integer(assay.line) || is.na(assay.line)){
    stop("The manifest file seems corrupted")
  }
  
  
  colNames <- readLines(manifestFile, n = assay.line+1L)[assay.line+1L]
  colNames <- strsplit(colNames, ",")[[1]]
  colClasses <- rep("character", length(colNames))
  names(colClasses) <- colNames
  names(colClasses) <- make.names(names(colClasses))
  
  if(minimal){ # To reduce the memory taken by the manifest
    colClasses[names(colClasses)] <- list(NULL)
    colClasses[c("Name","Infinium_Design_Type","Color_Channel")] <- "character"
  }else{
    colClasses[c("MAPINFO")] <- "integer"
  }
  colClasses[c("AddressA_ID","AddressB_ID")] <- "integer"
  
  # Reading the first part of the 
  manifest <- read.table(manifestFile, header = FALSE, col.names = names(colClasses),
                         sep = ",", comment.char = "", quote = "",
                         skip = assay.line + 1L, colClasses = colClasses,
                         nrows = control.line - assay.line - 2L)
  
  snpPositions=grep("^rs", manifest$Name)
  SNP=manifest[snpPositions,]
  CpG=manifest[-snpPositions,]
  
  byAddress <- function(tm){
    typeI=(tm$Infinium_Design_Type=="I")
    nTypeI=length(typeI)
    tm[typeI,]$Infinium_Design_Type=paste(tm[typeI,]$Infinium_Design_Type,tm[typeI,]$Color_Channel,sep="-")
    tm$Color_Channel=NULL
    colnames(tm)=replace(colnames(tm),which(colnames(tm)=="AddressA_ID"),"Address")
    
    duplicatedTypeI=tm[typeI,]
    duplicatedTypeI$Address=duplicatedTypeI$AddressB_ID
    duplicatedTypeI$AddressB_ID=NULL
    duplicatedTypeI$Infinium_Design_Type=paste(tm[typeI,]$Infinium_Design_Type,"B",sep="-")
    tm$AddressB_ID=NULL
    tm[typeI,]$Infinium_Design_Type=paste(tm[typeI,]$Infinium_Design_Type,"A",sep="-")
    orderFinal=order(c(tm$Address,duplicatedTypeI$Address))
    orderInitial=order(orderFinal)
    
    # To keep track of those pairs of probes that make a position
    listOfTypeIPairs=cbind(orderInitial[which(typeI)],orderInitial[nrow(tm)+1:sum(typeI)])
    colnames(listOfTypeIPairs)=c("probeA","probeB")
    
    list(table=rbind(tm,duplicatedTypeI)[orderFinal,],
         typeI=listOfTypeIPairs
        )
  }
  
  SNP=byAddress(SNP)
  CpG=byAddress(CpG)
  
  controls <- read.table(manifestFile, skip = control.line,
                         sep = ",", comment.char = "", quote = "",
                         colClasses = rep("character", 4))[,1:4]
  
  names(controls) <- c("Address", "Type", "Color", "ExtendedType")
  controls$Address=as.integer(controls$Address)
  
  if(chipType=="450K"){
    controls=controls[controls$Address!=21630339 & controls$Address!=24669308,]
  }
  controls[order(controls$Address),]
  
  list(chipType=chipType,CpG=CpG,SNP=SNP,controls=controls)
}

manifest450K=readManifestFile("450K")
manifest850K=readManifestFile("850K")














require("illuminaio")


readRedGreen <- function(basename,chipType,path="") {
  if(chipType=="850K"){
    manifest=manifest850K
  } else if (chipType=="450K"){
    manifest=manifest450K
  } else {
    stop(paste(chipType," is not a recognize manifest type, use 450K or 850K"))
  }
    
  basename <- sub("_Grn\\.idat$", "", basename)
  basename <- sub("_Red\\.idat$", "", basename)
  grnFile <- paste(basename, "_Grn.idat", sep = "")
  redFile <- paste(basename, "_Red.idat", sep = "")
  if(!all(file.exists(c(redFile,grnFile)))) {
    noexits <- G.files[!file.exists(c(redFile,grnFile))]
    stop("The following specified files do not exist:", paste(noexits, collapse = ", "))
  }
  
  #grnMean=readIDAT(grnFile)$Quants[, "Mean"]
  #redMean=readIDAT(redFile)$Quants[, "Mean"]
  
  list(chipType=chipType,
       grnMean=readIDAT(grnFile)$Quants[, "Mean"],
       redMean=readIDAT(redFile)$Quants[, "Mean"])
       #grnCpG = as.numeric( grnMean[as.character(manifest$CpG$Address)]),
       #grnCtrl= as.numeric( grnMean[as.character(manifest$controls$Address)]),
       #grnSNP = as.numeric( grnMean[as.character(manifest$SNP$Address)]),
       #redCpG = as.numeric( redMean[as.character(manifest$CpG$Address)]),
       #redCtrl= as.numeric( redMean[as.character(manifest$controls$Addrsefsa)]),
       #redSNP = as.numeric( redMean[as.character(manifest$SNP$Address)]),
       #typeIGrnBList=sapply(which(manifest$CpG$Infinium_Design_Type=="I-Grn-A"),function(x) which(manifest$CpG$Infinium_Design_Type=="I-Grn-B" & manifest$CpG$Name==manifest$CpG$Name[x])),
       #typeIRedBList=sapply(which(manifest$CpG$Infinium_Design_Type=="I-Red-A"),function(x) which(manifest$CpG$Infinium_Design_Type=="I-Red-B" & manifest$CpG$Name==manifest$CpG$Name[x])))
}




getTypeI <- function(data){
  if(data$chipType=="850K"){
    manifest=manifest850K
  } else if (data$chipType=="450K"){
    manifest=manifest450K
  }
  manifest$CpG$table[unlist(t(data.frame(manifest$CpG$typeI)),use.names = FALSE),]
}



getTypeII <- function(data){
  if(data$chipType=="850K"){
    manifest=manifest850K
  } else if (data$chipType=="450K"){
    manifest=manifest450K
  }
  manifest$CpG$table[which(manifest$CpG$table$Infinium_Design_Type=="II"),]
}
  
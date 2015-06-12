#use this to mount remote data locally on your linux machine
#system('sshfs aaa@bbb:/home/data1/share ~/mnt/grinux')
# to commit in github
#git remote add origin https://github.com/USERNAME/PROJECT.git
#git push -u origin master
library(pls)
########################################################
PATH_C=c("~/mnt/grinux/greenwood.group/PROJECTS/Marie_Hudson2/Methylation/09-08-2014/",
         "~/mnt/grinux/greenwood.group/PROJECTS/Ludmer_Van_Ijsendoorn/methylation/24-12-2014/",
         "~/mnt/grinux/greenwood.group/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/placenta_bloodcord/")
i=3
PATH=PATH_C[i]

library(pls)
#source("~/data/methylation/SCRIPT/funtoonormfunctions2.V2.R")

#####################
#Annotation file
#####################
annotation<-read.delim(paste(PATH,"RAW_DATA/Annotation_Illumina450k.txt",sep=""),as.is=T,sep="\t")

load(paste(PATH,"RAW_DATA/RAWSignalA.RData", sep=""))  # latest version Feb  2015  
rownames(datSigA) <- datSigA[,1]
cpg=intersect(rownames(datSigA),annotation$TargetID)
sigA <- datSigA[cpg,2:ncol(datSigA)]
load(paste(PATH,"RAW_DATA/RAWSignalB.RData", sep=""))  # latest version Feb  2015  
rownames(datSigB) <- datSigB[,1]
sigB <- datSigB[cpg,2:ncol(datSigB)]


# reduce annotation information to the same probe set and reorder
Annot<-data.frame(Index=seq(1,nrow(sigA),by=1),"probe"=annotation$TargetID, Type=annotation$INFINIUM_DESIGN_TYPE,
                  Color=annotation$COLOR_CHANNEL,Build=annotation$GENOME_BUILD, Chr=annotation$CHR,
                  Mapinfo=annotation$MAPINFO, stringsAsFactors=F)
rownames(Annot)=Annot$probe
Annot=Annot[cpg,]
# OPTIONs
########### !!!!!!!!!!!!!!!!
# number of PLS components
ncmp <- 4

# cleanup to save space
rm(datSigA)
rm(annotation)
rm(datSigB)


########################################################
#ClinicalData
SampleInfo<-read.csv(paste(PATH,"COVARIABLE_DATA/CovariableModified.csv",sep=""),header=T,  as.is=T) 
ID=rownames(SampleInfo)=paste(SampleInfo$Sentrix_ID, SampleInfo$Sentrix_Position, sep="_")
#only keep samples that have cell type information
cell_type=SampleInfo$tissue #change column name here
keep=rownames(SampleInfo)[!is.na(cell_type)]
SampleInfo=SampleInfo[keep,] 

colnames(sigA)= substr(colnames(sigA),2,(nchar(ID[1])+1))
colnames(sigB)= substr(colnames(sigB),2,(nchar(ID[1])+1))
cov.match <- intersect(colnames(sigA),keep)
SampleInfo2 <- SampleInfo[cov.match,]; sigA=sigA[ ,cov.match];sigB=sigB[,cov.match]
cell_type=SampleInfo2$tissue #change column name here

#check if order is matching
identical(colnames(sigA),paste(SampleInfo2$Sentrix_ID, SampleInfo2$Sentrix_Position, sep="_"))
identical(row.names(sigA), Annot$probe)

#################################################
# Control Probes
#################################################
CONTROL<-read.delim(paste(PATH,"/RAW_DATA/ControlProbeProfile.txt",sep=""),as.is=T, sep="\t")

# now match control data to sample data
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# this needs to be revised with care if control data are available for all sample
# apply  a log transformation to all signals  (this decision could be altered)
matgrn <- CONTROL[,seq(from=4, to=ncol(CONTROL), by=3)]
matred <- CONTROL[,seq(from=5, to=ncol(CONTROL), by=3)]
colnames(matgrn)= substr(colnames(matgrn),2,(nchar(ID[1])+1))
colnames(matred)= substr(colnames(matred),2,(nchar(ID[1])+1))
matgrn=matgrn[,cov.match]
matred=matred[,cov.match]
SampleInfo=SampleInfo[cov.match,]
cp.types <- CONTROL[,2]

# check that ID order matches between sigA, sigB and control probe data
identical(colnames(matgrn), colnames(sigA))
identical(colnames(matgrn), rownames(SampleInfo))

##### normalize data #
# this creates the fits required to do the normalization
sigAsample<-sigA[1:20000,]; sigBsample=sigB[1:20000,]; Annotsample<-Annot[1:20000,] 
save(sigAsample, sigBsample, Annotsample, matred, matgrn, cp.types, cell_type, file='./data/data.rda')
Annot<-default.Annot
devtools::use_data(Annot, internal = F, overwrite = T)
devtools::use_data(cp.types, internal = F, overwrite = T)

funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample, Annot=Annotsample, 
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type,
                       ncmp=4, save.quant=TRUE, save.loess=TRUE, apply.loess=TRUE, logit.quant=TRUE, validate=F)


funtoonormout <- funtoonorm(sigA=sigAsample, sigB=sigBsample,
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type,
                      ncmp=4, save.quant=TRUE, save.loess=TRUE, apply.loess=FALSE, logit.quant=FALSE, validate=5)
#convert -quality 100 -density 150 -sharpen 0x1.0 validationCurves.pdf valid.jpg


origBeta <- funtoonormout[[1]]
newBeta <- funtoonormout[[2]]
rownames(newBeta) <- rownames(origBeta)
colnames(newBeta) <- colnames(origBeta)
#save(newBeta, file=paste(PATH,"NEW_NORM_DATA/newBeta4.RData",sep=""))
#save(origBeta, file=paste(PATH,"NEW_NORM_DATA/origBeta.RData",sep=""))                                                                       

####  end of normalization algorithm
################################################
#BETA<-read.delim(paste(PATH,"RAW_DATA/AVG_Beta.txt", sep=""),sep="\t",as.is=T)
#colnames(BETA)=substr(colnames(BETA),2,18); rownames(BETA)=BETA[,1]
#AVG.BETA=BETA[cpg,keep]
#save(AVG.BETA,file=paste(PATH,"RAW_DATA/AVG.BETA.RData",sep="")
     
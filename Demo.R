#vignette 
##### normalize data #
# this creates the fits required to do the normalization
#sigAsample<-sigA[1:1000,]; sigBsample=sigB[1:1000,]; Annotsample<-Annot[1:1000,] 
library(pls)
save(sigAsample, sigBsample, Annotsample, matred, matgrn, cp.types, cell_type, file='data.Rda')


funtoonormout <- funtoonorm(ifquant=FALSE, sigA=sigAsample, sigB=sigBsample, Annot=Annotsample, 
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type,
                      ncmp = ncmp, save.loess=TRUE, applyloess=TRUE, logit.quant=TRUE)
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

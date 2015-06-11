#script to extract sygA, sygB redctrl,grnctrl for NewNorm packege to avoid using minfi
#this fonction read the idat file from path.BEFORE NORMALIZATION.
path='~' ; 
#load library

#LIB_METH = "~/share/greenwood.group/Rlibs/methylationR3.1"

require(minfi)
require(IlluminaHumanMethylation450kmanifest)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
require(matrixStats)
require(shinyMethyl)
require(shinyMethylData)
require(wateRmelon)

#read the sampleSheet and IDAT data
targets          <- read.450k.sheet(paste(path,"IDAT", sep= "/"))
RGset            <- read.450k.exp( targets = targets,  extended = TRUE)

#write the raw data in  Raw data folder
raw_data         <- preprocessRaw(RGset)                 #all the raw methylation data
raw_Bvalue       <- getBeta(raw_data, type = "Illumina") #only the raw Bvalue
raw_Mvalue     <- getM(raw_data)    				   #only the raw Mvalue
raw_snp          <- getSnpBeta(RGset)					   #only the SNP value
detection_pvalue <- detectionP(RGset)                    #produce the detection Pvalue
nbeadcount       <- beadcount(RGset)					   #produce the number of bead per samples		
colnames(nbeadcount)<- gsub("^X","",colnames(nbeadcount))

save(raw_Bvalue       , file = paste(path,"RAW_DATA/raw_Bvalue.RData", sep = "/"))
save(detection_pvalue , file = paste(path,"RAW_DATA/detection_pvalue.RData", sep = "/"))
save(raw_snp          , file = paste(path,"RAW_DATA/raw_snp.RData", sep = "/"))
save(RGset            , file = paste(path,"RAW_DATA/RGset.RData", sep = "/"))
save(nbeadcount       , file = paste(path,"RAW_DATA/nbeadcount.RData", sep = "/"))

#write the Fun_Norm data in normalize data folder
FunNorm_data_bckdyeAdj       <- preprocessFunnorm(RGset,sex = NULL,, verbose=T)           
FunNorm_data_NobckdyeAdj     <- preprocessFunnorm(RGset,sex = NULL, verbose=T, bgCorr = FALSE, dyeCorr = FALSE)          
FunNorm_Bvalue_bckdyeAdj     <- getBeta(FunNorm_data_bckdyeAdj)                         
FunNorm_Bvalue_NobckdyeAdj   <- getBeta(FunNorm_data_NobckdyeAdj)                         

save(FunNorm_data_bckdyeAdj     , file = paste(path,"FUN_NORM_DATA/FunNorm_data_bckdyeAdj.RData", sep = "/"))
save(FunNorm_Bvalue_bckdyeAdj   , file = paste(path,"FUN_NORM_DATA/FunNorm_Bvalue_bckdyeAdj.RData", sep = "/"))
save(FunNorm_data_NobckdyeAdj   , file = paste(path,"FUN_NORM_DATA/FunNorm_data.RData", sep = "/"))
save(FunNorm_Bvalue_NobckdyeAdj , file = paste(path,"FUN_NORM_DATA/FunNorm_Bvalue.RData", sep = "/"))

#write the annotation data in normalize data folder
annotation <- getAnnotation(FunNorm_data_bckdyeAdj)
location   <- getLocations(FunNorm_data_bckdyeAdj)
save(annotation, file = paste(path,"ANNOTATION_DATA/annotation.RData", sep = "/"))
save(location  , file = paste(path,"ANNOTATION_DATA/location.RData", sep = "/"))

#produce a boxplot to compare data no-normalized and normalized
pdf(paste(path,"RAW_DATA/boxplot_no_normData.pdf", sep= "/"))
par(mar = c(10,5,4,1))
boxplot(raw_Bvalue, main = "Boxplot of no_normalized_Bvalues \n", las = 2, ylab = "B_Values")
dev.off()

pdf(paste(path,"FUN_NORM_DATA/boxplot_normData.pdf", sep= "/"))
par(mar = c(10,5,4,1))
boxplot(FunNorm_Bvalue_bckdyeAdj, main = "Boxplot of normalized_Bvalues \n", las = 2, ylab = "B_Values")
dev.off()

pdf(paste(path,"FUN_NORM_DATA/boxplot_RawData_NormData.pdf", sep= "/"))
par(mar = c(10,5,4,1))
boxplot(cbind(raw_Bvalue,FunNorm_Bvalue_bckdyeAdj), main = "Boxplot of Bvalues \n", las = 2, ylab = "B_Values", col = c(rep("red", ncol(raw_Bvalue)), rep("blue",ncol(FunNorm_data_bckdyeAdj))))
legend("topright", inset=.05, c("raw data","functional normalization"), fill = c("red","blue"), horiz=F)
dev.off()

pdf(paste(path,"RAW_DATA/boxplot_RawData_NormData.pdf", sep= "/"))
par(mar = c(10,5,4,1), cex.lab=1)
boxplot(cbind(raw_Bvalue,FunNorm_Bvalue_bckdyeAdj), main = "Boxplot of Bvalues \n", las = 2, ylab = "B_Values", col = c(rep("red", ncol(raw_Bvalue)), rep("blue",ncol(FunNorm_Bvalue_bckdyeAdj))), cex.lab = 0.8)
legend("topright", c("raw data","functional normalization"), fill = c("red","blue"), horiz=F)
dev.off()

setwd("~")
#end of processionIdatData

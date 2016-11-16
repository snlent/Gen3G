###########################################################
##
##
##   850K: Gen3G Data Processing
##
##
##
##
##
##

## Check Version of R
R.version$version.string  

## BioC packages needed
methpackagesBioC <- c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                      "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
                      "minfi",
                      "IlluminaHumanMethylationEPICmanifest",
                      "sva")

## Install from BioC
toinstallBioC <- setdiff(methpackagesBioC, installed.packages()[,1])
if(length(toinstallBioC >= 1)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(toinstallBioC, suppressUpdates = T)
  cat("finished installing new packages from BioConductor\n")
} else cat("packages we need from BioConductor are already installed\n")


## Packcages needed
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)


## Set working directory for IDATs
setwd("/home/mh158/Gen3G_850K/IDATs")
baseDir <- "/home/mh158/Gen3G_850K/IDATs"
list.files(baseDir); length(unique(list.files(baseDir)))



## read EPIC sheet from Illumina (Harvard_Hivert_EPIC_2016.09_DNASummary.csv)
setwd("/home/mh158/Gen3G_850K/")
targets <- read.metharray.sheet()

## Basenames (file name witout the R01C01...)
targets$Basename<-gsub('_.*','',targets$Barcode)
length(unique(targets$Basename))  ## 141


files<-list.files(path = baseDir)
length(files)

## Read Methylation data from IDAT files
RGset<-read.metharray.exp(file.path(baseDir, files))
save(RGset,file="RGset.RData")


############################################################################

##
##  Re-Load Data
##  Otherwise system is very slow
##  Running on Orchestra: HMS




## Sex-prediction with genomic ranges for X-Y chromosomes
load("RGset.RData")
GRset <- mapToGenome(RGset)

## Plot sex-prediction for the entire dataset
tiff("Sex_Predicted.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

plotSex(getSex(GRset, cutoff = -2))

dev.off()


# GetSex predicts sex based on X and Y chromosome methylation intensity
estSex <- getSex(GRset)
save(estSex,file="estSex.RData")

## Generate MethylSet from Genomic Ratio Set
MSet.raw <- preprocessRaw(RGset)
save(MSet.raw,file="MSet_raw.RData")


## load sample sheet information (i.e. experimental set-up)
pd<-read.table("samps_2016_09_15.csv",sep=";",dec=",",as.is=T,header=T,na.string="")
## Label Row names with Sample_Names
rownames(pd)<-paste(pd$Sentrix_ID,pd$Sentrix_Position,sep="_")

## Same dimensions?
dim(pd);dim(MSet.raw)
all(rownames(pd)==colnames(MSet.raw)) ## FAlSE!


## Sort by Unique sample ID
pd<-pd[colnames(MSet.raw),]
all(rownames(pd)==colnames(MSet.raw)) ## True!
identical(rownames(pd),colnames(MSet.raw))


## Read Phenotype from Catherine
pheno <- read.csv("phenotypes_participants_850K_corr.csv")
save(pheno, file="pheno.RData")
load("pheno.RData")
## Merge experiment Info and Phenotype
pDat<-merge(pd,pheno,by.x="Sample_Name",by.y="Patient_No_etude",sort=F)
rownames(pDat)<-paste(pDat$Sentrix_ID,pDat$Sentrix_Position,sep="_")


# Align pDat with Estimated Sex by RowNames (i.e. sample names)
pDat<-pDat[rownames(estSex),]
identical(rownames(estSex),rownames(pDat))  #True
all(rownames(estSex)==rownames(pDat))       #True
table(pDat$Delivery_Sex,estSex$predictedSex)

## Subset to relevant columns
Sex_mismatch<-pDat[,c("Delivery_Sex","Sample_Name","Sample_ID",
                    "tissueType","Sentrix_ID","Sentrix_Position",
                    "Sample_Plate","Sample_Well")]


Sex_mismatch$predictedSex<-estSex[,c("predictedSex")]
table(Sex_mismatch$predictedSex,Sex_mismatch$Delivery_Sex)
Sex_mismatch$mismatch<-ifelse(Sex_mismatch$predictedSex=="F" &
                              Sex_mismatch$Delivery_Sex=="Female","",
                       ifelse(Sex_mismatch$predictedSex=="M" &
                              Sex_mismatch$Delivery_Sex=="Male","",1))

## 16 mismatch sex
table(Sex_mismatch$mismatch)



## Save as csv:
write.csv(Sex_mismatch, file = "Sex_mismatch.csv")


## Sex-prediction with sample name for the entire Dataset
tiff("Sex_Predicted_Sample_Names.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

plotSex(getSex(GRset, cutoff = -2),id=pDat$Sample_Name)

dev.off()




################## QC-Report for all of the samples #######################

## Reload Data
load("RGset.RData")
all(colnames(RGset)==rownames(pDat))
identical(colnames(RGset),rownames(pDat))



##  QC_report for all of the samples
pDat$ID<-rownames(pDat)
qcReport(RGset, sampNames=pDat$ID, sampGroups=pDat$tissueType, pdf="qcReport_CB_Placenta_ID.pdf") #


## Densities by Tissue Type
tiff("Tissue_Densities.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

densityPlot(RGset, sampGroups=pDat$tissueType) #

dev.off()


## Densities by sample plate
tiff("Tissue_Densities.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

densityPlot(RGset, sampGroups=pDat$tissueType) #

dev.off()


## Densities by sample Plate
tiff("Plate_Densities.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

densityPlot(RGset, sampGroups=pDat$Sample_Plate) #

dev.off()








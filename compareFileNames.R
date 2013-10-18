rm(list=ls())
phenoGen <- read.table(file="/Volumes/Laura/LxS Mice/Exon Arrays/PhenoGen/PublicILSXISSRIMice.arrayFiles.txt",sep="\t",skip=1)
phenoGen$fileName = unlist(lapply(strsplit(as.character(phenoGen$V1),split="/",fixed=TRUE),function(a) a[7]))
phenoGen = phenoGen[phenoGen$fileName!="M0546_LXS22_run19.CEL",]


origFile <- read.table(file="/Volumes/Laura/LxS Mice/Exon Arrays/Source/exonArray.noOutliers.fileListing.txt",sep="\t",header=TRUE)
origFile <- origFile[!(origFile$Strain %in% c("B6","DBA")),]

dim(phenoGen)
dim(origFile)
table(sort(phenoGen$fileName)==sort(origFile$cel_files))

phenoGen[!(phenoGen$fileName %in% origFile$cel_files),]

adjusted <- read.table(file="/Volumes/Laura/LxS Mice/Exon Arrays/Data/CoreTrans/Adjusted_rma.coreTrans.csv",sep="\t",header=TRUE,row.names=1)
sampleNames <- colnames(adjusted)
sampleNames <- sampleNames[-grep("B6",sampleNames)]
sampleNames <- sampleNames[-grep("DBA",sampleNames)]
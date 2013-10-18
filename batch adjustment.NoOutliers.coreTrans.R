rm(list=ls())

setwd("/Users/laurasaba/Documents/LXS")

fileListing <- read.table("Source/exonArray.noOutliers.fileListing.21May12.txt",sep="\t",header=TRUE)
fileListing$sampleName = paste(fileListing$Strain,fileListing$Mouse,fileListing$Batch,sep=".")

sampleInfo <- cbind(fileListing[,c("sampleName","cel_files")],fileListing$Batch,fileListing$Strain,fileListing$BredAt)
colnames(sampleInfo) <- c("Sample Name","FileName","Batch","strain","Bred At")
write.table(sampleInfo,file="Data/CoreTrans/sampleInfo.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

coreTrans <- read.table("Data/CoreTrans/rma.summary.txt",sep="\t",header=TRUE)
colnames(coreTrans)<-c("ProbeID",as.character(fileListing$sampleName))
write.table(coreTrans,file="Data/CoreTrans/rma.coreTrans.csv",sep=",",row.names=FALSE,col.names=TRUE)

source('Programs/ComBat.R')
setwd("Data/CoreTrans")
ComBat('rma.coreTrans.csv','sampleInfo.txt',skip=1,type='csv')

adjusted <- read.table(file="Adjusted_rma.coreTrans.csv",sep="\t",header=TRUE,row.names=1)
adjusted <- adjusted[1:16368,]

DABG <- read.table(file="dabg.summary.txt",sep="\t",header=TRUE,row.names=1)
DABG <- DABG[1:16368,]


present <- adjusted[rowSums(DABG>0.0001)==0,]


par(mfrow=c(1,1))
plot(hclust(as.dist(1-cor(adjusted)[fileListing$Strain %in% levels(fileListing$Strain)[5:14],fileListing$Strain %in% levels(fileListing$Strain)[5:14]])))

plot(hclust(as.dist(1-cor(present)[fileListing$Strain %in% levels(fileListing$Strain)[51:64],fileListing$Strain %in% levels(fileListing$Strain)[51:64]])))
plot(hclust(as.dist(1-cor(present)[fileListing$Strain %in% c("B6","DBA","ILS","ISS"), fileListing$Strain %in% c("B6","DBA","ILS","ISS")])))
plot(hclust(as.dist(1-cor(present)[fileListing$Strain %in% c("LXS114","LXS32","LXS22"), fileListing$Strain %in% c("LXS114","LXS32","LXS22")])))

plot(hclust(as.dist(1-cor(present))))


###  Create Data Sets For PhenoGen  ###
rm(list=ls())
setwd("/Users/laurasaba/Documents/LXS")

fileListing <- read.table("Source/exonArray.noOutliers.fileListing.21May12.txt",sep="\t",header=TRUE)

# downloadable data sets without the C57 and DBA mice
adjusted <- read.table(file="Data/CoreTrans/Adjusted_rma.coreTrans.csv",sep="\t",header=TRUE,row.names=1)
DABG <- read.table(file="Data/CoreTrans/dabg.summary.txt",sep="\t",header=TRUE,row.names=1)

adjusted <- adjusted[,!(fileListing$Strain %in% c("B6","DBA"))]
DABG <- DABG[,!(fileListing$Strain %in% c("B6","DBA"))]

colnames(adjusted) <- fileListing$cel_files[!(fileListing$Strain %in% c("B6","DBA"))]
colnames(DABG) <- fileListing$cel_files[!(fileListing$Strain %in% c("B6","DBA"))]

headLines <- read.table(file="Data/CoreTrans/rma.summary.txt",sep="\t",header=FALSE,comment.char="",nrows=1000,colClasses="character")
headLines <- headLines[grep("#",headLines[,1],fixed=TRUE),1]

headLines.dabg <- read.table(file="Data/CoreTrans/dabg.summary.txt",sep="\t",header=FALSE,comment.char="",nrows=1000,colClasses="character")
headLines.dabg <- headLines.dabg[grep("#",headLines.dabg[,1],fixed=TRUE),1]

write.table(headLines,file="Data/CoreTrans/rma.coreTrans.LXS.PhenoGen.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(cbind(probeset_id = rownames(adjusted),adjusted),file="Data/CoreTrans/rma.coreTrans.LXS.PhenoGen.txt",sep="\t",row.names=FALSE,quote=FALSE,append=TRUE)

write.table(headLines.dabg,file="Data/CoreTrans/dabg.coreTrans.LXS.PhenoGen.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(cbind(probeset_id = rownames(DABG),DABG),file="Data/CoreTrans/dabg.coreTrans.LXS.PhenoGen.txt",sep="\t",row.names=FALSE,quote=FALSE,append=TRUE)

# heritability for download and filtering

rsquare <- function(example){
	Y = matrix(example,nc=1)
	J = matrix(1,nr=length(example),nc=length(example))
	n = length(example)
	X = as.matrix(model.matrix(~ -1 + factor(strains)))
	b = solve(t(X)%*%X)%*%t(X)%*%Y
	SSTO <- t(Y)%*%Y - t(Y)%*%J%*%Y / n
	SSE = t(Y)%*%Y - t(b)%*%t(X)%*%Y
	rSquare = 1 - SSE/SSTO
	return(rSquare)
	}

strains <- fileListing$Strain[!(fileListing$Strain %in% c("B6","DBA"))]
r.square <- apply(adjusted,1,rsquare)

write.table(cbind(ProbeID=names(r.square),heritability=r.square),file="Data/CoreTrans/herits.coreTrans.LXS.Brain.txt",sep="\t",row.names=FALSE,quote=FALSE)


rm(list=ls())

fileListing <- read.table("~/LXSExonData/Source/exonArray.noOutliers.fileListing.21May12.txt",sep="\t",header=TRUE)
fileListing$sampleName = paste(fileListing$Strain,fileListing$Mouse,fileListing$Batch,sep=".")

sampleInfo <- cbind(fileListing[,c("sampleName","cel_files")],fileListing$Batch,fileListing$Strain,fileListing$BredAt)
colnames(sampleInfo) <- c("Sample Name","FileName","Batch","strain","Bred At")
write.table(sampleInfo,file="sampleInfo.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

fullPS <- read.table("rma.summary.txt",sep="\t",header=TRUE)
colnames(fullPS)<-c("ProbeID",as.character(fileListing$sampleName))
write.table(fullPS,file="rma.fullPS.csv",sep=",",row.names=FALSE,col.names=TRUE)

source('~/LXSExonData/Programs/ComBat.R')
ComBat('rma.fullPS.csv','sampleInfo.txt',skip=1,type='csv')

###  Create Data Sets For PhenoGen  ###
rm(list=ls())

fileListing <- read.table("~/LXSExonData/Source/exonArray.noOutliers.fileListing.21May12.txt",sep="\t",header=TRUE)

# downloadable data sets without the C57 and DBA mice
adjusted <- read.table(file="Adjusted_rma.fullPS.csv",sep="\t",header=TRUE,row.names=1)
DABG <- read.table(file="dabg.summary.txt",sep="\t",header=TRUE,row.names=1)

adjusted <- adjusted[,!(fileListing$Strain %in% c("B6","DBA"))]
DABG <- DABG[,!(fileListing$Strain %in% c("B6","DBA"))]

colnames(adjusted) <- fileListing$cel_files[!(fileListing$Strain %in% c("B6","DBA"))]
colnames(DABG) <- fileListing$cel_files[!(fileListing$Strain %in% c("B6","DBA"))]

headLines <- read.table(file="rma.summary.txt",sep="\t",header=FALSE,comment.char="",nrows=1000,colClasses="character")
headLines <- headLines[grep("#",headLines[,1],fixed=TRUE),1]

headLines.dabg <- read.table(file="dabg.summary.txt",sep="\t",header=FALSE,comment.char="",nrows=1000,colClasses="character")
headLines.dabg <- headLines.dabg[grep("#",headLines.dabg[,1],fixed=TRUE),1]

write.table(headLines,file="rma.fullPS.LXS.PhenoGen.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(cbind(probeset_id = rownames(adjusted),adjusted),file="rma.fullPS.LXS.PhenoGen.txt",sep="\t",row.names=FALSE,quote=FALSE,append=TRUE)

write.table(headLines.dabg,file="dabg.fullPS.LXS.PhenoGen.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(cbind(probeset_id = rownames(DABG),DABG),file="dabg.fullPS.LXS.PhenoGen.txt",sep="\t",row.names=FALSE,quote=FALSE,append=TRUE)

# heritability for download and filtering
#rm(list=ls())

#fileListing <- read.table("~/LXSExonData/Source/exonArray.noOutliers.fileListing.txt",sep="\t",header=TRUE)
#adjusted <- read.table(file="Adjusted_rma.fullPS.csv",sep="\t",header=TRUE,row.names=1)

#adjusted <- adjusted[,!(fileListing$Strain %in% c("B6","DBA"))]

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
#save(r.square, file="LXS.fullPS.herits.Rdata")
write.table(cbind(ProbeID=names(r.square),heritability=r.square),file="herits.fullPS.LXS.Brain.txt",sep="\t",row.names=FALSE,quote=FALSE)


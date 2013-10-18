cd /data/Tabastore3/LauraS/LXS/Source

R

rm(list=ls())
options(stringsAsFactors=FALSE)

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

generateFiles <- function(folder,fileSuffix){
	setwd("/data/Tabastore3/LauraS/LXS/Source")
	fileListing <- read.table("exonArray.noOutliers.fileListing.21May12.txt",sep="\t",header=TRUE)
	fileListing$sampleName = paste(fileListing$Strain,fileListing$Mouse,fileListing$Batch,sep=".")

	sampleInfo <- cbind(fileListing[,c("sampleName","cel_files")],fileListing$Batch,fileListing$Strain,fileListing$BredAt)
	colnames(sampleInfo) <- c("Sample Name","FileName","Batch","strain","Bred At")
	write.table(sampleInfo,file=paste(folder,"/sampleInfo.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

	coreTrans <- read.table(paste(folder,"/rma.summary.txt",sep=""),sep="\t",header=TRUE)
	colnames(coreTrans)<-c("ProbeID",as.character(fileListing$sampleName))
	write.table(coreTrans,file=paste(folder,"/rma.summary.csv",sep=""),sep=",",row.names=FALSE,col.names=TRUE)

	source('/data/Tabastore3/LauraS/LXS/Programs/ComBat.R')
	setwd(folder)
	ComBat('rma.summary.csv','sampleInfo.txt',skip=1,type='csv')

	adjusted <- read.table(file="Adjusted_rma.summary.csv",sep="\t",header=TRUE,row.names=1)
	DABG <- read.table(file="dabg.summary.txt",sep="\t",header=TRUE,row.names=1)

	# downloadable data sets without the C57 and DBA mice
	adjusted <- adjusted[,!(fileListing$Strain %in% c("B6","DBA"))]
	DABG <- DABG[,!(fileListing$Strain %in% c("B6","DBA"))]

	colnames(adjusted) <- fileListing$cel_files[!(fileListing$Strain %in% c("B6","DBA"))]
	colnames(DABG) <- fileListing$cel_files[!(fileListing$Strain %in% c("B6","DBA"))]

	headLines <- read.table(file="rma.summary.txt",sep="\t",header=FALSE,comment.char="",nrows=1000,colClasses="character")
	headLines <- headLines[grep("#",headLines[,1],fixed=TRUE),1]

	headLines.dabg <- read.table(file="dabg.summary.txt",sep="\t",header=FALSE,comment.char="",nrows=1000,colClasses="character")
	headLines.dabg <- headLines.dabg[grep("#",headLines.dabg[,1],fixed=TRUE),1]

	write.table(headLines,file=paste("/data/Tabastore3/LauraS/ForPhenoGen/LXS.Brain.",folder,"/rma.",fileSuffix,".txt",,sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
	write.table(cbind(probeset_id = rownames(adjusted),adjusted),file=paste("/data/Tabastore3/LauraS/ForPhenoGen/LXS.Brain.",folder,"/rma.",fileSuffix,".txt",,sep=""),sep="\t",row.names=FALSE,quote=FALSE,append=TRUE)

	write.table(headLines.dabg,file=paste("/data/Tabastore3/LauraS/ForPhenoGen/LXS.Brain.",folder,"/dabg.",fileSuffix,".txt",,sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
	write.table(cbind(probeset_id = rownames(DABG),DABG),file=paste("/data/Tabastore3/LauraS/ForPhenoGen/LXS.Brain.",folder,"/dabg.",fileSuffix,".txt",,sep=""),sep="\t",row.names=FALSE,quote=FALSE,append=TRUE)

	# heritability for download and filtering
	strains <- fileListing$Strain[!(fileListing$Strain %in% c("B6","DBA"))]
	r.square <- apply(adjusted,1,rsquare)

	write.table(cbind(ProbeID=names(r.square),heritability=r.square),file=paste("/data/Tabastore3/LauraS/ForPhenoGen/LXS.Brain.",folder,"/herits.",fileSuffix,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
}

generateFiles(folder="coreTrans.mm10",fileSuffix="coreTrans.LXS.PhenoGen.mm10")
generateFiles(folder="fullTrans.mm10",fileSuffix="fullTrans.LXS.PhenoGen.mm10")
generateFiles(folder="fullPS.mm10",fileSuffix="fullPS.LXS.PhenoGen.mm10")



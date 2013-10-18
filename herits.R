rm(list=ls())
options(stringsAsFactors=FALSE)

setwd("/Volumes/Laura/LxS Mice/Exon Arrays/Data")

exprData <- read.table(file="FullTrans/rma.fullTrans.LXS.PhenoGen.txt",sep="\t",header=TRUE,row.names=1)

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
#save(r.square, file="LXS.fullTrans.herits.Rdata")
write.table(cbind(ProbeID=rownames(r.square),heritability=r.square),file="herits.fullTrans.LXS.Brain.txt",sep="\t",row.names=FALSE,quote=FALSE)


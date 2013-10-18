rm(list=ls())
options(stringAsFactors=FALSE)

setwd("/Volumes/Laura/LxS Mice/Exon Arrays/")

eQTL <- read.table(file="eQTL/CoreTrans/Results/LXS.eQTL.coreTrans.01Jun12.txt",sep="\t",header=TRUE)

hotspots <- table(eQTL[eQTL$MaxLOD_pvalue<0.01,"MaxLOD_Locus"])

exprData <- read.table(file="Data/CoreTrans/rma.coreTrans.LXS.PhenoGen.txt",sep="\t",header=TRUE)

#reduced <- exprData[exprData$probeset_id %in% eQTL[eQTL$MaxLOD_pvalue<0.01 & eQTL$MaxLOD_Locus %in% c("rs46661326", "rs8237090"),"ProbeID"],]
reduced <- exprData[exprData$probeset_id %in% eQTL[eQTL$MaxLOD_pvalue<0.01,"ProbeID"],]
rownames(reduced) = reduced[,1]
reduced <- reduced[,-1]

strains <- as.factor(unlist(lapply(strsplit(colnames(reduced),split="_",fixed=TRUE),function(a) a[2])))

library(fastICA)

for(i in 100:125){
icaResults <- fastICA(t(reduced),n.comp=i,row.norm=TRUE)

sourceMatrix = icaResults$S
herits = apply(sourceMatrix,2,function(a) summary(lm(a ~ strains))$r.square)
print(paste(sum(herits>0.9)," source values above 0.9; ",i," sources estimated",sep=""))
}


for(i in (c(1:10)*10 + 125)){
icaResults <- fastICA(t(reduced),n.comp=i,row.norm=TRUE)

sourceMatrix = icaResults$S
herits = apply(sourceMatrix,2,function(a) summary(lm(a ~ strains))$r.square)
print(paste(sum(herits>0.9)," source values above 0.9; ",i," sources estimated",sep=""))
}

icaResults <- fastICA(t(reduced),n.comp=103,row.norm=TRUE)

sourceMatrix = icaResults$S
herits = apply(sourceMatrix,2,function(a) summary(lm(a ~ strains))$r.square)


geno <- read.table(file="eQTL/CoreTrans/QTLReaper Files/LXSgeno.txt",sep="\t",header=TRUE,comment.char="@")

colnames(geno) <- gsub("ILSXISS","LXS",colnames(geno))


genotype <- rep(NA,length(strains))
genotype[strains %in% colnames(snp)[snp=="L"]] = 1
genotype[strains %in% colnames(snp)[snp=="S"]] = 0

apply(sourceMatrix,2,function(a) summary(lm(a ~ genotype))$r.square)

strainMeans <- apply(sourceMatrix,2,function(a) aggregate(a,by=list(strains),mean)$x)
rownames(strainMeans) <- levels(strains)


strainMeans <- strainMeans[rownames(strainMeans) %in% colnames(geno),]
reducedGeno <- geno[,rownames(strainMeans)]
rownames(reducedGeno) = geno$Locus

rSquare.SNP = apply(strainMeans,2,function(b) max(apply(reducedGeno,1,function(a) summary(lm(as.numeric(b) ~ as.factor(a)))$r.square)))


a = as.character(reducedGeno[1,])
summary(lm(as.numeric(strainMeans[,1]) ~ as.factor(a)))$r.square





genotype <- rep(NA,nrow(sourceMatrix))
genotype[rownames(sourceMatrix) %in% colnames(snp)[snp=="L"]] = 1
genotype[rownames(sourceMatrix) %in% colnames(snp)[snp=="S"]] = 0


apply(sourceMatrix,2,function(a) summary(lm(a ~ genotype))$r.square)



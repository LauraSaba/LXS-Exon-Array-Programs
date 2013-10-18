rm(list=ls())

setwd("/Users/laurasaba/Documents/LXS")

###  Quality Control Checks  ###
qc.checks <- read.table(file="./Data/CoreTrans/rma.report.txt",sep="\t",header=TRUE)
numfiles <- nrow(qc.checks)

tmp <- unlist(strsplit(as.character(qc.checks$cel_files),"_"))

mouse <- tmp[grep("M",tmp)]
strain <- tmp[grep("M",tmp)+1]
run <- gsub(".CEL","",tmp[grep("un",tmp)])

qc.checks <- qc.checks[order(run,strain,mouse),]

###  Chips that deviate based on PM mean of raw intensity values  ###

qc.checks[(abs(mean(log2(qc.checks$pm_mean)) - log2(qc.checks$pm_mean))/sd(log2(qc.checks$pm_mean))) > 3, c("cel_files","pm_mean")]
qc.checks[(abs(mean(log2(qc.checks$bgrd_mean)) - log2(qc.checks$bgrd_mean))/sd(log2(qc.checks$bgrd_mean))) > 3, c("cel_files","bgrd_mean")]
qc.checks[(abs(mean(log2(qc.checks$mm_mean)) - log2(qc.checks$mm_mean))/sd(log2(qc.checks$mm_mean))) > 3, c("cel_files","mm_mean")]

###  Comparison of AUC values across samples  ###

range(qc.checks$pos_vs_neg_auc)
qc.checks[((mean(qc.checks$pos_vs_neg_auc) - qc.checks$pos_vs_neg_auc)/sd(qc.checks$pos_vs_neg_auc)) > 3,  c("cel_files","pos_vs_neg_auc")]

###  Plotting Mean Absolute Deviation of the Residuals  ###

qc.checks[((mean(qc.checks$all_probeset_mad_residual_mean) - qc.checks$all_probeset_mad_residual_mean)/sd(qc.checks$all_probeset_mad_residual_mean)) > 3, c("cel_files","all_probeset_mad_residual_mean") ]

###  Plotting Relative Log Expression  ###

qc.checks[(abs(mean(qc.checks$all_probeset_rle_mean) - qc.checks$all_probeset_rle_mean)/sd(qc.checks$all_probeset_rle_mean)) >3, c("cel_files","all_probeset_rle_mean")]
qc.checks[(abs(mean(qc.checks$all_probeset_rle_stdev) - qc.checks$all_probeset_rle_stdev)/sd(qc.checks$all_probeset_rle_stdev)) >3,c("cel_files","all_probeset_rle_stdev") ]


###  Comparison of BAC Spike Controls (hybridization and/or chip problems)  ###

qc.checks[(abs(mean(qc.checks$bac_spike_mean) - qc.checks$bac_spike_mean)/sd(qc.checks$bac_spike_mean)) > 3, c("cel_files","bac_spike_mean")]


###  Comparison of polyA Spike Controls (target prep)  ###

qc.checks[(abs(mean(qc.checks$polya_spike_mean) - qc.checks$polya_spike_mean)/sd(qc.checks$polya_spike_mean)) > 3,c("cel_files","polya_spike_mean") ]


###  Comparison of putative negative and positive controls (overall quality)  ###

qc.checks[(abs(mean(qc.checks$pos_control_mean) - qc.checks$pos_control_mean)/sd(qc.checks$pos_control_mean)) > 3, c("cel_files","pos_control_mean")]
qc.checks[(abs(mean(qc.checks$neg_control_mean) - qc.checks$neg_control_mean)/sd(qc.checks$neg_control_mean)) > 3, c("cel_files","neg_control_mean")]
qc.checks[(abs(mean(qc.checks$all_probeset_mean) - qc.checks$all_probeset_mean)/sd(qc.checks$all_probeset_mean)) > 3,c("cel_files","all_probeset_mean") ]


###  Quality Control Checks  ###
dabg.checks <- read.table(file="./Data/CoreTrans/dabg.report.txt",sep="\t",header=TRUE)

###  Comparison of DABG AUC values across samples  ###

range(dabg.checks$pos_vs_neg_auc)
dabg.checks[((mean(dabg.checks$pos_vs_neg_auc) - dabg.checks$pos_vs_neg_auc)/sd(dabg.checks$pos_vs_neg_auc)) > 3, ]

###  Comparison of putative negative and positive controls (overall quality)  ###

dabg.checks[(abs(mean(dabg.checks$pos_control_percent_called) - dabg.checks$pos_control_percent_called)/sd(dabg.checks$pos_control_percent_called)) > 3, c("cel_files","pos_control_percent_called")]
dabg.checks[(abs(mean(dabg.checks$neg_control_percent_called) - dabg.checks$neg_control_percent_called)/sd(dabg.checks$neg_control_percent_called)) > 3,c("cel_files","neg_control_percent_called") ]
dabg.checks[(abs(mean(dabg.checks$all_probeset_percent_called) - dabg.checks$all_probeset_percent_called)/sd(dabg.checks$all_probeset_percent_called)) > 3,c("cel_files","all_probeset_percent_called") ]


###  Comparison of polyA Spike Controls (target prep)  ###

dabg.checks[(abs(mean(dabg.checks$polya_spike_percent_called) - dabg.checks$polya_spike_percent_called)/sd(dabg.checks$polya_spike_percent_called)) > 3, c("cel_files","polya_spike_percent_called")]

###  Clustering Results Prior to Batch Effects Adjustment  ###

exprsData <- read.table(file="./Data/coreTrans/rma.summary.txt",sep="\t",header=TRUE,row.names=1)

plot(hclust(as.dist(1-cor(exprsData))))

#transData <- exprsData[substr(rownames(exprsData),1,1)=="7",]
#plot(hclust(as.dist(1-cor(transData))))

dabgData <- read.table(file="./Data/coreTrans/dabg.summary.txt",sep="\t",header=TRUE,row.names=1)
present <- exprsData[rowSums(dabgData>0.0001)==0,]
plot(hclust(as.dist(1-cor(present))))

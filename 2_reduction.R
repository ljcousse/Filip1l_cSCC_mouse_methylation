#duplicate samples in the two batched that were done
#combining B562 and O562

library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
#library(matrixStats)

setwd("/home/kevinr/oxRRBS/")
load("allRRBS.rdata")

####################
#combine replicates#
####################
tempB1<-allRRBS[,"B562_2"]
tempB2<-allRRBS[,"B56_2"]
RRBS_red<-allRRBS[, which(colnames(allRRBS)!="B56_2")]
RRBS_red<-RRBS_red[, which(colnames(RRBS_red)!="B562_2")]
totalReads(RRBS_red)[,"B562_1"]<-totalReads(RRBS_red)[,"B562_1"]+totalReads(tempB1)
methReads(RRBS_red)[,"B562_1"]<-methReads(RRBS_red)[,"B562_1"]+methReads(tempB1)
totalReads(RRBS_red)[,"B56_1"]<-totalReads(RRBS_red)[,"B56_1"]+totalReads(tempB2)
methReads(RRBS_red)[,"B56_1"]<-methReads(RRBS_red)[,"B56_1"]+methReads(tempB2)
colnames(RRBS_red)[5]<-"B562"
colnames(RRBS_red)[12]<-"B56"
rm(allRRBS)
allRRBS<-RRBS_red

tempB1<-allRRBS[,"O562_2"]
tempB2<-allRRBS[,"O56_2"]
RRBS_red<-allRRBS[, which(colnames(allRRBS)!="O56_2")]
RRBS_red<-RRBS_red[, which(colnames(RRBS_red)!="O562_2")]
totalReads(RRBS_red)[,"O562_1"]<-totalReads(RRBS_red)[,"O562_1"]+totalReads(tempB1)
methReads(RRBS_red)[,"O562_1"]<-methReads(RRBS_red)[,"O562_1"]+methReads(tempB1)
totalReads(RRBS_red)[,"O56_1"]<-totalReads(RRBS_red)[,"O56_1"]+totalReads(tempB2)
methReads(RRBS_red)[,"O56_1"]<-methReads(RRBS_red)[,"O56_1"]+methReads(tempB2)
colnames(RRBS_red)[21]<-"O562"
colnames(RRBS_red)[28]<-"O56"
rm(allRRBS)
allRRBS<-RRBS_red

#get rid of 168 as there were some issues with the sample
RRBS_red<-allRRBS[, which(colnames(allRRBS)!="B168")]
RRBS_red<-RRBS_red[, which(colnames(RRBS_red)!="B16")]
RRBS_red<-RRBS_red[, which(colnames(RRBS_red)!="O16")]
RRBS_red<-RRBS_red[, which(colnames(RRBS_red)!="O168")]
rm(allRRBS)
allRRBS<-RRBS_red

setwd("/home/kevinr/oxRRBS/reduced/")
save(allRRBS, file="allRRBS.rdata")
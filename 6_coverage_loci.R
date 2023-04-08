library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)

setwd("/home/kevinr/oxRRBS_qSCC/without_168/")
load("sigRRBS.rdata") 
#annotation for the dataset coverage should be checked
load("annot_siglocs.rdata")
#IDs for the loci of interest
annotresults<-as.matrix(annotresults)
siglocs<-unlist(annotresults[,1]) 

t<-totalReads(allRRBS)[siglocs,]
m<-methReads(allRRBS)[siglocs,]
cov<-totalReads(allRRBS)[siglocs,]+methReads(allRRBS)[siglocs,]
cov<-cbind(cov[unlist(annotresults[,"LocID"]),],annotresults)
write.table(cov, file="coverage_siglocs.txt",quote=F,col.names=T,row.names=T, sep="\t")

totalReads(RRBS_filterd)["7877957",]

annotresults[38,]

rownames(RRBS)[1:200]
totalReads(allRRBS)[7877950,]
methReads(allRRBS)[7877950,]

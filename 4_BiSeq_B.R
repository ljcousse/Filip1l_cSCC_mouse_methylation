#like 4_BISeq, but disregarding the O samples

library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
#library(matrixStats)

setwd("/home/kevinr/oxRRBS_qSCC/without_168/B/")
load("RRBS.rdata")

##########################
#filter by shared regions#
##########################
RRBS_filtered<-filterBySharedRegions(object=RRBS,perc.samples=1,minCov=10)
RRBS_rel<-rawToRel(RRBS_filtered)
save(RRBS_filtered, file="RRBS_filtered_cov10.rdata")

################
#data filtering#
################

#noise filtering
converror<-0.01 #conversion rate error
methFDR<-0.1
minmethlevel<-0.05
MR<-as.vector(methReads(RRBS_filtered)) 
TR<-as.vector(totalReads(RRBS_filtered)) 
MRTR<-cbind(MR,TR-MR)
reslist<-apply(MRTR,1,binom.test,p=converror,alternative="greater") 
getp<-function(reslist){return(reslist$p.value)} 
sigmeth<-unlist(lapply(reslist,getp))
sigmethfdr<-matrix(p.adjust(sigmeth,method="BH"),ncol=14) 
inclocsM<-rownames(RRBS_rel)[which(rowMax((sigmethfdr<methFDR)*(methLevel(RRBS_rel[,grep("B",colnames(RRBS_rel))])>minmethlevel))==1)] 
RRBS_red<-RRBS_filtered[inclocsM,] 
rm(RRBS)
RRBS<-RRBS_red
RRBS_rel<-rawToRel(RRBS)
save(RRBS, RRBS_rel, file="RRBS_reduced_cov10.rdata")

############
#clustering#
############
pdf("clustering_Bsamples_cov10.pdf")
Bclust<-plot(hclust(dist(t(methLevel(RRBS_rel)))),labels=colnames(methLevel(RRBS_rel)))
dev.off()

####################
#average methLevels#
####################
#used only for average methylation, skip for detection of significant loci
meth<-methLevel(RRBS_rel)

t<-colMeans(meth)[1:6]
mc<-colMeans(meth)[7:12]
avmethtest<-t.test(t, mc, paired=T)


name<-c("352", "444", "545", "562", "851", "874", "35", "44", "54", "56", "85", "87", "45", "77")
meth<-colMeans(meth, na.rm=T)
a<-meth[grep("B",names(meth))]
avmeth<-cbind(a)
rownames(avmeth)<-name
colnames(avmeth)<-"RRBS"
write.table(avmeth, file="average_meth.txt",quote=F,col.names=T,row.names=T, sep="\t")


############
#Statistics#
############
meth<-methLevel(RRBS_rel)
mC<-meth[,1:14]
name_mC<-c("mC 352", "mC 444", "mC 545", "mC 562", "mC 851", "mC 874", "mC 35", "mC 44", "mC 54", "mC 56", "mC 85", "mC 87", "mC 45", "mC 77")
colnames(mC)<-name_mC

#t-test methylation
mC_pV<-NULL
for(i in 1:nrow(mC)){
  print(100*i/nrow(mC))#progress
  mC_pV[i]<-t.test(mC[i,1:6], mC[i,7:12], paired=T)$p.value
}

s1<-mC_pV[mC_pV<=0.05]
s1<-s1[!is.na(s1)]
length(s1)

mC_adj_pV<-p.adjust(mC_pV, method="BH")
s2<-mC_adj_pV[mC_adj_pV<=0.05]
s2<-s2[!is.na(s2)]
length(s2)
sort(mC_adj_pV)[1:100]

#########################
#create methLevel matrix#
#########################
mC_mat<-cbind(mC, mC_adj_pV)
colnames(mC_mat)[15]<-"p-value TvsMC // mC"
meth_mat<-mC_mat
#save(meth_mat, file="methylation_cov10.rdata")

#significant loci
pV<-meth_mat[,15]
siglocs<-meth_mat[pV<=0.05,]




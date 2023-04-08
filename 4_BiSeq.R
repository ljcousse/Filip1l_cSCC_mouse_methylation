library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
#library(matrixStats)

setwd("/home/kevinr/oxRRBS/reduced/")
load("allRRBS.rdata")

##########################
#filter by shared regions#
##########################
RRBS_filtered<-filterBySharedRegions(object=allRRBS,perc.samples=1,minCov=10)
RRBS_rel<-rawToRel(RRBS_filtered)
#save(RRBS_filtered, file="allRRBS_filtered_cov10.rdata")

################
#data filtering#
################

#noise filtering
converror<-0.01 #conversion rate error
methFDR<-0.1
minmethlevel<-0.05
MR<-as.vector(methReads(RRBS_filtered)[,15:28]) 
TR<-as.vector(totalReads(RRBS_filtered)[,15:28]) 
MRTR<-cbind(MR,TR-MR)
reslist<-apply(MRTR,1,binom.test,p=converror,alternative="greater") 
getp<-function(reslist){return(reslist$p.value)} 
sigmeth<-unlist(lapply(reslist,getp))
sigmethfdr<-matrix(p.adjust(sigmeth,method="BH"),ncol=14) 
inclocsM<-rownames(RRBS_rel)[which(rowMax((sigmethfdr<methFDR)*(methLevel(RRBS_rel[,grep("O",colnames(RRBS_rel))])>minmethlevel))==1)] 
allRRBS_red<-RRBS_filtered[inclocsM,] 
rm(allRRBS)
allRRBS<-allRRBS_red
#save(allRRBS, file="allRRBS_reduced_cov10.rdata")

####################
#average methLevels#
####################
#used only for average methylation, skip for detection of significant loci
allRRBS_red<-rawToRel(allRRBS_red)
meth<-methLevel(allRRBS_red)

t<-colMeans(meth[,1:6])
mc<-colMeans(meth[,7:12])
avmethtest<-t.test(t, mc, paired=T)


name<-c("352", "444", "545", "562", "851", "874", "35", "44", "54", "56", "85", "87", "45", "77")
meth<-colMeans(meth, na.rm=T)
a<-meth[grep("B",names(meth))]
b<-meth[grep("O",names(meth))]
avmeth<-cbind(a,b)
rownames(avmeth)<-name
colnames(avmeth)<-c("RRBS","oxRRBS")
dif<-c(avmeth[,"RRBS"]-avmeth[,"oxRRBS"])
avmeth<-cbind(avmeth, dif)
#write.table(avmeth, file="average_meth.txt",quote=F,col.names=T,row.names=T, sep="\t")


##################################################
#separation of methylation and hydroxymethylation#
##################################################
meth<-methLevel(allRRBS_red)
Bmeth<-meth[,1:14]
Ometh<-meth[,15:28]
name_mC<-c("mC 352", "mC 444", "mC 545", "mC 562", "mC 851", "mC 874", "mC 35", "mC 44", "mC 54", "mC 56", "mC 85", "mC 87", "mC 45", "mC 77")
name_hmC<-c("hmC 352", "hmC 444", "hmC 545", "hmC 562", "hmC 851", "hmC 874", "hmC 35", "hmC 44", "hmC 54", "hmC 56", "hmC 85", "hmC 87", "hmC 45", "hmC 77")

#pure methylation
mC<-Ometh
#pure hydroxymethylation
hmC<-Bmeth-Ometh

colnames(mC)<-name_mC
colnames(hmC)<-name_hmC

############
#Statistics#
############

#t-test methylation
mC_pV<-NULL
for(i in 1:nrow(mC)){
  print(100*i/nrow(mC))#progress
  mC_pV[i]<-t.test(mC[i,1:6], mC[i,7:12], paired=T)$p.value
}
s1<-mC_pV[mC_pV<=0.01]
s1<-s1[!is.na(s1)]
length(s1)

mC_adj_pV<-p.adjust(mC_pV, method="BH")
s2<-mC_adj_pV[mC_adj_pV<=0.05]
s2<-s2[!is.na(s2)]
length(s2)
sort(mC_adj_pV)[1:100]

#t-test hydroxymethylation
hmC_pV<-NULL
for(i in 1:nrow(hmC)){
  print(100*i/nrow(hmC))#progress
  hmC_pV[i]<-t.test(hmC[i,1:6], hmC[i,7:12], paired=T)$p.value
}
s1<-hmC_pV[hmC_pV<=0.01]
s1<-s1[!is.na(s1)]
length(s1)

hmC_adj_pV<-p.adjust(hmC_pV, method="BH")
s2<-hmC_adj_pV[hmC_adj_pV<=0.05]
s2<-s2[!is.na(s2)]
length(s2)
sort(hmC_adj_pV)[1:100]

#########################
#create methLevel matrix#
#########################
mC_mat<-cbind(mC, mC_adj_pV)
colnames(mC_mat)[15]<-"p-value TvsMC // mC"
hmC_mat<-cbind(hmC, hmC_adj_pV)
colnames(hmC_mat)[15]<-"p-value TvsMC // hmC"
meth_mat<-cbind(mC_mat, hmC_mat)
#save(meth_mat, file="methylation_cov10.rdata")

#significant loci
pV<-meth_mat[,15]
siglocs<-meth_mat[pV<=0.05,]
siglocs<-siglocs[,1:15]




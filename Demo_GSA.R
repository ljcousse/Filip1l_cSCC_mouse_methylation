#set the directory with the tumour data
setwd("C:/Users/TimDeMeyer/Documents/UGent/Education/Vakken/AHTA/1617/PracticalSessions/Practicum 4/")

require(limma)
require(edgeR)

rnadatafull<-read.table("SalivaryGlandCancer.txt",header=T,sep="\t",quote="")
#removing duplicated rows
rnadatafull<-rnadatafull[duplicated(rnadatafull[,5:8])==F,]
#optionally:
rnadatafull<-rnadatafull[rowSums(rnadatafull[,9:14])>3,]
rownames(rnadatafull)<-1:nrow(rnadatafull)
rnadata<-rnadatafull[,9:14]

individuals<-factor(c("A","B","C","A","B","C"))
cvst<-factor(c("Control","Control","Control","Tumor","Tumor","Tumor"))

############
#limma-voom#
############
#cf. EdgeR
y<-DGEList(rnadata,genes=rnadatafull[,1:8])
y<-calcNormFactors(y)
y$samples
plotMDS(y)

design<-model.matrix(~individuals+cvst) #individuals as well
rownames(design)<-colnames(y)

#voom
v<-voom(y,design,plot=TRUE)
fit<-lmFit(v,design)
fit<-eBayes(fit)
topTable(fit,coef=4)
sum(topTable(fit,coef=4,nrow(fit))$adj.P.Val<0.05)
sum(topTable(fit,coef=4,nrow(fit))$adj.P.Val<0.1)
sum(topTable(fit,coef=4,nrow(fit))$adj.P.Val<0.2)

#gene set analysis
#First, obtain gene sets, get indices (typically Entrez ID based)
#cf. http://software.broadinstitute.org/gsea/msigdb/collections.jsp#H
download.file("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata","human_H_v5p2.rdata", mode = "wb")
load("human_H_v5p2.rdata")
H.indices<-ids2indices(Hs.H,y$genes$EntrezGene.ID)

#romer
r<-romer(v,index=H.indices,design=design,contrast=2)
topRomer(r,alternative="mixed") #Note: result not FDR adjusted


#goana
go<-goana(fit,coef=4,FDR=0.2,trend=T,plot=T)
topGO(go)
nrow(go)
kegg<-kegga(fit,coef=4,FDR=0.2,trend=T,plot=T)
topKEGG(kegg)
nrow(kegg)

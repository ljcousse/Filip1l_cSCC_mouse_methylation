oxRRBSannot<-function(siglocs){
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomicFeatures)
  library(BiSeq)
  library(biomaRt)
  
  #load("/home/kevinr/oxRRBS_qSCC/reduced/allrrbsrowranges.rdata") #loads "locs", consisting of the saved rowranges object (cf. my mail), i.e. perform the rowranges function on your full rrbs dataset to extract the location data and save it in an rdata object
  #siglocs<-locs[rownames(siglocs),]  #siglocs represent the significant loci, these are the names of each location (typically integers, but here put as characters)
  siglocs <- locs
  gtf_dir <- "/home/kevinr/oxRRBS_qSCC/annotation/Mus_musculus.GRCm38.84.gtf" # location GTF file on MIDAS
  txdb <- makeTxDbFromGFF(gtf_dir, format="gtf") #import gtf object
  
  prom<-promoters(txdb,upstream=2000,downstream=500,columns=c("GENEID","TXNAME")) #get CpGs located in promoter regions (2500 upstream, 500 nt downst) with associated geneIDs & transcriptIDs
  ex<-exons(txdb,columns=c("GENEID","EXONRANK")) #get geneID & rank for each CpG located in an exon
  ex$EXONRANK<-CharacterList(ex$EXONRANK) #some preproc
  gen<-genes(txdb,columns="GENEID") #get full gene (i.e. complement to intergenic regions)
  
  #here, the annotations loaded higher are overlapped with the significant loci, and for those loci additional columns are added to the siglocs object (different columns for different properties metnioned higher) 
  rrbsannot<-annotateGRanges(object=siglocs,regions=ex,name="Exonic_geneID",regionInfo="GENEID")
  rrbsannot<-annotateGRanges(object=rrbsannot,regions=prom,name="Promoter_geneID",regionInfo="GENEID")
  rrbsannot<-annotateGRanges(object=rrbsannot,regions=gen,name="geneID",regionInfo="GENEID")
  rrbsannot<-annotateGRanges(object=rrbsannot,regions=ex,name="Exon_rank",regionInfo="EXONRANK")
  rrbsannot<-annotateGRanges(object=rrbsannot,regions=prom,name="Promoter_transcriptID",regionInfo="TXNAME")
  
  #extract from this object the information I really need (many redundant columns)
  annotfile<-data.frame(chrom=seqnames(rrbsannot),ranges=ranges(rrbsannot),FuncAnnot=mcols(rrbsannot))
  rownames(annotfile)<-annotfile$ranges.name #add some rownamaes
  colnames(annotfile)[grep("ranges.start",colnames(annotfile))]<-"pos" #add some columninfo
  annotfile<-annotfile[,-c(grep("ranges",colnames(annotfile)))] #remove some uninteresting data
  
  annotresults<-vector("list",nrow(annotfile))
  print("Annotating results")
  
  #don't bother too much about the lines below, some advanced preprocessing
  #in a nutshell: higher, a single column is used per feature, whereas we want a single column with the annotation (and some additional columns with the extra information)
  for (i in 1:length(annotresults)){
    annot<-annotfile[i,]
    annotgenes<-unique(unlist(strsplit(grep("MUSG",annot,value=T),split=",")))
    annotgenes[is.null(annotgenes)]<-"Intergenic"
    numannot<-length(annotgenes)
    
    sname<-matrix(nrow=numannot,ncol=2)
    sname[,1]<-annotgenes
    for (j in 1:numannot){
      sname[j,2]<-c(NA,NA,"Exonic","Promoter","Intronic")[grep(sname[j,1],annot)[1]]
    }
    annotresults[[i]]<-cbind(rep(rownames(annot),numannot),rep(as.character(annot$chrom),numannot),rep(annot$pos,numannot),sname,matrix(rep(annot[3:length(annot)],numannot),nrow=numannot,byrow=T))
  }
  
  #Cleaning up: for some additional data columns, the information is not relevant (e.g. exonrank when considering it as a promoter)
  #Also, FunctionalAnnotation2 is created = in case CpGs are located in both promoter and exon (as these may overlap), exon is given in FunctionalAnnotation1, but promoter is given in FunctionalAnnotation2 (rest is the same)
  print("Cleaning up annotation")
  annotresults<-do.call(rbind,annotresults)
  colnames(annotresults)<-c("LocID",colnames(annotfile)[1:2],"GeneID","FunctionalAnnotation1",colnames(annot)[3:7])
  annotresults[annotresults[,"GeneID"]=="Intergenic","FunctionalAnnotation1"]<-"Intergenic"
  FunctionalAnnotation2<-annotresults[,"FunctionalAnnotation1"] #Some promoters are also exonic (interesting for 5'UTR, latter information is however unavailable in GFF T gUttata)
  FunctionalAnnotation2[FunctionalAnnotation2=="Exonic"&unlist(annotresults[,"FuncAnnot.Promoter_geneID"])==unlist(annotresults[,"GeneID"])]<-"Promoter"
  annotresults<-cbind(annotresults[,1:5],FunctionalAnnotation2,annotresults[,9:10])
  annotresults[annotresults[,"FunctionalAnnotation1"]=="Promoter"&annotresults[,"FunctionalAnnotation2"]=="Promoter","FuncAnnot.Exon_rank"]<-"NA"
  annotresults[annotresults[,"FunctionalAnnotation1"]=="Exonic"&annotresults[,"FunctionalAnnotation2"]=="Exonic","FuncAnnot.Promoter_transcriptID"]<-"NA"
  
  #Here, we use the biomart tool to extract gene symbols associated with Ensembl IDs
  print("Fetch Gene Symbols")
  mart <- useMart("ensembl")
  listDatasets(mart=mart) 
  mart2<-useDataset("mmusculus_gene_ensembl", mart = mart) 
  genes.table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),values=annotresults[,"GeneID"], mart= mart2)
  genes.table <- rbind(genes.table,c("Intergenic","NA"))
  #nosymbol<-unique(setdiff(annotresults[,"GeneID"],genes.table[,1]))
  rownames(genes.table)<-genes.table[,1]
  genes.table<-as.matrix(genes.table)
  
  #GeneSymbols are added to the outcome data
  GeneSymbol<-genes.table[unlist(annotresults[,4]),2]
  annotresults<-cbind(annotresults[,1:4],GeneSymbol,annotresults[,5:ncol(annotresults)])
  return(annotresults)
}

#annotresults from annotation script
annotresults<-as.matrix(annotresults)
siglocs<-unlist(annotresults[,1])
sigRRBS<-allRRBS[rownames(allRRBS)%in%siglocs,] 
sigRRBS<-sigRRBS[order(rownames(sigRRBS)),]
save(sigRRBS, file="sigRRBS.rdata")

load("annot_whole_set.rdata") #annotation for all loci
load("methylation_cov10.rdata") #meth_mat, containing methLevels and p-values for methylation and hydroximethylation
methylation<-cbind(meth_mat[unlist(annotresults[,"LocID"]),],annotresults)
save(methylation, file="methylation_cov10.rdata")
write.table(methylation, file="sigRRBS.txt",quote=F,col.names=T,row.names=T, sep="\t")

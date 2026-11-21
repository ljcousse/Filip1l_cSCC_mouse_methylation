##################
# Analysis Kevin #
##################

#duplicate samples in the two batched that were done
#combining B562 and O562

library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
library(biomaRt)
library(yaml)

# Directories are read from environment.yml (copy environment_example.yml and fill in the local paths)
if (!file.exists("environment.yml")){
  stop("environment.yml not found, copy environment_example.yml and fill in the local paths")
}
env <- read_yaml("environment.yml")

dir_work <- env$dir_work
setwd(dir_work)

# parameters
mincov <- 8
minmethreads <- 6
ct <- 0.01

load("allRRBS.rdata")

print(head(allRRBS))


allRRBS_repl <- allRRBS

# addition technical replicates
head(totalReads(allRRBS))

############
#clustering#
############

# Annotation with replicates
f1 <- c(rep('T',8),rep('C',10))
f2 <- c(rep(c(1,2,3,4,5,5,6,7),2),8,9)
f1 <- as.factor(f1)
f2 <- as.factor(f2)

sample_info <- t(rbind(colnames(allRRBS),as.character(f1),as.character(f2)))
write.table(sample_info,"sample_info.txt",col.names = F,row.names = F,sep="\t",quote = F)

# m_matrix <- matrix(as.matrix(methLevel(rrbs_repl.M)[apply(methLevel(rrbs_repl.M),1,var)>sort(apply(methLevel(rrbs_repl.M),1,var),decreasing = T)[10000],]),ncol=18)
m_matrix <- matrix(as.matrix(methLevel(rawToRel(allRRBS[,grep("B",colnames(allRRBS))]))),ncol=18)
m_matrix <- m_matrix[!is.na(rowSums(m_matrix)),]
cn <- paste(f1,f2,sep=",")
cn[duplicated(cn)] <- paste(cn[duplicated(cn)],c("*","*"),sep="")
colnames(m_matrix) <- cn
print(cn)

correl <- 1- abs(cov(m_matrix[apply(m_matrix,1,var)>sort(apply(m_matrix,1,var),decreasing = T)[10001],]))
distance <- as.dist(correl)
hclust_complete <- hclust(distance, method = "complete")
dend_complete <- as.dendrogram(hclust_complete)

jpeg("Clustering.jpg")
par(xpd=TRUE)
par(mar=c(2.5,2.5,2.5,11))
Bclust<-plot(hclust_complete)
# legend(20,0,legend=c(paste(colnames(RRBS_repl),"is",cn)),ncol = 1)
dev.off()
print("Clusterplot made")

# Merge technical replicates
totalReads(allRRBS)[,grep("_1",colnames(allRRBS))] <- totalReads(allRRBS)[,grep("_1",colnames(allRRBS))] + totalReads(allRRBS)[,grep("_2",colnames(allRRBS))]
methReads(allRRBS)[,grep("_1",colnames(allRRBS))] <- methReads(allRRBS)[,grep("_1",colnames(allRRBS))] + methReads(allRRBS)[,grep("_2",colnames(allRRBS))]
allRRBS <-allRRBS[,c(-grep("_2",colnames(allRRBS)))]
colnames(allRRBS) <- unlist(strsplit(colnames(allRRBS),"_1"))
head(totalReads(allRRBS))

# Annotation without replicates
f1 <- c(rep('T',7),rep('C',9))
f2 <- c(rep(1:7,2),8,9)
f1 <- as.factor(f1)
f2 <- as.factor(f2)

allRRBS_B <- allRRBS[,grep("B",colnames(allRRBS))]
allRRBS_O <- allRRBS[,grep("O",colnames(allRRBS))]

print(dim(allRRBS_B))
RRBS_B <- filterBySharedRegions(object=allRRBS_B,perc.samples=1,minCov=mincov)
print(dim(RRBS_B))


### 
# Coverage plot
### 

rrbs_B.rel_nof <- rawToRel(allRRBS_B)

data_lm_nof <- NULL
for(i in 1:ncol(methLevel(rrbs_B.rel_nof))){
  data_lm_nof <- rbind(data_lm_nof,cbind(mean(methLevel(rrbs_B.rel_nof)[,i],na.rm=T),f1[i],f2[i]))
}

data_lm_nof <- as.data.frame(data_lm_nof,stringsAsFactors=F)

colnames(data_lm_nof) <- c("y","f1","f2")
data_lm_nof$y <- as.numeric(data_lm_nof$y)
data_lm_nof$f1 <- as.factor(data_lm_nof$f1)
data_lm_nof$f2 <- as.factor(data_lm_nof$f2)

jpeg("Methylation_plot_nofilter.jpg", width = 2480, height = 1240)
par(mfrow=c(1,2))
boxplot(formula=y~f1,names=c("Control","Tumor"),data=data_lm_nof,xlab="Biopt",ylab="Methylation percentage")
boxplot(formula=y~f2,names=paste('id',levels(f2)),las=2,data=data_lm_nof,xlab="Biopt",ylab="Total reads")
dev.off()

data_boxplot_final <- data_lm_nof
data_boxplot_final$y <- data_boxplot_final$y*100
tiff("Methylation_plot_paper.tiff",height=5,width=4,units="in",res=150)
boxplot(formula=y~f1,names=c("Control","Tumor"),
        data=data_boxplot_final,xlab="",ylab="Fraction methylation (%)")
dev.off()

print("model 1, additive")
anova_diffmeth_1_nof <- aov(y~f1+f2,data=data_lm_nof)
print(summary(anova_diffmeth_1_nof))

print("model 1, only f1")
anova_diffmeth_2_nof <- aov(y~f1,data=data_lm_nof)
print(summary(anova_diffmeth_2_nof))

t.test(data_lm_nof$y[(data_lm_nof$f1==1) & !(data_lm_nof$f2%in%c(8,9))],
       data_lm_nof$y[(data_lm_nof$f1==2) & !(data_lm_nof$f2%in%c(8,9))],paired=T)



inclocsM <- names(which(rowSums(methReads(RRBS_B))>=minmethreads))
RRBS_B <- RRBS_B[inclocsM,]
print(dim(RRBS_B))

locs_B <- rowRanges(RRBS_B)
locs_O <- rowRanges(allRRBS_O)
locs_repl <- rowRanges(allRRBS_repl)

RRBS_O <- allRRBS_O[paste(locs_O@seqnames,locs_O@ranges@start,sep=",")%in%paste(locs_B@seqnames,locs_B@ranges@start,sep=","),]
RRBS_repl <- allRRBS_repl[paste(locs_repl@seqnames,locs_repl@ranges@start,sep=",")%in%paste(locs_B@seqnames,locs_B@ranges@start,sep=","),grep("B",colnames(allRRBS))]
print(dim(RRBS_B))
print(dim(RRBS_O))
print(dim(RRBS_repl))

rrbs_B.rel <- rawToRel(RRBS_B)
rrbs_O.rel <- rawToRel(RRBS_O)
rrbs_repl.rel <- rawToRel(RRBS_repl)


## Check Cmss1
# load("genes.Rda")
# genes[genes$mgi_symbol=="Cmss1",]
# 
# locs_B_df <- as.data.frame(locs_B)
# locs_B_df$seqnames <- as.character(locs_B_df$seqnames)
# 
# methlevel_cmss1 <- cbind(locs_B_df[(locs_B_df$seqnames==16) & (locs_B_df$start> 57302199) & (locs_B_df$start < 57606825),1:2],
#                          methLevel(rrbs_B.rel)[(locs_B_df$seqnames==16) & (locs_B_df$start> 57302199) & (locs_B_df$start < 57606825),])
# cvst <- c(rep('T',7),rep('C',9))
# ratio <- rowMeans(methlevel_cmss1[3:18][,grepl("C",cvst)])>rowMeans(methlevel_cmss1[3:18][,!grepl("C",cvst)])


write.table(methLevel(rrbs_B.rel),file="Betas_rrbs.txt",sep="\t",col.names=T,row.names = F,quote = F)
write.table(as.data.frame(locs_B)[,1:2],file="Annot_rrbs.txt",sep="\t",col.names=T,row.names = F,quote = F)

diff <- methLevel(rrbs_B.rel)-methLevel(rrbs_O.rel)

breaks <- seq(-1.05,1.05,0.1)
jpeg(file="RRBS_hydroxy.jpg",width = 1000, height = 1000, units = "px")
hist_diff <- hist(diff,breaks=breaks,xlab="Hydroxymethylation",cex.lab=2,cex.axis=2,cex.main=2)
dev.off()
## Hypomethylation mainly noise: continue with methylation files 

## M Values
###########

rrbs_B.M<-rrbs_B.rel 
for (i in 1:ncol(rrbs_B.rel)){
  ML<-methLevel(rrbs_B.rel)[,i]
  methLevel(rrbs_B.M)[,i]<-log2((ML+ct)/(1-ML+ct))
}
rrbs_repl.M<-rrbs_repl.rel 
for (i in 1:ncol(rrbs_repl.rel)){
  ML<-methLevel(rrbs_repl.rel)[,i]
  methLevel(rrbs_repl.M)[,i]<-log2((ML+ct)/(1-ML+ct))
}

## Beta values
##########

rrbs_B.betas <- methLevel(rrbs_B.rel)




############
#Annotation#
############

locs<-rowRanges(rrbs_B.rel)

if (!file.exists("genes.Rda")){
  ensembl = useMart("ensembl",host='jan2019.archive.ensembl.org')
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  genes <- getBM(attributes = c('ensembl_gene_id','mgi_symbol','chromosome_name',
                                'start_position','end_position','transcription_start_site'),mart = ensembl)
  colnames(genes) <- c("gene","gene_symbol","chr","start","end","TSS")
  dim(genes) 
  genes <- genes[!duplicated(genes$gene),]
  dim(genes)
  save(genes,file="genes.Rda")
} else {
  load("genes.Rda")
}


# genes <- getBM(attributes = c('ensembl_gene_id','entrezgene','mgi_symbol','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
# load("genes.Rda")
# colnames(genes) <- c("gene","entrez_gene","gene_symbol","chr","start","end","TSS")
# dim(genes) 
# genes <- genes[!duplicated(genes$gene),]
# dim(genes)

# Promoter & genic
sense <- abs(genes$start-genes$TSS)<abs(genes$end-genes$TSS)

new_start <- cbind(c(genes$start[sense]-1999,genes$end[!sense]-499),
                   c(genes$start[sense]-1999,genes$start[!sense]+1))
new_start <- apply(new_start,1,max)
new_end <- cbind(c(genes$start[sense]+499,genes$end[!sense]+1999),
                 c(genes$end[sense]-1,genes$end[!sense]+1999))
new_end <- apply(new_end,1,min)
prom_gr <- GRanges(data.frame(chr=c(genes$chr[sense],genes$chr[!sense]),
                              start=new_start,
                              end=new_end,
                              gene=c(genes$gene[sense],genes$gene[!sense]),
                              gene_symbol=c(genes$gene_symbol[sense],genes$gene_symbol[!sense])))
overlap_p <- subsetByOverlaps(locs,prom_gr)
promoter <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_p))

new_start <- cbind(c(genes$start[sense]+500,genes$start[!sense]),
                   c(genes$end[sense],genes$start[!sense]))
new_start <- apply(new_start,1,min)
new_end <- cbind(c(genes$end[sense],genes$end[!sense]-500),
                 c(genes$end[sense],genes$start[!sense]))
new_end <- apply(new_end,1,max)
genic_gr <- GRanges(data.frame(chr=c(genes$chr[sense],genes$chr[!sense]),
                               start=new_start,
                               end=new_end,
                               gene=c(genes$gene[sense],genes$gene[!sense]),
                               gene_symbol=c(genes$gene_symbol[sense],genes$gene_symbol[!sense])))
overlap_g <- subsetByOverlaps(locs,genic_gr)
genic <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_g))

genic[genic & promoter] <- FALSE

# CGI
cgi_data <- read.table("Mouse_CGI.txt",sep="\t",header=T,stringsAsFactors = F)
cgi_data$chrom <- unlist(strsplit(cgi_data$chrom,"chr"))[unlist(strsplit(cgi_data$chrom,"chr"))!=""]

cgi_gr <- GRanges(data.frame(chr=cgi_data$chrom,
                              start=cgi_data$chromStart,
                              end=cgi_data$chromEnd))
overlap_cgi <- subsetByOverlaps(locs,cgi_gr)
cgi <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_cgi))

shore_gr <- GRanges(data.frame(chr=c(cgi_data$chrom,cgi_data$chrom),
                               start=c(cgi_data$chromStart-2000,cgi_data$chromEnd+1),
                               end=c(cgi_data$chromStart-1,cgi_data$chromEnd+2000)))
overlap_shore <- subsetByOverlaps(locs,shore_gr)
shore <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_shore))

shelf_gr <- GRanges(data.frame(chr=c(cgi_data$chrom,cgi_data$chrom),
                               start=c(cgi_data$chromStart-4000,cgi_data$chromEnd+2001),
                               end=c(cgi_data$chromStart-2001,cgi_data$chromEnd+4000)))
overlap_shelf <- subsetByOverlaps(locs,shelf_gr)
shelf <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_shelf))

shore[shore & cgi] <- FALSE
shelf[shelf & cgi] <- FALSE
shelf[shelf & shore] <- FALSE
intergenic <- rep(T,length(cgi))
intergenic[cgi | shore | shelf] <- FALSE

# CGI context per CpG (column cgi_sort in output_rrbs.txt)
cgi_sort <- rep("open_sea",length(cgi))
cgi_sort[shelf] <- "shelf"
cgi_sort[shore] <- "shore"
cgi_sort[cgi] <- "island"


# Gene and gene symbol annotation
overlaps_genic <- findOverlaps(locs,genic_gr)
overlaps_prom <- findOverlaps(locs,prom_gr)

all_overlap <- c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene),as.character(prom_gr[subjectHits(overlaps_prom)]$gene))
all_overlap <- cbind(all_overlap,c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene_symbol),as.character(prom_gr[subjectHits(overlaps_prom)]$gene_symbol)))
all_overlap <- cbind(all_overlap,c(locs[queryHits(overlaps_genic)]@ranges@NAMES,locs[queryHits(overlaps_prom)]@ranges@NAMES))
colnames(all_overlap) <- c("gene","gene_symbol","rowname")


rn <- as.matrix(locs@ranges@NAMES,ncol=1)
colnames(rn) <- "rowname"

merger <- merge(rn,all_overlap,all=T)
merger <- merger[sort(as.numeric(as.character(merger$rowname)),index.return=T)$ix,]

merger$gene <- as.character(merger$gene)
merger$gene_symbol <- as.character(merger$gene_symbol)

gene <- as.character(merger$gene[!duplicated(merger$rowname)])
gene_symbol <- as.character(merger$gene_symbol[!duplicated(merger$rowname)])

max_n <- max(table(merger$rowname))
print(max_n)
gene <- cbind(gene,matrix(rep(NA,length(gene)*(max_n-1)),ncol=max_n-1))
gene_symbol <- cbind(gene_symbol,matrix(rep(NA,length(gene_symbol)*(max_n-1)),ncol=max_n-1))
for (i in 2:max_n){
  n_times <- merger[merger$rowname%in%names(table(merger$rowname)[table(merger$rowname)==i]),]
  if(nrow(n_times)==0){
    print(paste("no locus with ",i,"-fold annotation",sep=""))
    next
  }
  # print(head(n_times))
  # print(dim(n_times))
  # print(head(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]])))
  # print(length(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]])))
  # print(head(c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]))
  # print(length(c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]))
  # print(head(-c(c(1:(nrow(n_times)/i))*i-(i-1))))
  # print(tail(-c(c(1:(nrow(n_times)/i))*i-(i-1))))
  # print(length(-c(c(1:(nrow(n_times)/i))*i-(i-1))))
  
  if (i==2 | nrow(n_times)==i){
    # print(length(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]])))
    # print(length(gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i]))
    gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
    gene_symbol[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
  } else {
    # print(head(t(matrix(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]),nrow=i-1))))
    
    # print(dim(t(matrix(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]),nrow=i-1))))
    # print(dim(gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i]))
    gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
    gene_symbol[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
  }
  
}



################################
#Boxplot Methylation percentage#
################################

# Annotation without replicates
f1 <- c(rep('T',7),rep('C',9))
f2 <- c(rep(1:7,2),8,9)
f1 <- as.factor(f1)
f2 <- as.factor(f2)

print("Control of sample annotation")
print(rbind(colnames(RRBS_B),as.character(f1),as.character(f2)))

print("###############")
print("# Methylation #")
print("###############")

data_lm <- NULL
for(i in 1:ncol(methLevel(rrbs_B.rel))){
  data_lm <- rbind(data_lm,cbind(mean(methLevel(rrbs_B.rel)[,i]),f1[i],f2[i]))
}

data_lm <- as.data.frame(data_lm,stringsAsFactors=F)

colnames(data_lm) <- c("y","f1","f2")
data_lm$y <- as.numeric(data_lm$y)
data_lm$f1 <- as.factor(data_lm$f1)
data_lm$f2 <- as.factor(data_lm$f2)

jpeg("Methylation_plot.jpg", width = 2480, height = 1240)
par(mfrow=c(1,2))
boxplot(formula=y~f1,names=c("Control","Tumor"),data=data_lm,xlab="Biopt",ylab="Methylation percentage")
boxplot(formula=y~f2,names=paste('id',levels(f2)),las=2,data=data_lm,xlab="Biopt",ylab="Total reads")
dev.off()

tiff("Methylation_plot_paper_small.tiff",units="in", width=3, height=4,res=300)
boxplot(formula=y~f1,names=c("Control","Tumor"),data=data_lm,las=1,xlab="Biopt",
        ylab="Methylation percentage",cex.lab=0.6,cex.axis=0.5,cex.main=0.6,cex=0.5,lwd=0.5,bowlwd=0.5,
        xaxt="n",yaxt="n",axes=F)
axis(1, at=c(1:2),labels = c("Control","Tumor"),lwd.ticks=0.5,tck=-0.020,cex.axis=0.5,lwd=0.5)
axis(2, at=c(0.04,0.05,0.06,0.07,0.08),labels = c("0.04","0.05","0.06","0.07","0.08"),lwd.ticks=0.5,tck=-0.020,
     las=1,cex=0.5,lwd=0.5,cex.axis=0.5)
box(lwd=0.5)
dev.off()

print("model 1, additive")
anova_diffmeth_1 <- aov(y~f1+f2,data=data_lm)
print(summary(anova_diffmeth_1))

print("model 1, only f1")
anova_diffmeth_2 <- aov(y~f1,data=data_lm)
print(summary(anova_diffmeth_2))

t.test(data_lm$y[(data_lm$f1==1) & !(data_lm$f2%in%c(8,9))],data_lm$y[(data_lm$f1==2) & !(data_lm$f2%in%c(8,9))],paired=T)

print("############")
print("# Coverage #")
print("############")

data_lm <-cbind(colSums(totalReads(RRBS_B)),f1,f2)
data_lm <- as.data.frame(data_lm,stringsAsFactors=F)
colnames(data_lm) <- c("y","f1","f2")
data_lm$y <- as.numeric(data_lm$y)
data_lm$f1 <- as.factor(data_lm$f1)
data_lm$f2 <- as.factor(data_lm$f2)

jpeg("Coverage_plot.jpg", width = 2480, height = 1240)
par(mfrow=c(1,2))
boxplot(formula=y~f1,names=c("Control","Tumor"),data=data_lm,xlab="Biopt",ylab="Total reads")
boxplot(formula=y~f2,names=paste('id',levels(f2)),las=2,data=data_lm,xlab="Biopt",ylab="Total reads")
dev.off()

print("model 1, additive")
anova_diffcov_1 <- aov(y~f1+f2,data=data_lm)
print(summary(anova_diffcov_1))

print("model 1, only f1")
anova_diffcov_2 <- aov(y~f1,data=data_lm)
print(summary(anova_diffcov_2))

print("###############")
print("# Correlation #")
print("###############")

meth <- NULL
for(i in 1:ncol(methLevel(rrbs_B.rel))){
  meth <- c(meth,mean(methLevel(rrbs_B.rel)[,i]))
}
cover <- colSums(totalReads(RRBS_B))

correl <- cor.test(meth,cover)
print(correl)

######################
#Statistical analysis#
######################

coeff = 2

design<-model.matrix(~f1+f2)
fit<-lmFit(methLevel(rrbs_B.M),design)
fit <- eBayes(fit)
Bres <- topTable(fit, coef=coeff, n=nrow(methLevel(rrbs_B.M)))

print(head(Bres))

jpeg(paste("hist_pVal_cov",mincov,".jpg",sep=""))
hist(Bres$P.Value)
dev.off()

print(sum(Bres$adj.P.Val>0.05))

Bres <- Bres[sort(as.numeric(rownames(Bres)),index.return=T)$ix,]

B_obs <- cbind(apply(methLevel(rrbs_B.rel)[,f1=="T"],1,mean),
               apply(methLevel(rrbs_B.rel)[,f1=="T"],1,sd),
               apply(methLevel(rrbs_B.rel)[,f1=="C"],1,mean),
               apply(methLevel(rrbs_B.rel)[,f1=="C"],1,sd))


if (sum(rownames(rrbs_B.rel)==rownames(Bres))!=nrow(Bres)){
  stop("Annotation will be faulty!")
}
if (sum(rownames(rrbs_B.rel)==rownames(B_obs))!=nrow(Bres)){
  stop("Annotation will be faulty!")
}

func <- rep("intergenic",length(locs))
func[genic] <- "genic"
func[promoter] <- "promoter"

cn_ori <- colnames(Bres)

Bres_complete <- cbind(as.data.frame(locs)[1:2],gene[,1],gene_symbol[,1],func,Bres,B_obs,
                       gene[,2:ncol(gene)],gene_symbol[,2:ncol(gene_symbol)],cgi_sort)
colnames(Bres_complete) <- c("chr","locations","gene","gene_symbol","func",cn_ori,
                             "mean_T","sd_T","mean_C","sd_C",
                             paste("gene_alt",1:(ncol(gene)-1),sep="_"),
                             paste("gene_symbol_alt",1:(ncol(gene)-1),sep="_"),
                             "cgi_sort")

Bres_complete <- Bres_complete[sort(Bres_complete$P.Value,index.return=T)$ix,]
print(head(Bres_complete))
write.table(Bres_complete,"output_rrbs.txt",sep="\t",col.names = T,row.names = F, quote = F)


###################
# DMR analysis #
###################

dist_max <- 200
rrbs.clust<- clusterSites(object = RRBS_B, groups = f1, perc.samples = 0.75,
                          min.sites = 15, max.dist = dist_max, mc.cores = 20)
clusterSitesToGR(rrbs.clust)
predictedMeth <- predictMeth(object = rrbs.clust, mc.cores = 20)

betaResults <- betaRegression(formula = ~f1+f2, link = "probit", object = predictedMeth, 
                              type = "BR", mc.cores = 20)

vario <- makeVariogram(betaResults)
vario.sm <- smoothVariogram(vario, sill = 0.9)
locCor <- estLocCor(vario.sm)
clusters.rej <- testClusters(locCor,FDR.cluster = 0.1)
head(vario.sm$pValsList[[1]])
clusters.trimmed <- trimClusters(clusters.rej, FDR.loc = 0.05)

DMRs <- findDMRs(clusters.trimmed, max.dist = dist_max, diff.dir = TRUE)

cluster_tr <- GRanges(data.frame(chr=clusters.trimmed$chr,
                                 start=clusters.trimmed$pos,
                                 end=clusters.trimmed$pos))

overlaps_cpg <- findOverlaps(DMRs,cluster_tr)
DMRs@elementMetadata$am_CpG <- as.numeric(table(queryHits(overlaps_cpg)))


DMRs_sorted <- DMRs[sort(DMRs$median.p,index.return=T)$ix,]
DMRs_sorted <- data.frame(DMRs_sorted,stringsAsFactors=F)
DMRs_sorted$seqnames <- as.character(DMRs_sorted$seqnames)

print(head(DMRs_sorted))
print(dim(DMRs_sorted))


# Promoter & genic
DMR_locs <- GRanges(data.frame(chr=DMRs_sorted$seqnames,
                               start=DMRs_sorted$start,
                               end=DMRs_sorted$end))
DMR_locs@ranges@NAMES <- as.character(1:length(DMR_locs))

overlap_p <- subsetByOverlaps(DMR_locs,prom_gr)
overlap_g <- subsetByOverlaps(DMR_locs,genic_gr)

promoter_DMR <- rownames(as.data.frame(DMR_locs))%in%rownames(as.data.frame(overlap_p))
genic_DMR <- rownames(as.data.frame(DMR_locs))%in%rownames(as.data.frame(overlap_g))

genic_DMR[genic_DMR & promoter_DMR] <- FALSE


# Gene and gene symbol annotation
overlaps_DMR_genic <- findOverlaps(DMR_locs,genic_gr)
overlaps_DMR_prom <- findOverlaps(DMR_locs,prom_gr)

all_overlap_DMR <- c(as.character(genic_gr[subjectHits(overlaps_DMR_genic)]$gene),as.character(prom_gr[subjectHits(overlaps_DMR_prom)]$gene))
all_overlap_DMR <- cbind(all_overlap_DMR,c(as.character(genic_gr[subjectHits(overlaps_DMR_genic)]$gene_symbol),as.character(prom_gr[subjectHits(overlaps_DMR_prom)]$gene_symbol)))
all_overlap_DMR <- cbind(all_overlap_DMR,c(DMR_locs[queryHits(overlaps_DMR_genic)]@ranges@NAMES,DMR_locs[queryHits(overlaps_DMR_prom)]@ranges@NAMES))
colnames(all_overlap_DMR) <- c("gene","gene_symbol","rowname")


rn_DMR <- as.matrix(DMR_locs@ranges@NAMES,ncol=1)
colnames(rn_DMR) <- "rowname"

merger_DMR <- merge(rn_DMR,all_overlap_DMR,all=T)
merger_DMR <- merger_DMR[sort(as.numeric(as.character(merger_DMR$rowname)),index.return=T)$ix,]

merger_DMR$gene <- as.character(merger_DMR$gene)
merger_DMR$gene_symbol <- as.character(merger_DMR$gene_symbol)

gene_DMR <- as.character(merger_DMR$gene[!duplicated(merger_DMR$rowname)])
gene_symbol_DMR <- as.character(merger_DMR$gene_symbol[!duplicated(merger_DMR$rowname)])

max_n_DMR <- max(table(merger_DMR$rowname))
print(max_n_DMR)
gene_DMR <- cbind(gene_DMR,matrix(rep(NA,length(gene_DMR)*(max_n_DMR-1)),ncol=max_n_DMR-1))
gene_symbol_DMR <- cbind(gene_symbol_DMR,matrix(rep(NA,length(gene_symbol_DMR)*(max_n_DMR-1)),ncol=max_n_DMR-1))
for (i in 2:max_n_DMR){
  n_times_DMR <- merger_DMR[merger_DMR$rowname%in%names(table(merger_DMR$rowname)[table(merger_DMR$rowname)==i]),]
  if(nrow(n_times_DMR)==0){
    print(paste("no locus with ",i,"-fold annotation",sep=""))
    next
  }
  if (i==2 | nrow(n_times_DMR)==i){
    gene_DMR[DMR_locs@ranges@NAMES%in%unique(n_times_DMR$rowname),2:i] <- n_times_DMR$gene[c(1:((nrow(n_times_DMR)/i)*i))[-c(c(1:(nrow(n_times_DMR)/i))*i-(i-1))]]
    gene_symbol_DMR[DMR_locs@ranges@NAMES%in%unique(n_times_DMR$rowname),2:i] <- n_times_DMR$gene_symbol[c(1:((nrow(n_times_DMR)/i)*i))[-c(c(1:(nrow(n_times_DMR)/i))*i-(i-1))]]
  } else {
    gene_DMR[DMR_locs@ranges@NAMES%in%unique(n_times_DMR$rowname),2:i] <- t(matrix(n_times_DMR$gene[c(1:((nrow(n_times_DMR)/i)*i))[-c(c(1:(nrow(n_times_DMR)/i))*i-(i-1))]],nrow=i-1))
    gene_symbol_DMR[DMR_locs@ranges@NAMES%in%unique(n_times_DMR$rowname),2:i] <- t(matrix(n_times_DMR$gene_symbol[c(1:((nrow(n_times_DMR)/i)*i))[-c(c(1:(nrow(n_times_DMR)/i))*i-(i-1))]],nrow=i-1))
  }
  
}


type_DMR <- rep("intergenic",nrow(gene_DMR))
type_DMR[genic_DMR] <- "genic"
type_DMR[promoter_DMR] <- "promoter"
paste_genenames <- function(char_vec){
  paste(unique(char_vec[!is.na(char_vec)]),collapse=',')
}
DMRs_sorted_annot <- cbind(apply(gene_DMR,1,paste_genenames),apply(gene_symbol_DMR,1,paste_genenames),type_DMR,as.data.frame(DMRs_sorted))
colnames(DMRs_sorted_annot) <- c("genes","gene_symbols","type",colnames(DMRs_sorted))
print(head(DMRs_sorted_annot))

DMRs_sorted_annot$FDR <- p.adjust(DMRs_sorted_annot$median.p, method="BH")

write.table(DMRs_sorted_annot,file="DMR_analysis.txt",col.names=T,row.names = F,quote = F,sep="\t")

# #####################
# # Gene set analysis #
# #####################
# 
# # Entrez Gene id based!
# download.file("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c5_v5p2.rdata","mouse_c5_v5p2.rdata", mode = "wb")
# load("mouse_c5_v5p2.rdata")
# 
# indices <- ids2indices(Mm.c5,gene_entrez)
# 
# r<-romer(methLevel(rrbs_B.M),index=indices,design=design,contrast=2)
# write.table(topRomer(r,alternative="mixed",n=150),file="GSA.txt",col.names = T, row.names = T,sep="\t",quote=F)
# #Note: result not FDR adjusted
# 
# # setwd("C:/Users/Louis/Desktop/")
# # results <- read.table("output_rrbs.txt",header=T,sep="\t")
# # results <- results[!duplicated(results$gene_sort),]
# # gene.vector<-as.numeric(results$P.Value<0.005)
# # names(gene.vector)<-results$gene_sort
# # pwf<-nullp(gene.vector,"mm9","ensGene")
# # GO.wall<-goseq(pwf,"mm9","ensGene")
# # 
# # write.table(GO.wall,file="GSA_analysis_all.txt",col.names=T,row.names = F,quote = F,sep="\t")

############################################################################
## Human Infinium MethylationEPIC analysis and comparison with mouse RRBS ##
############################################################################
#
# Human data: Rodriguez-Paredes M, Bormann F, Raddatz G, et al. (Lyko lab).
# "Methylation profiling identifies two subclasses of squamous cell carcinoma
# related to distinct cells of origin." Nature Communications (2018).
# ArrayExpress E-MTAB-5738: healthy skin (young/old), actinic keratosis, cSCC.
#
# Risk analysis: Hervas-Marin et al., ArrayExpress E-MTAB-8542 (high- vs low-risk cSCC).


## Load additional functions
source("AdditionalFunctions.R")

library(wateRmelon)
library(ChAMPdata)
library(limma)
library(biomaRt)
library(gplots)
library(yaml)

# Directories are read from environment.yml (copy environment_example.yml and fill in the local paths)
if (!file.exists("environment.yml")){
  stop("environment.yml not found, copy environment_example.yml and fill in the local paths")
}
env <- read_yaml("environment.yml")

dir_work <- env$dir_work
setwd(dir_work)
dir.create("Plots", showWarnings = F)
dir.create("Risk_analysis", showWarnings = F)


##################
## Load in data ##
##################

## Load in right biomaRt objects
if (!file.exists("human.Rda")){
  human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host='jan2019.archive.ensembl.org')
  mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",host='jan2019.archive.ensembl.org')
  save(human,file="human.Rda")
  save(mouse,file="mouse.Rda")
} else {
  load("human.Rda")
  load("mouse.Rda")
}


## Load and merge the Infinium data
# can not be loaded together because idat files have different lengths
infdata <- readEPIC("idats/")
infdata_bis <- readEPIC("idats_bis/")
infdata <- infdata[rownames(infdata)%in%rownames(infdata_bis)]
infdata_all <- combo(infdata,infdata_bis)


## Filter on detection p-values
head(pfilter(infdata))


## Remove lines with NA's
print(dim(infdata_all))
infdata_all <- infdata_all[-which(rowSums(is.na(betas(infdata_all)))>0),]
print(dim(infdata_all))


## Load the annotation data in:
samp_annot <- read.table("E-MTAB-5738.sdrf.txt",header=T,sep="\t",stringsAsFactors=F)


# Rename the columns
c_names_new <- NULL
for (t in 1:length(colnames(infdata_all))){
  c_names_new <- c(c_names_new,samp_annot$Source.Name[unlist(strsplit(samp_annot$Assay.Name,"_Grn"))==colnames(infdata_all)[t]])
}
colnames(infdata_all) <- c_names_new


# Get annotation for infinium epic probes
if(!file.exists("probe.features_ens.Rda")){
  data(probe.features.epic)
  probe.features$ens_gene_id <- ""
  addgeneid <- getBM(c("ensembl_gene_id","chromosome_name","start_position","end_position"),values=TRUE, mart = human, uniqueRows=TRUE)
  for (i in 1:nrow(addgeneid)){
    if (i %% 5000 == 0){print(i)}
    probe.features$ens_gene_id[(probe.features$CHR==addgeneid$chromosome_name[i]) & 
                                 ((probe.features$MAPINFO > addgeneid$start_position[i]) & 
                                    (probe.features$MAPINFO < addgeneid$end_position[i]))] <- addgeneid$ensembl_gene_id[i]
  }
  save(probe.features,file="probe.features_ens.Rda")
} else {
  load("probe.features_ens.Rda")
}
annot <- probe.features


##################
## Make designs ##
##################

#   design 1: control vs tumor
#   design 2: Stemcell like vs keratinocyte like
groups <- as.factor(unlist(strsplit(colnames(infdata_all),"-"))[1:length(colnames(infdata_all))*2-1])
id_stemcell <- c("AK-01","AK-03","AK-04","AK-06","AK-07","AK-11","AK-14","SCC-02","SCC-03","SCC-04","SCC-05","SCC-06","SCC-19","SCC-22")
subtypes <- rep("",length(groups))
subtypes[grep("H",groups)] <- "control"
subtypes[colnames(infdata_all)%in%id_stemcell] <- "stemcell"
subtypes[subtypes==""] <- "keratinocyte"
design <- model.matrix(~0+groups)
design2 <- model.matrix(~0+subtypes)


#####################
## Create M-values ##
#####################

# Create M-values object
ct <- 0.01 # Constant value for calculation M-values from methylation degrees
Mval <- log2((betas(infdata_all)+ct)/(1-betas(infdata_all)+ct))


######################
## Methylation plot ##
######################

## Methylation lyko data
cmeans <- colMeans(betas(infdata_all))
state <- as.character(groups)
state[grep("H",as.character(groups))] <- "control"
dat_boxplots_state <- as.data.frame(cbind(cmeans,state))
dat_boxplots_state$cmeans <- as.numeric(as.character(dat_boxplots_state$cmeans)) 
dat_boxplots_state$state <- factor(dat_boxplots_state$state,c("control","AK","SCC"))
dat_boxplots_subtypes <- as.data.frame(cbind(cmeans,subtypes))
dat_boxplots_subtypes$cmeans <- as.numeric(as.character(dat_boxplots_subtypes$cmeans)) 
jpeg("Methylation_lyko_plot.jpg",height=720,width=1440)
par(mfrow=c(1,2),cex.axis=2)
boxplot(cmeans~state,dat_boxplots_state)
boxplot(cmeans~subtypes,dat_boxplots_subtypes)
dev.off()

dat_boxplot_final_state <- dat_boxplots_state
dat_boxplot_final_state$cmeans <- dat_boxplot_final_state$cmeans*100
dat_boxplot_final_subtypes <- dat_boxplots_subtypes
dat_boxplot_final_subtypes$cmeans <- dat_boxplot_final_subtypes$cmeans*100
tiff("Methylation_lyko_plot.tiff",height=8,width=4,units="in",res=150)
par(mfrow=c(2,1))
par(mar=c(2,4,2,2))
boxplot(cmeans~state,names=c("control","AK","cSCC"),data=dat_boxplot_final_state,
        xlab="",ylab="Fraction methylation (%)",cex.axis=0.8)
boxplot(cmeans~subtypes,names=c("control","keratinocyte","stem cell"),data=dat_boxplot_final_subtypes,
        xlab="",ylab="Fraction methylation (%)",cex.axis=0.8)
dev.off()


## Statistical analysis
anova_state <- aov(cmeans~state,dat_boxplots_state)
summary(anova_state)
TukeyHSD(anova_state)

anova_subtypes <- aov(cmeans~subtypes,dat_boxplots_subtypes)
summary(anova_subtypes)
TukeyHSD(anova_subtypes)


################
## Clustering ##
################

# correl <- 1- abs(cov(betas(infdata_all)))
correl <- 1- abs(cov(betas(infdata_all)[apply(betas(infdata_all),1,var)>sort(apply(betas(infdata_all),1,var),decreasing = T)[10001],]))
distance <- as.dist(correl)
hclust_complete <- hclust(distance, method = "complete")

hclust_complete$labels[subtypes=="keratinocyte"] <- paste(hclust_complete$labels[subtypes=="keratinocyte"],", K",sep="")
hclust_complete$labels[subtypes=="stemcell"] <- paste(hclust_complete$labels[subtypes=="stemcell"],", S",sep="")
jpeg("Clustering_lyko.jpg",height=700,width=1050)
par(xpd=TRUE)
par(mar=c(2.5,2.5,2.5,11))
Bclust<-plot(hclust_complete,col.lab="red")
dev.off()


##########################
## Statistical analysis ##
##########################

## normal vs cancer
cont.matrix <- makeContrasts((groupsHSo+groupsHSy)-(groupsAK+groupsSCC),levels=design)
fit <- lmFit(Mval,design)
fit <- contrasts.fit(fit,cont.matrix)
fit <- eBayes(fit)
jpeg("Plots/volcano_controlvscancer.jpg")
volcanoplot(fit)
dev.off()
jpeg("Plots/MA_controlvscancer.jpg")
limma::plotMA(fit)
dev.off()
LIMMAout <- topTable(fit,adjust="BH",number=nrow(infdata_all))
LIMMAout <- LIMMAout[!is.na(LIMMAout$adj.P.Val),]

## normal vs AK
cont.matrix_AK <- makeContrasts((groupsHSo+groupsHSy)-(groupsAK),levels=design)
fit_AK <- lmFit(Mval,design)
fit_AK <- contrasts.fit(fit_AK,cont.matrix_AK)
fit_AK <- eBayes(fit_AK)
jpeg("Plots/volcano_controlvsAK.jpg")
volcanoplot(fit_AK)
dev.off()
jpeg("Plots/MA_controlvsAK.jpg")
limma::plotMA(fit_AK)
dev.off()
LIMMAout_AK <- topTable(fit_AK,adjust="BH",number=nrow(infdata_all))
LIMMAout_AK <- LIMMAout_AK[!is.na(LIMMAout_AK$adj.P.Val),]

## normal vs SCC
cont.matrix_SCC <- makeContrasts((groupsHSo+groupsHSy)-(groupsSCC),levels=design)
fit_SCC <- lmFit(Mval,design)
fit_SCC <- contrasts.fit(fit_SCC,cont.matrix_SCC)
fit_SCC <- eBayes(fit_SCC)
jpeg("Plots/volcano_controlvsSCC.jpg")
volcanoplot(fit_SCC)
dev.off()
jpeg("Plots/MA_controlvsSCC.jpg")
limma::plotMA(fit_SCC)
dev.off()
LIMMAout_SCC <- topTable(fit_SCC,adjust="BH",number=nrow(infdata_all))
LIMMAout_SCC <- LIMMAout_SCC[!is.na(LIMMAout_SCC$adj.P.Val),]

## AK vs SCC
cont.matrix_AKSCC <- makeContrasts((groupsAK)-(groupsSCC),levels=design)
fit_AKSCC <- lmFit(Mval,design)
fit_AKSCC <- contrasts.fit(fit_AKSCC,cont.matrix_AKSCC)
fit_AKSCC <- eBayes(fit_AKSCC)
jpeg("Plots/volcano_AKvsSCC.jpg")
volcanoplot(fit_AKSCC)
dev.off()
jpeg("Plots/MA_AKvsSCC.jpg")
limma::plotMA(fit_AKSCC)
dev.off()
LIMMAout_AKSCC <- topTable(fit_AKSCC,adjust="BH",number=nrow(infdata_all))
LIMMAout_AKSCC <- LIMMAout_AKSCC[!is.na(LIMMAout_AKSCC$adj.P.Val),]

## Control vs stemcell-like
cont.matrix_Stem <- makeContrasts((subtypescontrol)-(subtypesstemcell),levels=design2)
fit_Stem <- lmFit(Mval,design2)
fit_Stem <- contrasts.fit(fit_Stem,cont.matrix_Stem)
fit_Stem <- eBayes(fit_Stem)
jpeg("Plots/volcano_controlvsstemcell.jpg")
volcanoplot(fit_Stem)
dev.off()
jpeg("Plots/MA_controlvsstemcell.jpg")
limma::plotMA(fit_Stem)
dev.off()
LIMMAout_Stem <- topTable(fit_Stem,adjust="BH",number=nrow(infdata_all))
LIMMAout_Stem <- LIMMAout_Stem[!is.na(LIMMAout_Stem$adj.P.Val),]

## Control vs keratinocyte
cont.matrix_Ker <- makeContrasts((subtypescontrol)-(subtypeskeratinocyte),levels=design2)
fit_Ker <- lmFit(Mval,design2)
fit_Ker <- contrasts.fit(fit_Ker,cont.matrix_Ker)
fit_Ker <- eBayes(fit_Ker)
jpeg("Plots/volcano_controlvskeratinocyte.jpg")
volcanoplot(fit_Ker)
dev.off()
jpeg("Plots/MA_controlvskeratinocyte.jpg")
limma::plotMA(fit_Ker)
dev.off()
LIMMAout_Ker <- topTable(fit_Ker,adjust="BH",number=nrow(infdata_all))
LIMMAout_Ker <- LIMMAout_Ker[!is.na(LIMMAout_Ker$adj.P.Val),]

## Stemcell-like vs keratinocyte
cont.matrix_StemKer <- makeContrasts((subtypesstemcell)-(subtypeskeratinocyte),levels=design2)
fit_StemKer <- lmFit(Mval,design2)
fit_StemKer <- contrasts.fit(fit_StemKer,cont.matrix_StemKer)
fit_StemKer <- eBayes(fit_StemKer)
jpeg("Plots/volcano_stemcellvskeratinocyte.jpg")
volcanoplot(fit_StemKer)
dev.off()
jpeg("Plots/MA_stemcellvskeratinocyte.jpg")
limma::plotMA(fit_StemKer)
dev.off()
LIMMAout_StemKer <- topTable(fit_StemKer,adjust="BH",number=nrow(infdata_all))
LIMMAout_StemKer <- LIMMAout_StemKer[!is.na(LIMMAout_StemKer$adj.P.Val),]


#############
## Results ##
#############

## Summary results
head(LIMMAout)
diff_betas <- betas(infdata)
dim(LIMMAout)
sum(LIMMAout$adj.P.Val<0.05)

head(LIMMAout_AK)
sum(LIMMAout_AK$adj.P.Val<0.05)

head(LIMMAout_SCC)
sum(LIMMAout_SCC$adj.P.Val<0.05)

head(LIMMAout_AKSCC)
sum(LIMMAout_AKSCC$adj.P.Val<0.05)

head(LIMMAout_Stem)
sum(LIMMAout_Stem$adj.P.Val<0.05)

head(LIMMAout_Ker)
sum(LIMMAout_Ker$adj.P.Val<0.05)

head(LIMMAout_StemKer)
sum(LIMMAout_StemKer$adj.P.Val<0.05)

dat_chisq <- rbind(c(sum(LIMMAout_StemKer$adj.P.Val<0.05),nrow(LIMMAout_StemKer)-sum(LIMMAout_StemKer$adj.P.Val<0.05)),
                   c(sum(LIMMAout$adj.P.Val<0.05),nrow(LIMMAout)-sum(LIMMAout$adj.P.Val<0.05)))
dat_chisq/rowSums(dat_chisq)
chisq.test(dat_chisq)


## Get annotation for resultsobject
annot <- annot[rownames(annot)%in%rownames(LIMMAout),]

LIMMAout <- LIMMAout[sort(rownames(LIMMAout),index.return=T)$ix,]
index <- sort(LIMMAout$P.Value,index.return=T)$ix
LIMMAout <- LIMMAout[sort(LIMMAout$P.Value,index.return=T)$ix,]

annot_LIMMAout <- annot[index,]
paste(length(unique(annot_LIMMAout$gene[LIMMAout$adj.P.Val<0.05]))-1,"out of",length(unique(annot_LIMMAout$gene))-1)

LIMMAout_AK <- LIMMAout_AK[sort(rownames(LIMMAout_AK),index.return=T)$ix,]
index_AK <- sort(LIMMAout_AK$P.Value,index.return=T)$ix
LIMMAout_AK <- LIMMAout_AK[sort(LIMMAout_AK$P.Value,index.return=T)$ix,]

annot_LIMMAout_AK <- annot[index_AK,]

LIMMAout_SCC <- LIMMAout_SCC[sort(rownames(LIMMAout_SCC),index.return=T)$ix,]
index_SCC <- sort(LIMMAout_SCC$P.Value,index.return=T)$ix
LIMMAout_SCC <- LIMMAout_SCC[sort(LIMMAout_SCC$P.Value,index.return=T)$ix,]

annot_LIMMAout_SCC <- annot[index_SCC,]

LIMMAout_AKSCC <- LIMMAout_AKSCC[sort(rownames(LIMMAout_AKSCC),index.return=T)$ix,]
index_AKSCC <- sort(LIMMAout_AKSCC$P.Value,index.return=T)$ix
LIMMAout_AKSCC <- LIMMAout_AKSCC[sort(LIMMAout_AKSCC$P.Value,index.return=T)$ix,]

annot_LIMMAout_AKSCC <- annot[index_AKSCC,]

LIMMAout_Stem <- LIMMAout_Stem[sort(rownames(LIMMAout_Stem),index.return=T)$ix,]
index_Stem <- sort(LIMMAout_Stem$P.Value,index.return=T)$ix
LIMMAout_Stem <- LIMMAout_Stem[sort(LIMMAout_Stem$P.Value,index.return=T)$ix,]

annot_LIMMAout_Stem <- annot[index_Stem,]

LIMMAout_Ker <- LIMMAout_Ker[sort(rownames(LIMMAout_Ker),index.return=T)$ix,]
index_Ker <- sort(LIMMAout_Ker$P.Value,index.return=T)$ix
LIMMAout_Ker <- LIMMAout_Ker[sort(LIMMAout_Ker$P.Value,index.return=T)$ix,]

annot_LIMMAout_Ker <- annot[index_Ker,]

LIMMAout_StemKer <- LIMMAout_StemKer[sort(rownames(LIMMAout_StemKer),index.return=T)$ix,]
index_StemKer <- sort(LIMMAout_StemKer$P.Value,index.return=T)$ix
LIMMAout_StemKer <- LIMMAout_StemKer[sort(LIMMAout_StemKer$P.Value,index.return=T)$ix,]

annot_LIMMAout_StemKer <- annot[index_StemKer,]

dim(LIMMAout)
sum(rownames(LIMMAout)==rownames(annot_LIMMAout))
dim(LIMMAout_AK)
sum(rownames(LIMMAout_AK)==rownames(annot_LIMMAout_AK))
dim(LIMMAout_SCC)
sum(rownames(LIMMAout_SCC)==rownames(annot_LIMMAout_SCC))
dim(LIMMAout_Ker)
sum(rownames(LIMMAout_Ker)==rownames(annot_LIMMAout_Ker))
dim(LIMMAout_Stem)
sum(rownames(LIMMAout_Stem)==rownames(annot_LIMMAout_Stem))


## Make object for Lyko probes in DM genes of analysis Kevin
lyko <- cbind(annot_LIMMAout,
              rowMeans(betas(infdata_all)[,!grepl("H",groups)])[index],
              apply(betas(infdata_all)[,!grepl("H",groups)],1,sd)[index],
              rowMeans(betas(infdata_all)[,grepl("H",groups)])[index],
              apply(betas(infdata_all)[,grepl("H",groups)],1,sd)[index],
              LIMMAout)
colnames(lyko) <- c(colnames(annot_LIMMAout),"Mean_T","Sd_T","Mean_C","Sd_C",colnames(LIMMAout))

lyko_Ker <- cbind(annot_LIMMAout_Ker,
                  rowMeans(betas(infdata_all)[,subtypes=="keratinocyte"])[index_Ker],
                  apply(betas(infdata_all)[,subtypes=="keratinocyte"],1,sd)[index_Ker],
                  rowMeans(betas(infdata_all)[,grepl("H",groups)])[index_Ker],
                  apply(betas(infdata_all)[,grepl("H",groups)],1,sd)[index_Ker],
                  LIMMAout_Ker)
colnames(lyko_Ker) <- c(colnames(annot_LIMMAout_Ker),"Mean_T","Sd_T","Mean_C","Sd_C",colnames(LIMMAout_Ker))
write.table(lyko_Ker,file="results_lyko_Ker.txt",col.names=T,row.names = T,sep="\t",quote=F)

lyko_Stem <- cbind(annot_LIMMAout_Stem,
                   rowMeans(betas(infdata_all)[,subtypes=="stemcell"])[index_Stem],
                   apply(betas(infdata_all)[,subtypes=="stemcell"],1,sd)[index_Stem],
                   rowMeans(betas(infdata_all)[,grepl("H",groups)])[index_Stem],
                   apply(betas(infdata_all)[,grepl("H",groups)],1,sd)[index_Stem],
                   LIMMAout_Stem)
colnames(lyko_Stem) <- c(colnames(annot_LIMMAout_Stem),"Mean_T","Sd_T","Mean_C","Sd_C",colnames(LIMMAout_Stem))
write.table(lyko_Stem,file="results_lyko_Stem.txt",col.names=T,row.names = T,sep="\t",quote=F)

# Numbers for paper
tot_genes <- length(unique(lyko$ens_gene_id)[unique(lyko$ens_gene_id)!=""])
table_DM <- NULL
table_DM <- rbind(table_DM,c(sum(lyko$adj.P.Val<0.05),
                             sum(lyko_Ker$adj.P.Val<0.05),
                             sum(lyko_Stem$adj.P.Val<0.05)))
table_DM <- rbind(table_DM,table_DM[1,]/nrow(lyko))
table_DM <- rbind(table_DM,c(length(unique(lyko$ens_gene_id[lyko$adj.P.Val<0.05])),
                             length(unique(lyko_Ker$ens_gene_id[lyko_Ker$adj.P.Val<0.05])),
                             length(unique(lyko_Stem$ens_gene_id[lyko_Stem$adj.P.Val<0.05]))))
table_DM <- rbind(table_DM,table_DM[3,]/tot_genes)
rownames(table_DM) <- c("DM probes","DM probes rel",
                        "Genes with a DM probe","Genes with a DM probe rel")
colnames(table_DM) <- c("tumor","Ker","Stem")
t(table_DM)


#######################
## Load RRBS results ##
#######################

## Load in data
output_rrbs <- read.table("output_rrbs.txt",header=T,sep="\t",stringsAsFactors = F)
betas_rrbs <- read.table("Betas_rrbs.txt",header=T,sep="\t")
annot_rrbs <- read.table("Annot_rrbs.txt",header=T,sep="\t")

output_rrbs$gene_symbol_all <- ""
output_rrbs$gene_all <- ""
for (i in 1:nrow(output_rrbs)){
  gene_symbols <- output_rrbs[i,grepl("gene_symbol$",colnames(output_rrbs)) | grepl("gene_symbol_alt",colnames(output_rrbs))]
  gene_symbols <- sort(unique(gene_symbols[!is.na(gene_symbols)]))
  gene <- output_rrbs[i,grepl("gene$",colnames(output_rrbs)) | grepl("gene_alt",colnames(output_rrbs))]
  gene <- sort(unique(gene[!is.na(gene)]))
  if(length(gene_symbols)==0){
    output_rrbs$gene_symbol_all[i] <- NA
  } else {
    output_rrbs$gene_symbol_all[i] <- paste(gene_symbols,collapse=",")
  }
  if(length(gene)==0){
    output_rrbs$gene_all[i] <- NA
  } else {
    output_rrbs$gene_all[i] <- paste(gene,collapse=",")
  }
}
write.table(output_rrbs,file="output_rrbs_suppl.txt",col.names=T,row.names = F,quote = F,sep="\t")

######################
## Extract DE genes ##
######################


## Extract DE genes
# lyko: Control vs cancer
if (!file.exists("all_genes_lyko.Rda")){
  all_genes_lyko <- lyko$ens_gene_id
  all_genes_lyko <- sort(toupper(unique(all_genes_lyko)))[sort(toupper(unique(all_genes_lyko)))!=""]
  all_genes_lyko <- data.frame(gene=all_genes_lyko,stringsAsFactors = F)
  all_genes_lyko$aveLogFC <- getAveLogFC(all_genes_lyko$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("logFC",colnames(lyko)))
  all_genes_lyko$gene_symbol <- getGeneSymbol(all_genes_lyko$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("gene$",colnames(lyko)))
  save(all_genes_lyko,file="all_genes_lyko.Rda")
  
  DE_genes_lyko <- lyko$ens_gene_id[lyko$adj.P.Val<0.05]
  DE_genes_lyko <- sort(toupper(unique(DE_genes_lyko)))[sort(toupper(unique(DE_genes_lyko)))!=""]
  DE_genes_lyko <- data.frame(gene=DE_genes_lyko,stringsAsFactors = F)
  DE_genes_lyko$aveLogFC <- getAveLogFC(DE_genes_lyko$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("logFC",colnames(lyko)))
  DE_genes_lyko$gene_symbol <- getGeneSymbol(DE_genes_lyko$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("gene$",colnames(lyko)))
  save(DE_genes_lyko,file="DE_genes_lyko.Rda")
  
  DE_genes_lyko_prom <- lyko$ens_gene_id[(lyko$adj.P.Val<0.05) & grepl("TSS",lyko$feature)]
  DE_genes_lyko_prom <- sort(toupper(unique(DE_genes_lyko_prom)))[sort(toupper(unique(DE_genes_lyko_prom)))!=""]
  DE_genes_lyko_prom <- data.frame(gene=DE_genes_lyko_prom,stringsAsFactors = F)
  DE_genes_lyko_prom$aveLogFC <- getAveLogFC(DE_genes_lyko_prom$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("logFC",colnames(lyko)),
                                             grepl("TSS",lyko$feature))
  DE_genes_lyko_prom$gene_symbol <- getGeneSymbol(DE_genes_lyko_prom$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("gene$",colnames(lyko)),
                                                  grepl("TSS",lyko$feature))
  save(DE_genes_lyko_prom,file="DE_genes_lyko_prom.Rda")
  
  # lyko: Control vs keratinocyte-like
  DE_genes_lyko_Ker <- lyko_Ker$ens_gene_id[lyko_Ker$adj.P.Val<0.05]
  DE_genes_lyko_Ker <- sort(toupper(unique(DE_genes_lyko_Ker)))[sort(toupper(unique(DE_genes_lyko_Ker)))!=""]
  DE_genes_lyko_Ker <- data.frame(gene=DE_genes_lyko_Ker,stringsAsFactors = F)
  DE_genes_lyko_Ker$aveLogFC <- getAveLogFC(DE_genes_lyko_Ker$gene,lyko_Ker,grep("ens_gene_id",colnames(lyko_Ker)),grep("logFC",colnames(lyko_Ker)))
  DE_genes_lyko_Ker$gene_symbol <- getGeneSymbol(DE_genes_lyko_Ker$gene,lyko_Ker,grep("ens_gene_id",colnames(lyko_Ker)),grep("gene$",colnames(lyko_Ker)))
  save(DE_genes_lyko_Ker,file="DE_genes_lyko_Ker.Rda")
  
  DE_genes_lyko_Ker_prom <- lyko_Ker$ens_gene_id[(lyko_Ker$adj.P.Val<0.05) & grepl("TSS",lyko_Ker$feature)]
  DE_genes_lyko_Ker_prom <- sort(toupper(unique(DE_genes_lyko_Ker_prom)))[sort(toupper(unique(DE_genes_lyko_Ker_prom)))!=""]
  DE_genes_lyko_Ker_prom <- data.frame(gene=DE_genes_lyko_Ker_prom,stringsAsFactors = F)
  DE_genes_lyko_Ker_prom$aveLogFC <- getAveLogFC(DE_genes_lyko_Ker_prom$gene,lyko_Ker,grep("ens_gene_id",colnames(lyko_Ker)),grep("logFC",colnames(lyko_Ker)),
                                                 grepl("TSS",lyko_Ker$feature))
  DE_genes_lyko_Ker_prom$gene_symbol <- getGeneSymbol(DE_genes_lyko_Ker_prom$gene,lyko_Ker,grep("ens_gene_id",colnames(lyko_Ker)),grep("gene$",colnames(lyko_Ker)),
                                                  grepl("TSS",lyko_Ker$feature))
  save(DE_genes_lyko_Ker_prom,file="DE_genes_lyko_Ker_prom.Rda")
  
  # lyko: Control vs stemcell-like
  DE_genes_lyko_Stem <- lyko_Stem$ens_gene_id[lyko_Stem$adj.P.Val<0.05]
  DE_genes_lyko_Stem <- sort(toupper(unique(DE_genes_lyko_Stem)))[sort(toupper(unique(DE_genes_lyko_Stem)))!=""]
  DE_genes_lyko_Stem <- data.frame(gene=DE_genes_lyko_Stem,stringsAsFactors = F)
  DE_genes_lyko_Stem$aveLogFC <- getAveLogFC(DE_genes_lyko_Stem$gene,lyko_Stem,grep("ens_gene_id",colnames(lyko_Stem)),grep("logFC",colnames(lyko_Stem)))
  DE_genes_lyko_Stem$gene_symbol <- getGeneSymbol(DE_genes_lyko_Stem$gene,lyko_Stem,grep("ens_gene_id",colnames(lyko_Stem)),grep("gene$",colnames(lyko_Stem)))
  save(DE_genes_lyko_Stem,file="DE_genes_lyko_Stem.Rda")
  
  DE_genes_lyko_Stem_prom <- lyko_Stem$ens_gene_id[(lyko_Stem$adj.P.Val<0.05) & grepl("TSS",lyko_Stem$feature)]
  DE_genes_lyko_Stem_prom <- sort(toupper(unique(DE_genes_lyko_Stem_prom)))[sort(toupper(unique(DE_genes_lyko_Stem_prom)))!=""]
  DE_genes_lyko_Stem_prom <- data.frame(gene=DE_genes_lyko_Stem_prom,stringsAsFactors = F)
  DE_genes_lyko_Stem_prom$aveLogFC <- getAveLogFC(DE_genes_lyko_Stem_prom$gene,lyko_Stem,grep("ens_gene_id",colnames(lyko_Stem)),grep("logFC",colnames(lyko_Stem)),
                                                 grepl("TSS",lyko_Stem$feature))
  DE_genes_lyko_Stem_prom$gene_symbol <- getGeneSymbol(DE_genes_lyko_Stem_prom$gene,lyko_Stem,grep("ens_gene_id",colnames(lyko_Stem)),grep("gene$",colnames(lyko_Stem)),
                                                  grepl("TSS",lyko_Stem$feature))
  save(DE_genes_lyko_Stem_prom,file="DE_genes_lyko_Stem_prom.Rda")
  
  # Kevin: Control vs cancer
  # DE_genes_rrbs_mat <- unlist(output_rrbs[output_rrbs$adj.P.Val<0.10,
  #                                   grepl("gene",colnames(output_rrbs)) & !(grepl("gene_symbol",colnames(output_rrbs)))])
  # DE_genes_rrbs <- DE_genes_rrbs[!is.na(DE_genes_rrbs)]
  DE_genes_rrbs <- output_rrbs$gene[output_rrbs$adj.P.Val<0.10]
  DE_genes_rrbs <- sort(toupper(unique(DE_genes_rrbs)))[sort(toupper(unique(DE_genes_rrbs)))!=""]
  DE_genes_rrbs <- data.frame(gene=DE_genes_rrbs,stringsAsFactors = F)
  DE_genes_rrbs$aveLogFC <- getAveLogFC(DE_genes_rrbs$gene,output_rrbs,grep("gene$",colnames(output_rrbs)),grep("logFC",colnames(output_rrbs)))
  DE_genes_rrbs$gene_symbol <- getGeneSymbol(DE_genes_rrbs$gene,output_rrbs,grep("gene$",colnames(output_rrbs)),grep("gene_symbol_all",colnames(output_rrbs)))
  save(DE_genes_rrbs,file="DE_genes_rrbs.Rda")
  
  # DE_genes_rrbs_mat <- unlist(output_rrbs[output_rrbs$adj.P.Val<0.10,
  #                                   grepl("gene",colnames(output_rrbs)) & !(grepl("gene_symbol",colnames(output_rrbs)))])
  # DE_genes_rrbs <- DE_genes_rrbs[!is.na(DE_genes_rrbs)]
  DE_genes_rrbs_prom <- output_rrbs$gene[output_rrbs$adj.P.Val<0.10 & output_rrbs$func=="promoter"]
  DE_genes_rrbs_prom <- sort(toupper(unique(DE_genes_rrbs_prom)))[sort(toupper(unique(DE_genes_rrbs_prom)))!=""]
  DE_genes_rrbs_prom <- data.frame(gene=DE_genes_rrbs_prom,stringsAsFactors = F)
  DE_genes_rrbs_prom$aveLogFC <- getAveLogFC(DE_genes_rrbs_prom$gene,output_rrbs,grep("gene$",colnames(output_rrbs)),grep("logFC",colnames(output_rrbs)),
                                             output_rrbs$func=="promoter")
  DE_genes_rrbs_prom$gene_symbol <- getGeneSymbol(DE_genes_rrbs_prom$gene,output_rrbs,grep("gene$",colnames(output_rrbs)),grep("gene_symbol_all",colnames(output_rrbs)),
                                                  output_rrbs$func=="promoter")
  save(DE_genes_rrbs_prom,file="DE_genes_rrbs_prom.Rda")
} else {
  load("all_genes_lyko.Rda")
  load("DE_genes_lyko.Rda")
  load("DE_genes_lyko_prom.Rda")
  
  load("DE_genes_lyko_Stem.Rda")
  load("DE_genes_lyko_Stem_prom.Rda")
  
  load("DE_genes_lyko_Ker.Rda")
  load("DE_genes_lyko_Ker_prom.Rda")
  
  load("DE_genes_rrbs.Rda")
  load("DE_genes_rrbs_prom.Rda")
}

## Overlap
length(all_genes_lyko$gene)
length(DE_genes_lyko$gene)
sum(DE_genes_lyko$gene%in%all_genes_lyko$gene)/length(all_genes_lyko$gene)

length(DE_genes_rrbs$gene)


#####################
## Mouse orthologs ##
#####################

## Make ortholog list
if (!file.exists("orth.mouse.Rda")){
  attributes = c("ensembl_gene_id","mmusculus_homolog_ensembl_gene","mmusculus_homolog_perc_id_r1","mmusculus_homolog_orthology_confidence")
  orth.mouse_1 = getBM(attributes=attributes,filters="with_mmusculus_homolog",values=TRUE, mart = human, uniqueRows=TRUE)
  orth.mouse_2 = getBM(c("ensembl_gene_id","hgnc_symbol"),filters="with_mmusculus_homolog",values=TRUE, mart = human, uniqueRows=TRUE)
  mgi = getBM(c("ensembl_gene_id","mgi_symbol"),values=TRUE, mart = mouse, uniqueRows=TRUE)
  orth.mouse <- merge(orth.mouse_1,orth.mouse_2)
  save(orth.mouse,file="orth.mouse.Rda")
  save(mgi,file="mgi.Rda")
} else {
  load("orth.mouse.Rda")
  load("mgi.Rda")
}
dim(orth.mouse)
sum(duplicated(orth.mouse$ensembl_gene_id))
orth.mouse_new <- orth.mouse[!duplicated(orth.mouse$ensembl_gene_id) & !base::duplicated(orth.mouse$ensembl_gene_id,fromLast=T),]
rest <- unique(orth.mouse$ensembl_gene_id)[!(unique(orth.mouse$ensembl_gene_id)%in%orth.mouse_new$ensembl_gene_id)]
teller <- 0
teller2 <- 0
for(i in 1:length(rest)){
  print(i)
  sub <- orth.mouse[orth.mouse$ensembl_gene_id==rest[i],]
  mgi_sub <- mgi[mgi$ensembl_gene_id%in%sub$mmusculus_homolog_ensembl_gene,]
  if(sum(sub$mmusculus_homolog_orthology_confidence,na.rm = T)==0 ){
    best <- sub[sub$mmusculus_homolog_perc_id_r1==max(sub$mmusculus_homolog_perc_id_r1),]
  } else if(sum(sub$mmusculus_homolog_orthology_confidence)==1){
    best <- sub[sub$mmusculus_homolog_orthology_confidence==1,]
  } else if(sum(sub$mmusculus_homolog_perc_id_r1==max(sub$mmusculus_homolog_perc_id_r1))==1){
    best <- sub[sub$mmusculus_homolog_perc_id_r1==max(sub$mmusculus_homolog_perc_id_r1),]
  } else {
    sub <- sub[sub$mmusculus_homolog_perc_id_r1==max(sub$mmusculus_homolog_perc_id_r1),]
    best <- sub[sub$hgnc_symbol%in%toupper(mgi_sub$mgi_symbol),]
  }
  print(best)
  orth.mouse_new <- rbind(orth.mouse_new,best)
}

DE_genes_rrbs <- getOrthologues(DE_genes_rrbs)
DE_genes_rrbs_prom <- getOrthologues(DE_genes_rrbs_prom)


####################################
## Make overlap between both sets ##
####################################

  
table_overlaps <- matrix(rep(NA,6*3),ncol=3)
colnames(table_overlaps) <- c("Cancer","keratinocyte-like","Stemcell-like")
rownames(table_overlaps) <- c("All N with ortholog","All N present in both","All N same direction",
                              "Prom N with ortholog","Prom N present in both","Prom N same direction")

table_overlaps[1,] <- sum(DE_genes_rrbs$human_ortholog!="")
table_overlaps[4,] <- sum(DE_genes_rrbs_prom$human_ortholog!="")

table_overlaps_GS <- matrix(rep(NA,6*3),ncol=3)
colnames(table_overlaps_GS) <- c("Cancer","keratinocyte-like","Stemcell-like")
rownames(table_overlaps_GS) <- c("All N with ortholog","All N present in both","All N same direction",
                              "Prom N with ortholog","Prom N present in both","Prom N same direction")

table_overlaps_GS[1,] <- length(unique(unlist(strsplit(DE_genes_rrbs$gene_symbol,","))))
table_overlaps_GS[4,] <- length(unique(unlist(strsplit(DE_genes_rrbs_prom$gene_symbol,","))))

## Control vs cancer
# all
Control_vs_cancer <- overlapSets(DE_genes_lyko,DE_genes_rrbs,lyko,output_rrbs)
table_overlaps[2,1] <- nrow(Control_vs_cancer)
table_overlaps[3,1] <- sum(sign(Control_vs_cancer$lyko_aveLogFC)!=sign(Control_vs_cancer$roth_aveLogFC))
Control_vs_cancer_GS <- overlapSetsGS(DE_genes_lyko,DE_genes_rrbs,lyko,output_rrbs)
table_overlaps_GS[2,1] <- nrow(Control_vs_cancer_GS)
table_overlaps_GS[3,1] <- sum(sign(Control_vs_cancer_GS$lyko_aveLogFC)!=sign(Control_vs_cancer_GS$roth_aveLogFC))

# Prom 
Control_vs_cancer_prom <- overlapSets(DE_genes_lyko_prom,DE_genes_rrbs_prom,lyko,output_rrbs)
table_overlaps[5,1] <- nrow(Control_vs_cancer_prom)
table_overlaps[6,1] <- sum(sign(Control_vs_cancer_prom$lyko_aveLogFC)!=sign(Control_vs_cancer_prom$roth_aveLogFC))
Control_vs_cancer_prom_GS <- overlapSetsGS(DE_genes_lyko_prom,DE_genes_rrbs_prom,lyko,output_rrbs)
table_overlaps_GS[5,1] <- nrow(Control_vs_cancer_prom_GS)
table_overlaps_GS[6,1] <- sum(sign(Control_vs_cancer_prom_GS$lyko_aveLogFC)!=sign(Control_vs_cancer_prom_GS$roth_aveLogFC))

## Control vs keratinocyte-like
# all 
Control_vs_Ker <- overlapSets(DE_genes_lyko_Ker,DE_genes_rrbs,lyko_Ker,output_rrbs)
table_overlaps[2,2] <- nrow(Control_vs_Ker)
table_overlaps[3,2] <- sum(sign(Control_vs_Ker$lyko_aveLogFC)!=sign(Control_vs_Ker$roth_aveLogFC))
Control_vs_Ker_GS <- overlapSetsGS(DE_genes_lyko_Ker,DE_genes_rrbs,lyko_Ker,output_rrbs)
table_overlaps_GS[2,2] <- nrow(Control_vs_Ker_GS)
table_overlaps_GS[3,2] <- sum(sign(Control_vs_Ker_GS$lyko_aveLogFC)!=sign(Control_vs_Ker_GS$roth_aveLogFC))

# Prom 
Control_vs_Ker <- overlapSets(DE_genes_lyko_Ker_prom,DE_genes_rrbs_prom,lyko_Ker,output_rrbs)
table_overlaps[5,2] <- nrow(Control_vs_Ker)
table_overlaps[6,2] <- sum(sign(Control_vs_Ker$lyko_aveLogFC)!=sign(Control_vs_Ker$roth_aveLogFC))
Control_vs_Ker_GS <- overlapSetsGS(DE_genes_lyko_Ker_prom,DE_genes_rrbs_prom,lyko_Ker,output_rrbs)
table_overlaps_GS[5,2] <- nrow(Control_vs_Ker_GS)
table_overlaps_GS[6,2] <- sum(sign(Control_vs_Ker_GS$lyko_aveLogFC)!=sign(Control_vs_Ker_GS$roth_aveLogFC))

## Control vs stemcell-like
# all 
Control_vs_Stem <- overlapSets(DE_genes_lyko_Stem,DE_genes_rrbs,lyko_Stem,output_rrbs)
table_overlaps[2,3] <- nrow(Control_vs_Stem)
table_overlaps[3,3] <- sum(sign(Control_vs_Stem$lyko_aveLogFC)!=sign(Control_vs_Stem$roth_aveLogFC))
Control_vs_Stem_GS <- overlapSetsGS(DE_genes_lyko_Stem,DE_genes_rrbs,lyko_Stem,output_rrbs)
table_overlaps_GS[2,3] <- nrow(Control_vs_Stem_GS)
table_overlaps_GS[3,3] <- sum(sign(Control_vs_Stem_GS$lyko_aveLogFC)!=sign(Control_vs_Stem_GS$roth_aveLogFC))

# Prom 
Control_vs_Stem_prom <- overlapSets(DE_genes_lyko_Stem_prom,DE_genes_rrbs_prom,lyko_Stem,output_rrbs)
table_overlaps[5,3] <- nrow(Control_vs_Stem_prom)
table_overlaps[6,3] <- sum(sign(Control_vs_Stem_prom$lyko_aveLogFC)!=sign(Control_vs_Stem_prom$roth_aveLogFC))
Control_vs_Stem_prom_GS <- overlapSetsGS(DE_genes_lyko_Stem_prom,DE_genes_rrbs_prom,lyko_Stem,output_rrbs)
table_overlaps_GS[5,3] <- nrow(Control_vs_Stem_prom_GS)
table_overlaps_GS[6,3] <- sum(sign(Control_vs_Stem_prom_GS$lyko_aveLogFC)!=sign(Control_vs_Stem_prom_GS$roth_aveLogFC))



# Statistics

c(table_overlaps[2,]/table_overlaps[1,])/table_DM[4,]
pvals <- c(chisq.test(rbind(c(table_overlaps[2,1],table_overlaps[1,1]-table_overlaps[2,1]),
                            c(table_DM[3,1],tot_genes-table_DM[3,1])))$p.value,
           chisq.test(rbind(c(table_overlaps[2,2],table_overlaps[1,2]-table_overlaps[2,2]),
                            c(table_DM[3,2],tot_genes-table_DM[3,2])))$p.value,
           chisq.test(rbind(c(table_overlaps[2,3],table_overlaps[1,3]-table_overlaps[2,3]),
                            c(table_DM[3,3],tot_genes-table_DM[3,3])))$p.value)


c(table_overlaps[3,]/table_overlaps[2,])/table_DM[4,]
pvals <- c(binom.test(c(table_overlaps[3,1],table_overlaps[2,1]-table_overlaps[3,1]),c(0.5,0.5))$p.value,
           binom.test(c(table_overlaps[3,2],table_overlaps[2,2]-table_overlaps[3,2]),c(0.5,0.5))$p.value,
           binom.test(c(table_overlaps[3,3],table_overlaps[2,3]-table_overlaps[3,3]),c(0.5,0.5))$p.value)
           
chisq.test(cbind(c(table_overlaps[3,2],table_overlaps[2,2]-table_overlaps[3,2]),
                 c(table_overlaps[3,3],table_overlaps[2,3]-table_overlaps[3,3])))

###################################### 
## Create Genes of Interest GOI-set ##
######################################

## Extract DE genes
if (!file.exists("DE_genes_rrbs_top.Rda")){
  # Kevin: Control vs cancer
  # All
  # DE_genes_rrbs_top_mat <- unlist(output_rrbs[(output_rrbs$adj.P.Val<0.10) & (abs(output_rrbs$mean_T-output_rrbs$mean_C)>0.30),
  #                                   grepl("gene",colnames(output_rrbs)) & !(grepl("gene_symbol",colnames(output_rrbs)))])
  # DE_genes_rrbs_top <- DE_genes_rrbs_top_mat[!is.na(DE_genes_rrbs_top_mat)]
  DE_genes_rrbs_top <- output_rrbs$gene[(output_rrbs$adj.P.Val<0.10) & (abs(output_rrbs$mean_T-output_rrbs$mean_C)>0.30)]
  DE_genes_rrbs_top <- sort(toupper(unique(DE_genes_rrbs_top)))[sort(toupper(unique(DE_genes_rrbs_top)))!=""]
  DE_genes_rrbs_top <- data.frame(gene=DE_genes_rrbs_top,stringsAsFactors = F)
  DE_genes_rrbs_top$aveLogFC <- getAveLogFC(DE_genes_rrbs_top$gene,output_rrbs,grep("gene$",colnames(output_rrbs)),grep("logFC",colnames(output_rrbs)))
  DE_genes_rrbs_top$gene_symbol <- getGeneSymbol(DE_genes_rrbs_top$gene,output_rrbs,grep("gene$",colnames(output_rrbs)),grep("gene_symbol_all",colnames(output_rrbs)))
  DE_genes_rrbs_top <- getOrthologues(DE_genes_rrbs_top)
  save(DE_genes_rrbs_top,file="DE_genes_rrbs_top.Rda")
  
  # Prom
  # DE_genes_rrbs_prom_top_mat <- unlist(output_rrbs[(output_rrbs$adj.P.Val<0.10) & (abs(output_rrbs$mean_T-output_rrbs$mean_C)>0.30) & (output_rrbs$func=="promoter"),
  #                                             grepl("gene",colnames(output_rrbs)) & !(grepl("gene_symbol",colnames(output_rrbs)))])
  # DE_genes_rrbs_prom_top <- DE_genes_rrbs_prom_top_mat[!is.na(DE_genes_rrbs_prom_top_mat)]
  DE_genes_rrbs_prom_top <- output_rrbs$gene[(output_rrbs$adj.P.Val<0.10) & (abs(output_rrbs$mean_T-output_rrbs$mean_C)>0.30) & (output_rrbs$func=="promoter")]
  DE_genes_rrbs_prom_top <- sort(toupper(unique(DE_genes_rrbs_prom_top)))[sort(toupper(unique(DE_genes_rrbs_prom_top)))!=""]
  DE_genes_rrbs_prom_top <- data.frame(gene=DE_genes_rrbs_prom_top,stringsAsFactors = F)
  DE_genes_rrbs_prom_top$aveLogFC <- getAveLogFC(DE_genes_rrbs_prom_top$gene,output_rrbs,grep("gene$",colnames(output_rrbs)),grep("logFC",colnames(output_rrbs)),
                                                 output_rrbs$func=="promoter")
  DE_genes_rrbs_prom_top$gene_symbol <- getGeneSymbol(DE_genes_rrbs_prom_top$gene,output_rrbs,grep("gene$",colnames(output_rrbs)),grep("gene_symbol_all",colnames(output_rrbs)),
                                                      output_rrbs$func=="promoter")
  DE_genes_rrbs_prom_top <- getOrthologues(DE_genes_rrbs_prom_top)
  save(DE_genes_rrbs_prom_top,file="DE_genes_rrbs_prom_top.Rda")
  
  # Lyko: Control vs cancer
  # All
  DE_genes_lyko_top <- lyko$ens_gene_id[(lyko$adj.P.Val <0.05) & (abs(lyko$Mean_T-lyko$Mean_C)>0.30)]
  DE_genes_lyko_top <- sort(toupper(unique(DE_genes_lyko_top)))[sort(toupper(unique(DE_genes_lyko_top)))!=""]
  DE_genes_lyko_top <- data.frame(gene=DE_genes_lyko_top,stringsAsFactors = F)
  DE_genes_lyko_top$aveLogFC <- getAveLogFC(DE_genes_lyko_top$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("logFC",colnames(lyko)))
  DE_genes_lyko_top$gene_symbol <- getGeneSymbol(DE_genes_lyko_top$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("gene$",colnames(lyko)))
  save(DE_genes_lyko_top,file="DE_genes_lyko_top.Rda")
  
  # Prom
  DE_genes_lyko_prom_top <- lyko$ens_gene_id[(lyko$adj.P.Val <0.05) & (abs(lyko$Mean_T-lyko$Mean_C)>0.30) & grepl("TSS",lyko$feature)]
  DE_genes_lyko_prom_top <- sort(toupper(unique(DE_genes_lyko_prom_top)))[sort(toupper(unique(DE_genes_lyko_prom_top)))!=""]
  DE_genes_lyko_prom_top <- data.frame(gene=DE_genes_lyko_prom_top,stringsAsFactors = F)
  DE_genes_lyko_prom_top$aveLogFC <- getAveLogFC(DE_genes_lyko_prom_top$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("logFC",colnames(lyko)),
                                                 grepl("TSS",lyko$feature))
  DE_genes_lyko_prom_top$gene_symbol <- getGeneSymbol(DE_genes_lyko_prom_top$gene,lyko,grep("ens_gene_id",colnames(lyko)),grep("gene$",colnames(lyko)),
                                                 grepl("TSS",lyko$feature))
  save(DE_genes_lyko_prom_top,file="DE_genes_lyko_prom_top.Rda")
  
  # Lyko: Keratinocyte-like vs cancer
  # all
  DE_genes_lyko_Ker_top <- lyko_Ker$ens_gene_id[(lyko_Ker$adj.P.Val <0.05) & (abs(lyko_Ker$Mean_T-lyko_Ker$Mean_C)>0.30)]
  DE_genes_lyko_Ker_top <- sort(toupper(unique(DE_genes_lyko_Ker_top)))[sort(toupper(unique(DE_genes_lyko_Ker_top)))!=""]
  DE_genes_lyko_Ker_top <- data.frame(gene=DE_genes_lyko_Ker_top,stringsAsFactors = F)
  DE_genes_lyko_Ker_top$aveLogFC <- getAveLogFC(DE_genes_lyko_Ker_top$gene,lyko_Ker,grep("ens_gene_id",colnames(lyko_Ker)),grep("logFC",colnames(lyko_Ker)))
  DE_genes_lyko_Ker_top$gene_symbol <- getGeneSymbol(DE_genes_lyko_Ker_top$gene,lyko_Ker,grep("ens_gene_id",colnames(lyko_Ker)),grep("gene$",colnames(lyko_Ker)))
  save(DE_genes_lyko_Ker_top,file="DE_genes_lyko_Ker_top.Rda")
  
  # Prom
  DE_genes_lyko_Ker_prom_top <- lyko_Ker$ens_gene_id[(lyko_Ker$adj.P.Val <0.05) & (abs(lyko_Ker$Mean_T-lyko_Ker$Mean_C)>0.30) & grepl("TSS",lyko_Ker$feature)]
  DE_genes_lyko_Ker_prom_top <- sort(toupper(unique(DE_genes_lyko_Ker_prom_top)))[sort(toupper(unique(DE_genes_lyko_Ker_prom_top)))!=""]
  DE_genes_lyko_Ker_prom_top <- data.frame(gene=DE_genes_lyko_Ker_prom_top,stringsAsFactors = F)
  DE_genes_lyko_Ker_prom_top$aveLogFC <- getAveLogFC(DE_genes_lyko_Ker_prom_top$gene,lyko_Ker,grep("ens_gene_id",colnames(lyko_Ker)),grep("logFC",colnames(lyko_Ker)),
                                                     grepl("TSS",lyko_Ker$feature))
  DE_genes_lyko_Ker_prom_top$gene_symbol <- getGeneSymbol(DE_genes_lyko_Ker_prom_top$gene,lyko_Ker,grep("ens_gene_id",colnames(lyko_Ker)),grep("gene$",colnames(lyko_Ker)),
                                                          grepl("TSS",lyko_Ker$feature))
  save(DE_genes_lyko_Ker_prom_top,file="DE_genes_lyko_Ker_prom_top.Rda")
  
  # Lyko: Stemcell-like vs cancer
  # all
  DE_genes_lyko_Stem_top <- lyko_Stem$ens_gene_id[(lyko_Stem$adj.P.Val <0.05) & (abs(lyko_Stem$Mean_T-lyko_Stem$Mean_C)>0.30)]
  DE_genes_lyko_Stem_top <- sort(toupper(unique(DE_genes_lyko_Stem_top)))[sort(toupper(unique(DE_genes_lyko_Stem_top)))!=""]
  DE_genes_lyko_Stem_top <- data.frame(gene=DE_genes_lyko_Stem_top,stringsAsFactors = F)
  DE_genes_lyko_Stem_top$aveLogFC <- getAveLogFC(DE_genes_lyko_Stem_top$gene,lyko_Stem,grep("ens_gene_id",colnames(lyko_Stem)),grep("logFC",colnames(lyko_Stem)))
  DE_genes_lyko_Stem_top$gene_symbol <- getGeneSymbol(DE_genes_lyko_Stem_top$gene,lyko_Stem,grep("ens_gene_id",colnames(lyko_Stem)),grep("gene$",colnames(lyko_Stem)))
  save(DE_genes_lyko_Stem_top,file="DE_genes_lyko_Stem_top.Rda")
  
  # Prom
  DE_genes_lyko_Stem_prom_top <- lyko_Stem$ens_gene_id[(lyko_Stem$adj.P.Val <0.05) & (abs(lyko_Stem$Mean_T-lyko_Stem$Mean_C)>0.30) & grepl("TSS",lyko_Stem$feature)]
  DE_genes_lyko_Stem_prom_top <- sort(toupper(unique(DE_genes_lyko_Stem_prom_top)))[sort(toupper(unique(DE_genes_lyko_Stem_prom_top)))!=""]
  DE_genes_lyko_Stem_prom_top <- data.frame(gene=DE_genes_lyko_Stem_prom_top,stringsAsFactors = F)
  DE_genes_lyko_Stem_prom_top$aveLogFC <- getAveLogFC(DE_genes_lyko_Stem_prom_top$gene,lyko_Stem,grep("ens_gene_id",colnames(lyko_Stem)),grep("logFC",colnames(lyko_Stem)),
                                                      grepl("TSS",lyko_Stem$feature))
  DE_genes_lyko_Stem_prom_top$gene_symbol <- getGeneSymbol(DE_genes_lyko_Stem_prom_top$gene,lyko_Stem,grep("ens_gene_id",colnames(lyko_Stem)),grep("gene$",colnames(lyko_Stem)),
                                                          grepl("TSS",lyko_Stem$feature))
  save(DE_genes_lyko_Stem_prom_top,file="DE_genes_lyko_Stem_prom_top.Rda")
} else {
  load("DE_genes_rrbs_top.Rda")
  load("DE_genes_rrbs_prom_top.Rda")
  
  load("DE_genes_lyko_top.Rda")
  load("DE_genes_lyko_prom_top.Rda")
  
  load("DE_genes_lyko_Ker_top.Rda")
  load("DE_genes_lyko_Ker_prom_top.Rda")
  
  load("DE_genes_lyko_Stem_top.Rda")
  load("DE_genes_lyko_Stem_prom_top.Rda")
}


## Tables_overlap
table_overlaps_GOI <- matrix(rep(NA,6*3),ncol=3)
colnames(table_overlaps_GOI) <- c("Cancer","keratinocyte-like","Stemcell-like")
rownames(table_overlaps_GOI) <- c("All N with ortholog","All N present in both","All N same direction",
                                  "Prom N with ortholog","Prom N present in both","Prom N same direction")

table_overlaps_GOI[1,] <- sum(DE_genes_rrbs_top$human_ortholog!="")
table_overlaps_GOI[4,] <- sum(DE_genes_rrbs_prom_top$human_ortholog!="")

table_overlaps_GOI_GS <- matrix(rep(NA,6*3),ncol=3)
colnames(table_overlaps_GOI_GS) <- c("Cancer","keratinocyte-like","Stemcell-like")
rownames(table_overlaps_GOI_GS) <- c("All N with ortholog","All N present in both","All N same direction",
                                     "Prom N with ortholog","Prom N present in both","Prom N same direction")

table_overlaps_GOI_GS[1,] <- length(unique(unlist(strsplit(DE_genes_rrbs_top$gene_symbol,","))))
table_overlaps_GOI_GS[4,] <- length(unique(unlist(strsplit(DE_genes_rrbs_prom_top$gene_symbol,","))))

## Control vs cancer
# all
GOI <- overlapSets(DE_genes_lyko_top,DE_genes_rrbs_top,lyko,output_rrbs)
table_overlaps_GOI[2,1] <- nrow(GOI)
table_overlaps_GOI[3,1] <- sum(sign(GOI$lyko_aveLogFC)!=sign(GOI$roth_aveLogFC))
GOI_GS <- overlapSetsGS(DE_genes_lyko_top,DE_genes_rrbs_top,lyko,output_rrbs)
table_overlaps_GOI_GS[2,1] <- nrow(GOI_GS)
table_overlaps_GOI_GS[3,1] <- sum(sign(GOI_GS$lyko_aveLogFC)!=sign(GOI_GS$roth_aveLogFC))

# prom
GOI_prom <- overlapSets(DE_genes_lyko_prom_top,DE_genes_rrbs_prom_top,lyko,output_rrbs)
table_overlaps_GOI[5,1] <- nrow(GOI_prom)
table_overlaps_GOI[6,1] <- sum(sign(GOI_prom$lyko_aveLogFC)!=sign(GOI_prom$roth_aveLogFC))
GOI_prom_GS <- overlapSetsGS(DE_genes_lyko_prom_top,DE_genes_rrbs_prom_top,lyko,output_rrbs)
table_overlaps_GOI_GS[5,1] <- nrow(GOI_prom_GS)
table_overlaps_GOI_GS[6,1] <- sum(sign(GOI_prom_GS$lyko_aveLogFC)!=sign(GOI_prom_GS$roth_aveLogFC))


## Ker vs cancer
# all
GOI_Ker <- overlapSets(DE_genes_lyko_Ker_top,DE_genes_rrbs_top,lyko_Ker,output_rrbs)
table_overlaps_GOI[2,2] <- nrow(GOI_Ker)
table_overlaps_GOI[3,2] <- sum(sign(GOI_Ker$lyko_aveLogFC)!=sign(GOI_Ker$roth_aveLogFC))
GOI_Ker_GS <- overlapSetsGS(DE_genes_lyko_Ker_top,DE_genes_rrbs_top,lyko_Ker,output_rrbs)
table_overlaps_GOI_GS[2,2] <- nrow(GOI_Ker_GS)
table_overlaps_GOI_GS[3,2] <- sum(sign(GOI_Ker_GS$lyko_aveLogFC)!=sign(GOI_Ker_GS$roth_aveLogFC))

# prom
GOI_Ker_prom <- overlapSets(DE_genes_lyko_Ker_prom_top,DE_genes_rrbs_prom_top,lyko_Ker,output_rrbs)
table_overlaps_GOI[5,2] <- nrow(GOI_Ker_prom)
table_overlaps_GOI[6,2] <- sum(sign(GOI_Ker_prom$lyko_aveLogFC)!=sign(GOI_Ker_prom$roth_aveLogFC))
GOI_Ker_prom_GS <- overlapSetsGS(DE_genes_lyko_Ker_prom_top,DE_genes_rrbs_prom_top,lyko_Ker,output_rrbs)
table_overlaps_GOI_GS[5,2] <- nrow(GOI_Ker_prom_GS)
table_overlaps_GOI_GS[6,2] <- sum(sign(GOI_Ker_prom_GS$lyko_aveLogFC)!=sign(GOI_Ker_prom_GS$roth_aveLogFC))


## Control vs Stem
# all
GOI_Stem <- overlapSets(DE_genes_lyko_Stem_top,DE_genes_rrbs_top,lyko_Stem,output_rrbs)
table_overlaps_GOI[2,3] <- nrow(GOI_Stem)
table_overlaps_GOI[3,3] <- sum(sign(GOI_Stem$lyko_aveLogFC)!=sign(GOI_Stem$roth_aveLogFC))
GOI_Stem_GS <- overlapSetsGS(DE_genes_lyko_Stem_top,DE_genes_rrbs_top,lyko_Stem,output_rrbs)
table_overlaps_GOI_GS[2,3] <- nrow(GOI_Stem_GS)
table_overlaps_GOI_GS[3,3] <- sum(sign(GOI_Stem_GS$lyko_aveLogFC)!=sign(GOI_Stem_GS$roth_aveLogFC))

# prom
GOI_Stem_prom <- overlapSets(DE_genes_lyko_Stem_prom_top,DE_genes_rrbs_prom_top,lyko_Stem,output_rrbs)
table_overlaps_GOI[5,3] <- nrow(GOI_Stem_prom)
table_overlaps_GOI[6,3] <- sum(sign(GOI_Stem_prom$lyko_aveLogFC)!=sign(GOI_Stem_prom$roth_aveLogFC))
GOI_Stem_prom_GS <- overlapSetsGS(DE_genes_lyko_Stem_prom_top,DE_genes_rrbs_prom_top,lyko_Stem,output_rrbs)
table_overlaps_GOI_GS[5,3] <- nrow(GOI_Stem_prom_GS)
table_overlaps_GOI_GS[6,3] <- sum(sign(GOI_Stem_prom_GS$lyko_aveLogFC)!=sign(GOI_Stem_prom_GS$roth_aveLogFC))


table_overlaps
table_overlaps_GS
table_overlaps_GOI
table_overlaps_GOI_GS

GOI
GOI_Stem[!(GOI_Stem$lyko_gene%in%GOI_Ker$lyko_gene),]
GOI_Ker[!(GOI_Ker$lyko_gene%in%GOI_Ker$lyko_gene),]

order_Ker <- c(4,1,2,3,6,5,7,8)
order_Stem_1 <- c(7,1,3,6,10,8,12,14)
order_Stem_2 <- c(2,9,5,13,15)

cbind(GOI_Ker$lyko_gene_symbol[order_Ker],round(GOI_Ker$roth_aveLogFC[order_Ker],2),
      round(-GOI_Ker$lyko_aveLogFC[order_Ker],2),round(-GOI_Stem$lyko_aveLogFC[order_Stem_1],2))

cbind(GOI_Stem$lyko_gene_symbol[order_Stem_2],round(GOI_Stem$roth_aveLogFC[order_Stem_2],2),
      round(-GOI_Stem$lyko_aveLogFC[order_Stem_2],2))

GOI_am_CpG_Stem <- cbind(lyko_Stem[(lyko_Stem$adj.P.Val <0.05) & (abs(lyko_Stem$Mean_T-lyko_Stem$Mean_C)>0.30) & (lyko_Stem$ens_gene_id%in%GOI_Stem$lyko_gene),c(1,2,5)],
                         lyko_Stem[(lyko_Stem$adj.P.Val <0.05) & (abs(lyko_Stem$Mean_T-lyko_Stem$Mean_C)>0.30) & (lyko_Stem$ens_gene_id%in%GOI_Stem$lyko_gene),13] - 
                           lyko_Stem[(lyko_Stem$adj.P.Val <0.05) & (abs(lyko_Stem$Mean_T-lyko_Stem$Mean_C)>0.30) & (lyko_Stem$ens_gene_id%in%GOI_Stem$lyko_gene),15])
colnames(GOI_am_CpG_Stem) <- c("chr","pos","gene","D % meth")
table(as.character(GOI_am_CpG_Stem$gene))
GOI_am_CpG_Ker <- cbind(lyko_Ker[(lyko_Ker$adj.P.Val <0.05) & (abs(lyko_Ker$Mean_T-lyko_Ker$Mean_C)>0.30) & (lyko_Ker$ens_gene_id%in%GOI_Ker$lyko_gene),c(1,2,5)],
                         lyko_Ker[(lyko_Ker$adj.P.Val <0.05) & (abs(lyko_Ker$Mean_T-lyko_Ker$Mean_C)>0.30) & (lyko_Ker$ens_gene_id%in%GOI_Ker$lyko_gene),13] - 
                           lyko_Ker[(lyko_Ker$adj.P.Val <0.05) & (abs(lyko_Ker$Mean_T-lyko_Ker$Mean_C)>0.30) & (lyko_Ker$ens_gene_id%in%GOI_Ker$lyko_gene),15])
colnames(GOI_am_CpG_Ker) <- c("chr","pos","gene","D % meth")
table(as.character(GOI_am_CpG_Ker$gene))

mask_am_CpG_rrbs <- rep(FALSE,nrow(output_rrbs))
for (i in 1:nrow(GOI_Stem)){
  mask_am_CpG_rrbs <- mask_am_CpG_rrbs + rowSums((output_rrbs[,grep("gene",colnames(output_rrbs))]==GOI_Stem$roth_gene[i]) & !is.na(output_rrbs[,grep("gene",colnames(output_rrbs))]))
}
mask_am_CpG_rrbs <- mask_am_CpG_rrbs>0

GOI_am_CpG_rrbs <- cbind(output_rrbs[(output_rrbs$adj.P.Val<0.10) & (abs(output_rrbs$mean_T-output_rrbs$mean_C)>0.30) & mask_am_CpG_rrbs,c("chr","locations","gene_symbol_all")],
                         output_rrbs[(output_rrbs$adj.P.Val<0.10) & (abs(output_rrbs$mean_T-output_rrbs$mean_C)>0.30) & mask_am_CpG_rrbs,"mean_T"] - 
                           output_rrbs[(output_rrbs$adj.P.Val<0.10) & (abs(output_rrbs$mean_T-output_rrbs$mean_C)>0.30) & mask_am_CpG_rrbs,"mean_C"])
colnames(GOI_am_CpG_rrbs) <- c("chr","pos","gene","D % meth")
table(as.character(GOI_am_CpG_rrbs$gene))

write.table("Mouse",file="CpG_in_GOI.txt",col.names = F,row.names = F,quote=F)
write.table(GOI_am_CpG_rrbs[sort(as.character(GOI_am_CpG_rrbs$gene),index.return=T)$ix,],file="CpG_in_GOI.txt",col.names = T,quote = F,row.names = F,sep="\t",append=T)
write.table("Human (Stem)",file="CpG_in_GOI.txt",col.names = F,row.names = F,quote=F,append=T)
write.table(GOI_am_CpG_Stem[sort(as.character(GOI_am_CpG_Stem$gene),index.return=T)$ix,],file="CpG_in_GOI.txt",col.names = T,quote = F,row.names = F,sep="\t",append=T)
write.table("Human (Ker)",file="CpG_in_GOI.txt",col.names = F,row.names = F,quote=F,append=T)
write.table(GOI_am_CpG_Ker[sort(as.character(GOI_am_CpG_Ker$gene),index.return=T)$ix,],file="CpG_in_GOI.txt",col.names = T,quote = F,row.names = F,sep="\t",append=T)

# FILIP1L checks (interactive): significant promoter (TSS) probes in human, CpGs in mouse
# lyko_Stem[(lyko_Stem$adj.P.Val <0.05) & (abs(lyko_Stem$Mean_T-lyko_Stem$Mean_C)>0.15) & grepl("TSS",lyko_Stem$feature) & (lyko_Stem$gene=="FILIP1L"),]
# lyko_Ker[(lyko_Ker$adj.P.Val <0.05) & (abs(lyko_Ker$Mean_T-lyko_Ker$Mean_C)>0.15) & grepl("TSS",lyko_Ker$feature) & (lyko_Ker$gene=="FILIP1L"),]
# 
# output_rrbs$gene[(output_rrbs$adj.P.Val<0.10) & (abs(output_rrbs$mean_T-output_rrbs$mean_C)>0.20) & (output_rrbs$func=="promoter") & (grepl("FILIP1L",toupper(output_rrbs$gene_symbol_all)))]
# output_rrbs$gene[(output_rrbs$func=="promoter") & (grepl("FILIP1L",toupper(output_rrbs$gene_symbol_all)))]
# output_rrbs[(grepl("FILIP1L",toupper(output_rrbs$gene_symbol_all))),1:14]


################
# makeBedGraph #
################

output_rrbs <- output_rrbs[!grepl("GL",output_rrbs$chr),]
output_rrbs <- output_rrbs[!grepl("JH",output_rrbs$chr),]
output_rrbs <- output_rrbs[!grepl("MT",output_rrbs$chr),]
output_rrbs <- output_rrbs[!grepl("PATCH",output_rrbs$chr),]
# bedGraph positions are 0-based and the end is not included: one CpG runs from locations-1 to locations
output_rrbs$start_bed <- output_rrbs$locations - 1
makeBedGraph(output_rrbs,c("chr","start_bed","locations","mean_T"),"RRBS_tumour","Meth_perc",col=c("blue","orange"),range=c(0,1),isSign = (output_rrbs$adj.P.Val<0.10))
makeBedGraph(output_rrbs,c("chr","start_bed","locations","mean_C"),"RRBS_control","Meth_perc",col=c("green","red"),range=c(0,1),isSign = (output_rrbs$adj.P.Val<0.10))




#################################################
## Clustering Stemcell-like vs keratinocyte-like ##
#################################################

# Keratin cluster
keratin <- rownames(annot_LIMMAout)[annot_LIMMAout[,1]==19 & (as.numeric(annot_LIMMAout[,2])>38800000 & as.numeric(annot_LIMMAout[,2])<39900000)]
Mval_keratin <- Mval[rownames(Mval)%in%keratin,]
betas_keratin <- betas(infdata_all)[rownames(infdata_all)%in%keratin,]
distr <- dist(t(Mval_keratin))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
plot(ddr)
heatmap.2(Mval_keratin)
heatmap.2(t(betas_keratin),tracecol=NA)

# TP63
tp63 <- rownames(annot_LIMMAout)[annot_LIMMAout[,1]==3 & (as.numeric(annot_LIMMAout[,2])>189631416 & as.numeric(annot_LIMMAout[,2])<189897279 )]
Mval_tp63 <- Mval[rownames(Mval)%in%tp63,]
betas_tp63 <- betas(infdata_all)[rownames(infdata_all)%in%tp63,]
distr <- dist(t(Mval_tp63))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
plot(ddr)
heatmap.2(Mval_tp63)
heatmap.2(t(betas_tp63),tracecol=NA)
sort()



# KRT5: chr12:52908359-52914243 chr12:52514575-52520687 gencode v19 52908359-52914471
# krt5 <- rownames(annot_LIMMAout)[annot_LIMMAout[,1]==12 & (as.numeric(annot_LIMMAout[,2])>52908359 & as.numeric(annot_LIMMAout[,2])<52914471)]
krt5 <- rownames(annot_LIMMAout)[annot_LIMMAout[,5]=="KRT5"]
Mval_krt5 <- Mval[rownames(Mval)%in%krt5,]
betas_krt5 <- betas(infdata_all)[rownames(infdata_all)%in%krt5,]
distr <- dist(t(Mval_krt5))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
# plot(ddr)
# heatmap.2(Mval_krt5)


# KRT14: 17:39738531-39743173:-1
# krt14 <- rownames(annot_LIMMAout)[annot_LIMMAout[,1]==17 & (as.numeric(annot_LIMMAout[,2])>39738531 & as.numeric(annot_LIMMAout[,2])<39743173)]
krt14 <- rownames(annot_LIMMAout)[annot_LIMMAout[,5]=="KRT14"]
Mval_krt14 <- Mval[rownames(Mval)%in%krt14,]
betas_krt14 <- betas(infdata_all)[rownames(infdata_all)%in%krt14,]
distr <- dist(t(Mval_krt14))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
# plot(ddr)
# heatmap.2(Mval_krt14)


# KRT15: 17:39669995-39678781:-1
# krt15 <- rownames(annot_LIMMAout)[annot_LIMMAout[,1]==17 & (as.numeric(annot_LIMMAout[,2])>39669995 & as.numeric(annot_LIMMAout[,2])<39678781)]
krt15 <- rownames(annot_LIMMAout)[annot_LIMMAout[,5]=="KRT15"]
Mval_krt15 <- Mval[rownames(Mval)%in%krt15,]
betas_krt15 <- betas(infdata_all)[rownames(infdata_all)%in%krt15,]
distr <- dist(t(Mval_krt15))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
# plot(ddr)
# heatmap.2(Mval_krt15)


# KRT80: 12:52562780-52585784:-1
# krt80 <- rownames(annot_LIMMAout)[annot_LIMMAout[,1]==12 & (as.numeric(annot_LIMMAout[,2])>52562780 & as.numeric(annot_LIMMAout[,2])<52585784)]
krt80 <- rownames(annot_LIMMAout)[annot_LIMMAout[,5]=="KRT80"]
Mval_krt80 <- Mval[rownames(Mval)%in%krt80,]
betas_krt80 <- betas(infdata_all)[rownames(infdata_all)%in%krt80,]
distr <- dist(t(Mval_krt80))
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
# plot(ddr)
# heatmap.2(Mval_krt15)

jpeg("KRT5.jpg")
heatmap.2(t(betas_krt5),tracecol=NA,col="bluered")
dev.off()
jpeg("KRT14.jpg")
heatmap.2(t(betas_krt14),tracecol=NA,col="bluered")
dev.off()
jpeg("KRT15.jpg")
heatmap.2(t(betas_krt15),tracecol=NA,col="bluered")
dev.off()
jpeg("KRT80.jpg")
heatmap.2(t(betas_krt80),tracecol=NA,col="bluered")
dev.off()




#KRT5 mouse: 15:101707070-101712898:-1
betas_krt5_rrbs <- betas_rrbs[annot_rrbs[,1]==15 & (as.numeric(annot_rrbs[,2])>101707070 & as.numeric(annot_rrbs[,2])<101712898),]
betas_krt5_rrbs <- betas_krt5_rrbs[!is.na(rowSums(betas_krt5_rrbs)),]
jpeg("KRT5_rrbs.jpg")
heatmap.2(t(betas_krt5_rrbs),tracecol=NA,col="bluered")
dev.off()

#KRT14 mouse: 11:100203162-100207548:-1
betas_krt14_rrbs <- betas_rrbs[annot_rrbs[,1]==11 & (as.numeric(annot_rrbs[,2])>100203162 & as.numeric(annot_rrbs[,2])<100207548),]
betas_krt14_rrbs <- betas_krt14_rrbs[!is.na(rowSums(betas_krt14_rrbs)),]
jpeg("KRT14_rrbs.jpg")
heatmap.2(t(betas_krt14_rrbs),tracecol=NA,col="bluered")
dev.off()

#KRT15 mouse: 11:100131758-100135928:-1
betas_krt15_rrbs <- betas_rrbs[annot_rrbs[,1]==11 & (as.numeric(annot_rrbs[,2])>100131758 & as.numeric(annot_rrbs[,2])<100135928),]
betas_krt15_rrbs <- betas_krt15_rrbs[!is.na(rowSums(betas_krt15_rrbs)),]
jpeg("KRT15_rrbs.jpg")
heatmap.2(t(betas_krt15_rrbs),tracecol=NA,col="bluered")
dev.off()

#KRT80 mouse: 15: 101,347,444-101,370,162
betas_krt80_rrbs <- betas_rrbs[annot_rrbs[,1]==15 & (as.numeric(annot_rrbs[,2])>101347444 & as.numeric(annot_rrbs[,2])<101370162),]
betas_krt80_rrbs <- betas_krt80_rrbs[!is.na(rowSums(betas_krt80_rrbs)),]
jpeg("KRT80_rrbs.jpg")
heatmap.2(t(betas_krt80_rrbs),tracecol=NA,col="bluered")
dev.off()

#keratin mouse: 15:101343355-102046669
betas_keratin_rrbs <- betas_rrbs[annot_rrbs[,1]==15 & (as.numeric(annot_rrbs[,2])>101343355 & as.numeric(annot_rrbs[,2])<102046669),]
betas_keratin_rrbs <- betas_keratin_rrbs[!is.na(rowSums(betas_keratin_rrbs)),]
jpeg("keratin_rrbs.jpg")
heatmap.2(t(betas_keratin_rrbs),tracecol=NA,col="bluered")
dev.off()


################
## CpG panel ##
################

CpG_panel <- read.csv("CpG_panel.csv",header=T,stringsAsFactors = F)

lyko_Ker_meant_CpGpanel <- lyko_Ker[rownames(lyko_Ker)%in%CpG_panel$IlmnID,c("Mean_T","Sd_T")]
lyko_Stem_meant_CpGpanel <- lyko_Stem[rownames(lyko_Stem)%in%CpG_panel$IlmnID,c("Mean_T","Sd_T")]

CpG_panel_sorted <- CpG_panel[sort(CpG_panel$IlmnID,index.return=T)$ix,]
lyko_Ker_meant_CpGpanel_sorted <- lyko_Ker_meant_CpGpanel[sort(rownames(lyko_Ker_meant_CpGpanel),index.return=T)$ix,]
lyko_Stem_meant_CpGpanel_sorted <- lyko_Stem_meant_CpGpanel[sort(rownames(lyko_Stem_meant_CpGpanel),index.return=T)$ix,]

dat_cpg_panel <- cbind(CpG_panel_sorted[-which(!(CpG_panel_sorted$IlmnID%in%rownames(lyko_Ker_meant_CpGpanel_sorted))),c(1,4,5,2,3)],
                       lyko_Ker_meant_CpGpanel_sorted,lyko_Stem_meant_CpGpanel_sorted)
colnames(dat_cpg_panel) <- c("ID","CHR","MAPINFO","COEF","CHANGE",
                             "MEAN_KER","SD_KER","MEAN_STEM","SD_STEM")

dim(lyko_Stem_meant_CpGpanel)
sum(rownames(lyko_Stem_meant_CpGpanel_sorted)==rownames(lyko_Ker_meant_CpGpanel_sorted))
sum(rownames(lyko_Stem_meant_CpGpanel_sorted)==dat_cpg_panel$ID)

dat_cpg_panel$DIFF <- dat_cpg_panel$MEAN_STEM-dat_cpg_panel$MEAN_KER
dat_cpg_panel$MAX_SD <- apply(rbind(dat_cpg_panel$SD_STEM,dat_cpg_panel$SD_KER),1,max)

mask_realdiff <- abs(dat_cpg_panel$DIFF)>dat_cpg_panel$MAX_SD

sum((dat_cpg_panel$CHANGE[mask_realdiff]=="UP") & (dat_cpg_panel$DIFF[mask_realdiff]>0))



lyko_Ker_inPanel <- lyko_Ker[rownames(lyko_Ker)%in%CpG_panel$IlmnID] 
CpG_panel[grepl("FILIP",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("TSPAN",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("HOXD3",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("ABR",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("GPC6",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("GAL3ST3",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("SEPT9",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("LBRA",CpG_panel$UCSC_RefGene_Name),]

CpG_panel[grepl("TCFL5",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("RGMA",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("GALNT13",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("SLC2A10",CpG_panel$UCSC_RefGene_Name),]
CpG_panel[grepl("R3HDM4",CpG_panel$UCSC_RefGene_Name),]



infdata_Sandoval <- readEPIC("idats_Sandoval/")


pfilter(infdata_Sandoval)


## Remove lines with NA's
print(dim(infdata_Sandoval))
infdata_Sandoval <- infdata_Sandoval[-which(rowSums(is.na(betas(infdata_Sandoval)))>0),]
print(dim(infdata_Sandoval))

annot_Sandoval <- probe.features[rownames(probe.features)%in%rownames(betas(infdata_Sandoval)),]

sum(rownames(betas(infdata_Sandoval))==rownames(annot_Sandoval))
length(rownames(annot_Sandoval))


## Load the annotation data in:
samp_annot_Sandoval_part1 <- read.table("E-MTAB-8542.sdrf.txt",header=T,sep="\t",stringsAsFactors=F)
samp_annot_Sandoval_part1 <- samp_annot_Sandoval_part1[!grepl("Grn",samp_annot_Sandoval_part1$Array.Data.File),]
samp_annot_Sandoval_part1 <- samp_annot_Sandoval_part1[c(1:5,23,7:10,12:16,6,11,17,19,20,18,21,22),]
# samples should be in this order: c(10:14,9,15,17,18:24,1,2,3,5,6,4,7,8)

samp_annot_Sandoval_part2 <- read.table("Sample_info_Sandoval.txt",header=T,sep="\t",stringsAsFactors=F)

samp_annot_Sandoval <- data.frame(File=unique(gsub("\\_[Red]+.idat$","",samp_annot_Sandoval_part1$Array.Data.File)),
                                  Sample_analysis=samp_annot_Sandoval_part1$Extract.Name,
                                  Sample_paper=samp_annot_Sandoval_part2$Sample,
                                  Age_analysis=samp_annot_Sandoval_part1$Characteristics.age.,
                                  Age_paper=samp_annot_Sandoval_part2$Age.at.diagnosis,
                                  Risk=samp_annot_Sandoval_part2$Type.of.lesion,
                                  stringsAsFactors = F)



samp_annot_Sandoval


# samp_annot_Sandoval$Risk <- NA
# samp_annot_Sandoval$Risk[grepl("HIGH",samp_annot_Sandoval_part2$Type.of.lesion)] <- "HIGH"
# samp_annot_Sandoval$Risk[grepl("LOW",samp_annot_Sandoval_part2$Type.of.lesion)] <- "LOW"


samp_annot_Sandoval_sorted <- NULL
for (i in 1:nrow(samp_annot_Sandoval)){
  samp_annot_Sandoval_sorted <- rbind(samp_annot_Sandoval_sorted,samp_annot_Sandoval[samp_annot_Sandoval$File==colnames(infdata_Sandoval)[i],])
}

samp_annot_Sandoval_sorted


unlink(paste("Risk_analysis/*",sep=""))
if(!dir.exists("Risk_analysis")) dir.create("Risk_analysis")
goi <- c("FILIP1L","TSPAN9","HOXD3","ABR","GPC6","GAL3ST3","SEPT9","MAB21L1","TCFL5","RGMA","GALNT13","SLC2A10","R3HDM4")
for (i in 1:length(goi)){
  print(goi[i])
  
  # betas_filip1l <- betas(infdata_Sandoval[annot_Sandoval$ens_gene_id=="ENSG00000168386",])
  # annot_Sandoval_filip1l <- annot_Sandoval[annot_Sandoval$ens_gene_id=="ENSG00000168386",]
  
  betas_goi <- betas(infdata_Sandoval[grepl(goi[i],annot_Sandoval$gene),])
  annot_Sandoval_goi <- annot_Sandoval[grepl(goi[i],annot_Sandoval$gene),]
  
  dim(betas_goi)
  
  # rownames(lyko_Ker)[(lyko_Ker$adj.P.Val<0.05) & (abs(lyko_Ker$Mean_T-lyko_Ker$Mean_C)>0.30)]
  
  dat_goi <- data.frame(pos=annot_Sandoval_goi$MAPINFO,
                        high=rowMeans(betas_goi[,(samp_annot_Sandoval_sorted$Risk=="HIGH") & !is.na(samp_annot_Sandoval_sorted$Risk)]),
                        low=rowMeans(betas_goi[,(samp_annot_Sandoval_sorted$Risk=="LOW") & !is.na(samp_annot_Sandoval_sorted$Risk)]),
                        stringsAsFactors = F)
  
  mask_spec_CpGs_Ker <- rep(1,nrow(dat_goi))
  mask_spec_CpGs_Stem <- rep(1,nrow(dat_goi))
  if (sum(grepl(goi[i],CpG_panel$UCSC_RefGene_Name))>0){
    CpG_panel_goi <- CpG_panel$IlmnID[grepl(goi[i],CpG_panel$UCSC_RefGene_Name)]
    print(CpG_panel_goi)
    mask_spec_CpGs_Ker[rownames(dat_goi)%in%CpG_panel_goi] <- 2
    mask_spec_CpGs_Stem[rownames(dat_goi)%in%CpG_panel_goi] <- 2
  }
  

  rn_lyko_oi_Ker_pos <- rownames(lyko_Ker)[grepl(goi[i],lyko_Ker$gene) & (lyko_Ker$adj.P.Val<0.05) & (lyko_Ker$Mean_T-lyko_Ker$Mean_C>0.30)]
  rn_lyko_oi_Ker_neg <- rownames(lyko_Ker)[grepl(goi[i],lyko_Ker$gene) & (lyko_Ker$adj.P.Val<0.05) & (lyko_Ker$Mean_T-lyko_Ker$Mean_C<(-0.30))]
  rn_lyko_oi_Stem_pos <- rownames(lyko_Stem)[grepl(goi[i],lyko_Stem$gene) & (lyko_Stem$adj.P.Val<0.05) & (lyko_Stem$Mean_T-lyko_Stem$Mean_C>0.30)]
  rn_lyko_oi_Stem_neg <- rownames(lyko_Stem)[grepl(goi[i],lyko_Stem$gene) & (lyko_Stem$adj.P.Val<0.05) & (lyko_Stem$Mean_T-lyko_Stem$Mean_C<(-0.30))]
  
  if (sum(rownames(lyko_Ker)%in%rn_lyko_oi_Ker_pos)>0){
    print(cbind(rn_lyko_oi_Ker_pos,lyko_Ker$Mean_T[rownames(lyko_Ker)%in%rn_lyko_oi_Ker_pos]-lyko_Ker$Mean_C[rownames(lyko_Ker)%in%rn_lyko_oi_Ker_pos]))
  }
  if (sum(rownames(lyko_Ker)%in%rn_lyko_oi_Ker_neg)>0){
    print(cbind(rn_lyko_oi_Ker_neg,lyko_Ker$Mean_T[rownames(lyko_Ker)%in%rn_lyko_oi_Ker_neg]-lyko_Ker$Mean_C[rownames(lyko_Ker)%in%rn_lyko_oi_Ker_neg]))
  }
  if (sum(rownames(lyko_Stem)%in%rn_lyko_oi_Stem_pos)>0){
    print(cbind(rn_lyko_oi_Stem_pos,lyko_Stem$Mean_T[rownames(lyko_Stem)%in%rn_lyko_oi_Stem_pos]-lyko_Stem$Mean_C[rownames(lyko_Stem)%in%rn_lyko_oi_Stem_pos]))
  }
  if (sum(rownames(lyko_Stem)%in%rn_lyko_oi_Stem_neg)>0){
    print(cbind(rn_lyko_oi_Stem_neg,lyko_Stem$Mean_T[rownames(lyko_Stem)%in%rn_lyko_oi_Stem_neg]-lyko_Stem$Mean_C[rownames(lyko_Stem)%in%rn_lyko_oi_Stem_neg]))
  }
  
  mask_spec_CpGs_Ker[rownames(dat_goi)%in%rn_lyko_oi_Ker_pos] <- 3
  mask_spec_CpGs_Ker[rownames(dat_goi)%in%rn_lyko_oi_Ker_neg] <- 4

  mask_spec_CpGs_Stem[rownames(dat_goi)%in%rn_lyko_oi_Stem_pos] <- 3
  mask_spec_CpGs_Stem[rownames(dat_goi)%in%rn_lyko_oi_Stem_neg] <- 4
  
  
  jpeg(paste("Risk_analysis/Risk_analysis_",goi[i],".jpg",sep=""),width=1024+256,height=512)
  par(mfrow=c(1,2))
  par(mar=c(4,4,4,10))
  plot(dat_goi$pos,dat_goi$high-dat_goi$low,
       col=c("black","red","blue","green")[mask_spec_CpGs_Ker],
       pch=c(1,16,16,16)[mask_spec_CpGs_Ker],
       main="Keratinocyte like",xlab="Position",ylab="Diff meth (high - low risk)")
  abline(h=0,lty=2)
  par(mar=c(4,4,4,10))
  plot(dat_goi$pos,dat_goi$high-dat_goi$low,
       col=c("black","red","blue","green")[mask_spec_CpGs_Stem],
       pch=c(1,16,16,16)[mask_spec_CpGs_Stem],
       main="Stemcell like",xlab="Position",ylab="Diff meth (high - low risk)")
  abline(h=0,lty=2)
  par(xpd=T)
  legend(par("usr")[2]+500,par("usr")[4],c("sign high vs low","pos diff lyko (>30%)","neg diff lyko (>30%)"),
         pch=16,col=c("red","blue","green"))
  dev.off()

}


## DM analysis
##
#########
#########

samp_annot_Sandoval_sorted$Risk_simple <- gsub("\\-[A-Z\\ ]+$","",samp_annot_Sandoval_sorted$Risk)

samp_annot_Sandoval_sorted_risk <- samp_annot_Sandoval_sorted[!grepl("ACT",samp_annot_Sandoval_sorted$Risk_simple),]
infdata_Sandoval_risk <- infdata_Sandoval[,!grepl("ACT",samp_annot_Sandoval_sorted$Risk_simple)]

#   design 1: control vs tumor
#   design 2: Stemcell like vs keratinocyte like
risk <- as.factor(samp_annot_Sandoval_sorted_risk$Risk_simple)
design3 <- model.matrix(~0+risk)

Mval_Sandoval <- log2((betas(infdata_Sandoval_risk)+ct)/(1-betas(infdata_Sandoval_risk)+ct))

## normal vs cancer
cont.matrix <- makeContrasts((riskHIGH)-(riskLOW),levels=design3)
fit <- lmFit(Mval_Sandoval,design3)
fit <- contrasts.fit(fit,cont.matrix)
fit <- eBayes(fit)




jpeg("Plots/MA_highvslowrisk.jpg")
limma::plotMA(fit)
dev.off()

LIMMAout_highvslow <- topTable(fit,adjust="BH",number=nrow(infdata_Sandoval))
LIMMAout_highvslow <- LIMMAout_highvslow[!is.na(LIMMAout_highvslow$adj.P.Val),]

jpeg("Plots/MA_highvslowrisk_DE.jpg")
plot(LIMMAout_highvslow$AveExpr,LIMMAout_highvslow$logFC,pch=20,xlab="Mean Log Expression",ylab="Log Fold Change")
points(LIMMAout_highvslow$AveExpr[(LIMMAout_highvslow$adj.P.Val<0.05)],
       LIMMAout_highvslow$logFC[(LIMMAout_highvslow$adj.P.Val<0.05)],
       col="red",pch=20)
abline(h=0,col="grey")
legend("topright",c("non-DE genes","DE genes"),pch=20,col=c("black","red"))
dev.off()

jpeg("Plots/Pvalue_dist.jpg")
hist(LIMMAout_highvslow$P.Value)
dev.off()

infdata_Sandoval_risk <- infdata_Sandoval_risk[rownames(annot_Sandoval)%in%rownames(LIMMAout_highvslow),]
annot_Sandoval <- annot_Sandoval[rownames(annot_Sandoval)%in%rownames(LIMMAout_highvslow),]

LIMMAout_highvslow <- LIMMAout_highvslow[sort(rownames(LIMMAout_highvslow),index.return=T)$ix,]
index_highvslow  <- sort(LIMMAout_highvslow$P.Value,index.return=T)$ix
LIMMAout_highvslow <- LIMMAout_highvslow[sort(LIMMAout_highvslow$P.Value,index.return=T)$ix,]
annot_LIMMAout_highvslow <- annot_Sandoval[index_highvslow ,]

paste(length(unique(annot_LIMMAout_highvslow$gene[LIMMAout_highvslow$adj.P.Val<0.05]))-1,"out of",length(unique(annot_LIMMAout_highvslow$gene))-1)


sandoval <- cbind(annot_LIMMAout_highvslow,
                  rowMeans(betas(infdata_Sandoval_risk)[,grepl("HIGH",risk)])[index_highvslow],
                  apply(betas(infdata_Sandoval_risk)[,grepl("HIGH",risk)],1,sd)[index_highvslow],
                  rowMeans(betas(infdata_Sandoval_risk)[,!grepl("HIGH",risk)])[index_highvslow],
                  apply(betas(infdata_Sandoval_risk)[,!grepl("HIGH",risk)],1,sd)[index_highvslow],
                  LIMMAout_highvslow)
colnames(sandoval) <- c(colnames(annot_LIMMAout_highvslow),"Mean_HIGH","Sd_HIGH","Mean_LOW","Sd_LOW",colnames(LIMMAout_highvslow))

sum((sandoval$adj.P.Val<0.05) & (grepl("FILIP1L",sandoval$gene)))


sign_CpG_filip1l <- rownames(sandoval)[(sandoval$adj.P.Val<0.05) & (grepl("FILIP1L",sandoval$gene))]

print("FILIP1L")

betas_filip1l <- betas(infdata_Sandoval[rownames(annot_Sandoval)%in%sign_CpG_filip1l,])
annot_Sandoval_filip1l <- annot_Sandoval[rownames(annot_Sandoval)%in%sign_CpG_filip1l,]

dim(betas_filip1l)

dat_goi <- data.frame(pos=annot_Sandoval_filip1l$MAPINFO,
                      high=rowMeans(betas_filip1l[,(samp_annot_Sandoval_sorted$Risk_simple=="HIGH") & !is.na(samp_annot_Sandoval_sorted$Risk)]),
                      low=rowMeans(betas_filip1l[,(samp_annot_Sandoval_sorted$Risk_simple=="LOW") & !is.na(samp_annot_Sandoval_sorted$Risk)]),
                      stringsAsFactors = F)

mask_spec_CpGs_Ker <- rep(1,nrow(dat_goi))
mask_spec_CpGs_Stem <- rep(1,nrow(dat_goi))


rn_lyko_oi_Ker_pos <- rownames(lyko_Ker)[grepl("FILIP1L",lyko_Ker$gene) & (lyko_Ker$adj.P.Val<0.05) & (lyko_Ker$Mean_T-lyko_Ker$Mean_C>0.30)]
rn_lyko_oi_Ker_neg <- rownames(lyko_Ker)[grepl("FILIP1L",lyko_Ker$gene) & (lyko_Ker$adj.P.Val<0.05) & (lyko_Ker$Mean_T-lyko_Ker$Mean_C<(-0.30))]
rn_lyko_oi_Stem_pos <- rownames(lyko_Stem)[grepl("FILIP1L",lyko_Stem$gene) & (lyko_Stem$adj.P.Val<0.05) & (lyko_Stem$Mean_T-lyko_Stem$Mean_C>0.30)]
rn_lyko_oi_Stem_neg <- rownames(lyko_Stem)[grepl("FILIP1L",lyko_Stem$gene) & (lyko_Stem$adj.P.Val<0.05) & (lyko_Stem$Mean_T-lyko_Stem$Mean_C<(-0.30))]

if (sum(rownames(lyko_Ker)%in%rn_lyko_oi_Ker_pos)>0){
  print(cbind(rn_lyko_oi_Ker_pos,lyko_Ker$Mean_T[rownames(lyko_Ker)%in%rn_lyko_oi_Ker_pos]-lyko_Ker$Mean_C[rownames(lyko_Ker)%in%rn_lyko_oi_Ker_pos]))
}
if (sum(rownames(lyko_Ker)%in%rn_lyko_oi_Ker_neg)>0){
  print(cbind(rn_lyko_oi_Ker_neg,lyko_Ker$Mean_T[rownames(lyko_Ker)%in%rn_lyko_oi_Ker_neg]-lyko_Ker$Mean_C[rownames(lyko_Ker)%in%rn_lyko_oi_Ker_neg]))
}
if (sum(rownames(lyko_Stem)%in%rn_lyko_oi_Stem_pos)>0){
  print(cbind(rn_lyko_oi_Stem_pos,lyko_Stem$Mean_T[rownames(lyko_Stem)%in%rn_lyko_oi_Stem_pos]-lyko_Stem$Mean_C[rownames(lyko_Stem)%in%rn_lyko_oi_Stem_pos]))
}
if (sum(rownames(lyko_Stem)%in%rn_lyko_oi_Stem_neg)>0){
  print(cbind(rn_lyko_oi_Stem_neg,lyko_Stem$Mean_T[rownames(lyko_Stem)%in%rn_lyko_oi_Stem_neg]-lyko_Stem$Mean_C[rownames(lyko_Stem)%in%rn_lyko_oi_Stem_neg]))
}

mask_spec_CpGs_Ker[rownames(dat_goi)%in%rn_lyko_oi_Ker_pos] <- 3
mask_spec_CpGs_Ker[rownames(dat_goi)%in%rn_lyko_oi_Ker_neg] <- 4

mask_spec_CpGs_Stem[rownames(dat_goi)%in%rn_lyko_oi_Stem_pos] <- 3
mask_spec_CpGs_Stem[rownames(dat_goi)%in%rn_lyko_oi_Stem_neg] <- 4


jpeg(paste("Risk_analysis_FILIP1L_sign.jpg",sep=""),width=1024,height=512)
# par(mfrow=c(1,2))
# par(mar=c(4,4,4,10))
# plot(dat_goi$pos,dat_goi$high-dat_goi$low,
#      col=c("black","red","blue","green")[mask_spec_CpGs_Ker],
#      pch=c(1,16,16,16)[mask_spec_CpGs_Ker],
#      main="Keratinocyte like",xlab="Position",ylab="Diff meth (high - low risk)")
# abline(h=0,lty=2)
par(mar=c(4,4,4,12))
plot(dat_goi$pos,dat_goi$high-dat_goi$low,
     col=c("black","red","blue","green")[mask_spec_CpGs_Stem],
     pch=c(1,16,16,16)[mask_spec_CpGs_Stem],
     xlab="Position",ylab="Diff meth (high - low risk)")
abline(h=0,lty=2)
par(xpd=T)
legend(par("usr")[2]+2500,par("usr")[4],c("hypermethylated (>30%)","hypomethylated (>30%)"),
       pch=16,col=c("blue","green"))
dev.off()


lyko_Ker[rownames(lyko_Ker)%in%c(rownames(dat_goi)[mask_spec_CpGs_Stem>1,])]
lyko[rownames(lyko)%in%c(rownames(dat_goi)[mask_spec_CpGs_Stem>1,])]
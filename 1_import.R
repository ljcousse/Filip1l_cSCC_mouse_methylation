library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
library(yaml)

# Directories are read from environment.yml (copy environment_example.yml and fill in the local paths)
if (!file.exists("environment.yml")){
  stop("environment.yml not found, copy environment_example.yml and fill in the local paths")
}
env <- read_yaml("environment.yml")

dir_pool20 <- env$dir_pool20
dir_pool16 <- env$dir_pool16
dir_work <- env$dir_work


##############
#Import data #
##############

###B files###
##grep all files individually for both folders

#first folder
setwd(dir_pool20)
getwd()
only_B_20<-grep("B", dir())
only_B_files_20<-dir()[only_B_20]
only_B_files_20
only_bis_20<-grep("bismark", only_B_files_20)
only_bis_files_20<-only_B_files_20[only_bis_20]
pool_20<-only_bis_files_20
pool_20<-pool_20[grep(".cov", pool_20)]
b_16<-pool_20[grep("16_", pool_20)]
b_168<-pool_20[grep("168", pool_20)]
b_35<-pool_20[grep("35_", pool_20)]
b_352<-pool_20[grep("352", pool_20)]
b_44<-pool_20[grep("B44_", pool_20)]
b_444<-pool_20[grep("444", pool_20)]
b_54<-pool_20[grep("54_", pool_20)]
b_545<-pool_20[grep("545", pool_20)]
b_56_1<-pool_20[grep("56_", pool_20)]
b_562_1<-pool_20[grep("562", pool_20)]

#2nd folder
setwd(dir_pool16)
getwd()
only_B_16<-grep("B", dir())
only_B_files_16<-dir()[only_B_16]
only_B_files_16
only_bis_16<-grep("bismark", only_B_files_16)
only_bis_files_16<-only_B_files_16[only_bis_16]
pool_16<-only_bis_files_16
pool_16<-pool_16[grep(".cov", pool_16)]
b_45<-pool_16[grep("45_", pool_16)]
b_77<-pool_16[grep("77_", pool_16)]
b_85<-pool_16[grep("85_", pool_16)]
b_851<-pool_16[grep("851", pool_16)]
b_87<-pool_16[grep("87_", pool_16)]
b_874<-pool_16[grep("874", pool_16)]
b_56_2<-pool_16[grep("56_", pool_16)]
b_562_2<-pool_16[grep("562", pool_16)]

#combine in the correct order

### different folders, reading files in folder by folder, then combining
tumor_RRBS_samples_20<-c(b_168, b_352, b_444, b_545, b_562_1)
tumor_RRBS_samples_20new<-paste(dir_pool20, tumor_RRBS_samples_20, sep="/")
tumor_RRBS_samples_16<-c(b_562_2, b_851, b_874)
tumor_RRBS_samples_16new<-paste(dir_pool16, tumor_RRBS_samples_16, sep="/")
tumor_RRBS_samples<-c(tumor_RRBS_samples_20new, tumor_RRBS_samples_16new)
tumor_RRBS_samples

matcon_RRBS_samples_20<-c(b_16, b_35, b_44, b_54, b_56_1)
matcon_RRBS_samples_20new<-paste(dir_pool20, matcon_RRBS_samples_20, sep="/")
matcon_RRBS_samples_16<-c(b_56_2, b_85, b_87)
matcon_RRBS_samples_16new<-paste(dir_pool16, matcon_RRBS_samples_16, sep="/")
matcon_RRBS_samples<-c(matcon_RRBS_samples_20new, matcon_RRBS_samples_16new)
matcon_RRBS_samples

ihs_RRBS_samples<-paste(dir_pool16, c(b_45, b_77), sep="/")
ihs_RRBS_samples

all_RRBS_samples<-c(tumor_RRBS_samples, matcon_RRBS_samples, ihs_RRBS_samples)
all_RRBS_samples

#no writing permission... creation colData
#Matched control
b16<-c("B16", "20", "mc", "n")
b35<-c("B35", "20", "mc", "n")
b44<-c("B44", "20", "mc", "n")
b54<-c("B54", "20", "mc", "n")
b56_1<-c("B56_1", "20", "mc", "n")
b56_2<-c("B56_2", "16", "mc", "n")
b85<-c("B85", "16", "mc", "n")
b87<-c("B87", "16", "mc", "n")
#Tumor
b168<-c("B168", "20", "T", "n")
b352<-c("B352", "20", "T", "n")
b444<-c("B444", "20", "T", "n")
b545<-c("B545", "20", "T", "n")
b562_1<-c("B562_1", "20", "T", "n")
b562_2<-c("B562_2", "16", "T", "n")
b851<-c("B851", "16", "T", "n")
b874<-c("B874", "16", "T", "n")
#indipendant healthy skin
b45<-c("B45", "16", "iHS", "n")
b77<-c("B77", "16", "iHS", "n")

###OX files###
##grep all files individually for both folders

#first folder
setwd(dir_pool20)
getwd()
only_ox_20<-grep("O", dir())
only_ox_files_20<-dir()[only_ox_20]
only_ox_files_20
only_bis_20<-grep("bismark", only_ox_files_20)
only_bis_files_20<-only_ox_files_20[only_bis_20]
pool_20<-only_bis_files_20
pool_20<-pool_20[grep(".cov", pool_20)]
o_16<-pool_20[grep("16_", pool_20)]
o_168<-pool_20[grep("168", pool_20)]
o_35<-pool_20[grep("35_", pool_20)]
o_352<-pool_20[grep("352", pool_20)]
o_44<-pool_20[grep("O44_", pool_20)]
o_444<-pool_20[grep("444", pool_20)]
o_54<-pool_20[grep("54_", pool_20)]
o_545<-pool_20[grep("545", pool_20)]
o_56_1<-pool_20[grep("56_", pool_20)]
o_562_1<-pool_20[grep("562", pool_20)]


#2nd folder
setwd(dir_pool16)
getwd()
only_ox_16<-grep("O", dir())
only_ox_files_16<-dir()[only_ox_16]
only_ox_files_16
only_bis_16<-grep("bismark", only_ox_files_16)
only_bis_files_16<-only_ox_files_16[only_bis_16]
pool_16<-only_bis_files_16
pool_16<-pool_16[grep(".cov", pool_16)]
o_45<-pool_16[grep("45_", pool_16)]
o_77<-pool_16[grep("77_", pool_16)]
o_85<-pool_16[grep("85_", pool_16)]
o_851<-pool_16[grep("851", pool_16)]
o_87<-pool_16[grep("87_", pool_16)]
o_874<-pool_16[grep("874", pool_16)]
o_56_2<-pool_16[grep("56_", pool_16)]
o_562_2<-pool_16[grep("562", pool_16)]

#combine in the correct order

### different folders, reading files in folder by folder, then combining
tumor_oxRRBS_samples_20<-c(o_168, o_352, o_444, o_545, o_562_1)
tumor_oxRRBS_samples_20new<-paste(dir_pool20, tumor_oxRRBS_samples_20, sep="/")
tumor_oxRRBS_samples_16<-c(o_562_2, o_851, o_874)
tumor_oxRRBS_samples_16new<-paste(dir_pool16, tumor_oxRRBS_samples_16, sep="/")
tumor_oxRRBS_samples<-c(tumor_oxRRBS_samples_20new, tumor_oxRRBS_samples_16new)
tumor_oxRRBS_samples

matcon_oxRRBS_samples_20<-c(o_16, o_35, o_44, o_54, o_56_1)
matcon_oxRRBS_samples_20new<-paste(dir_pool20, matcon_oxRRBS_samples_20, sep="/")
matcon_oxRRBS_samples_16<-c(o_56_2, o_85, o_87)
matcon_oxRRBS_samples_16new<-paste(dir_pool16, matcon_oxRRBS_samples_16, sep="/")
matcon_oxRRBS_samples<-c(matcon_oxRRBS_samples_20new, matcon_oxRRBS_samples_16new)
matcon_oxRRBS_samples

ihs_oxRRBS_samples<-paste(dir_pool16, c(o_45, o_77), sep="/")
ihs_oxRRBS_samples

all_oxRRBS_samples<-c(tumor_oxRRBS_samples, matcon_oxRRBS_samples, ihs_oxRRBS_samples)
all_oxRRBS_samples

#no writing permission... creation colData
#Matched control
o16<-c("O16", "20", "mc", "y")
o35<-c("O35", "20", "mc", "y")
o44<-c("O44", "20", "mc", "y")
o54<-c("O54", "20", "mc", "y")
o56_1<-c("O56_1", "20", "mc", "y")
o56_2<-c("O56_2", "16", "mc", "y")
o85<-c("O85", "16", "mc", "y")
o87<-c("O87", "16", "mc", "y")
#Tumor
o168<-c("O168", "20", "T", "y")
o352<-c("O352", "20", "T", "y")
o444<-c("O444", "20", "T", "y")
o545<-c("O545", "20", "T", "y")
o562_1<-c("O562_1", "20", "T", "y")
o562_2<-c("O562_2", "16", "T", "y")
o851<-c("O851", "16", "T", "y")
o874<-c("O874", "16", "T", "y")
#indipendant healthy skin
o45<-c("O45", "16", "iHS", "y")
o77<-c("O77", "16", "iHS", "y")

#DataFrame
annot<-c("ID", "pool", "type", "Ox?", b168, b352, b444, b545, b562_1, b562_2, b851, b874, b16, b35, b44, b54, b56_1, b56_2, b85, b87, b45, b77,
                                      o168, o352, o444, o545, o562_1, o562_2, o851, o874, o16, o35, o44, o54, o56_1, o56_2, o85, o87, o45, o77)
annot<-matrix(annot,nrow=37,ncol=4,byrow=T)
colnames(annot)<-annot[1,]
annot<-annot[-1,]
rownames(annot)<-annot[,1]
annot_df<-as.data.frame(annot)
annot_df<-annot_df[c(2,3,4)]
annot_df

#Bismark file import
all<-c(all_RRBS_samples, all_oxRRBS_samples)
allRRBS<-readBismark(all,colData=annot_df)

#save
setwd(dir_work)
save(allRRBS, file="allRRBS.rdata")
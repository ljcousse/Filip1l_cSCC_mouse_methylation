library(yaml)

# Directories are read from environment.yml (copy environment_example.yml and fill in the local paths)
if (!file.exists("environment.yml")){
  stop("environment.yml not found, copy environment_example.yml and fill in the local paths")
}
env <- read_yaml("environment.yml")

dir_work <- env$dir_work
setwd(dir_work)

links <- c("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3634nnn/GSM3634100/suppl/GSM3634100_WT1_H3K27Ac_peaks.txt.gz",
           "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3634nnn/GSM3634103/suppl/GSM3634103_WT2_H3K27Ac_peaks.txt.gz",
           "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3634nnn/GSM3634102/suppl/GSM3634102_WT1_H3K4me3_peaks.txt.gz",
           "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3634nnn/GSM3634105/suppl/GSM3634105_WT2_H3K4me3_peaks.txt.gz",
           "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3634nnn/GSM3634101/suppl/GSM3634101_WT1_H3K27me3_peaks.txt.gz",
           "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3634nnn/GSM3634104/suppl/GSM3634104_WT2_H3K27me3_peaks.txt.gz")

for (i in 1:length(links)){
  system(paste("wget ",links[i],sep=""))
}
system("wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz")
system("pigz -d *.gz")



file_names <- c("GSM3634100_WT1_H3K27Ac_peaks.txt","GSM3634103_WT2_H3K27Ac_peaks.txt",
                "GSM3634102_WT1_H3K4me3_peaks.txt","GSM3634105_WT2_H3K4me3_peaks.txt",
                "GSM3634101_WT1_H3K27me3_peaks.txt","GSM3634104_WT2_H3K27me3_peaks.txt")
chroms <- paste(rep("chr",length(c(1:19,"X","Y"))),c(1:19,"X","Y"),sep="")
for (i in 1:length(file_names)){
  peak_dat <- read.table(file_names[i],header=T,stringsAsFactors = F)
  sample_name <- gsub("\\_[A-Z,0-9]+\\_[A-Z,a-z,0-9]+\\_peaks.txt","",file_names[i])
  modification <- gsub("_peaks.txt$","",file_names[i])
  modification <- gsub("^[A-Z,0-9]+\\_[A-Z,0-9]+\\_","",modification)
  header_track <- paste("track type=narrowPeak name=",sample_name,"_",modification,sep="")
  print(file_names[i])
  print(sample_name)
  print(modification)
  print(header_track)
  
  narrowpeak_file <- cbind(paste(rep("chr",length(peak_dat$chr)),peak_dat$chr,sep=""),
                           peak_dat$start,peak_dat$end,peak_dat$name,
                           rep(0,length(peak_dat$chr)),
                           rep(".",length(peak_dat$chr)),
                           peak_dat$fold_enrichment,
                           peak_dat$'X.log10.pvalue.',
                           peak_dat$'X.log10.qvalue.',
                           rep(-1,length(peak_dat$chr)))
  colnames(narrowpeak_file) <- c("chrom","chromStart","chromEnd","name",
                                 "score","strand","signalValue","pValue",
                                 "qValue","peak")   
  narrowpeak_file <- data.frame(narrowpeak_file,stringsAsFactors = F)
  
  print(dim(narrowpeak_file))
  narrowpeak_file <- narrowpeak_file[narrowpeak_file$chrom%in%chroms,]
  print(dim(narrowpeak_file))
  print(table(narrowpeak_file$chrom))
  
  write.table(header_track,paste(sample_name,"_mm10_narrowPeak.txt",sep=""),col.names = F,row.names = F,quote = F)
  write.table(narrowpeak_file,paste(sample_name,"_mm10_narrowPeak.txt",sep=""),col.names = F,row.names = F,quote = F,sep="\t",append=T)
  
  chr <- narrowpeak_file$chrom
  pos1 <- narrowpeak_file$chromStart
  pos2 <- narrowpeak_file$chromEnd
  
  toliftover <- cbind(chr,pos1,pos2) 
  write.table(toliftover,file="PosToLiftover.txt",sep="\t",col.names = F,row.names = F,quote = F)
  
  system("liftOver PosToLiftover.txt mm10ToMm9.over.chain PosFromLiftover.txt PosNoLiftover.txt")
  
  new_pos <- read.table("PosFromLiftover.txt",sep="\t",header=F,stringsAsFactors = F)
  rem_pos <- read.table("PosNoLiftover.txt",sep="\t",header=F,stringsAsFactors = F)
  
  system("rm PosFromLiftover.txt")
  system("rm PosNoLiftover.txt")
  system("rm PosToLiftover.txt")
  
  # get pos that had a succesfull liftover
  pos_found <- !(paste(chr,pos1,sep="_")%in%paste(rem_pos$V1,rem_pos$V2,sep="_"))
  
  narrowpeak_file_mm9 <- narrowpeak_file[pos_found,]
  narrowpeak_file_mm9$chrom <- new_pos$V1
  narrowpeak_file_mm9$chromStart <- new_pos$V2
  narrowpeak_file_mm9$chromEnd <- new_pos$V3
  
  write.table(header_track,paste(sample_name,"_mm9_narrowPeak.txt",sep=""),col.names = F,row.names = F,quote = F)
  write.table(narrowpeak_file_mm9,paste(sample_name,"_mm9_narrowPeak.txt",sep=""),col.names = F,row.names = F,quote = F,sep="\t",append=T)
}

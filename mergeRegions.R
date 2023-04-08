sample_name <- commandArgs(TRUE)

bedgraph_sorted <- read.table(paste("chipseq_tracks/",sample_name,"_mm10_sorted.bedGraph",sep=""),header=F,stringsAsFactors = F,sep="\t")

options(warn=2)

to_merge <- which(bedgraph_sorted$V2[2:nrow(bedgraph_sorted)]<bedgraph_sorted$V3[1:(nrow(bedgraph_sorted)-1)] &
                    bedgraph_sorted$V1[2:nrow(bedgraph_sorted)]==bedgraph_sorted$V1[1:(nrow(bedgraph_sorted)-1)])
while(length(to_merge)!=0){
  print(length(to_merge))
  to_remove <- NULL
  rows2add <- NULL
  overall_counter <- 1
  perc <- 0.05
  for (i in to_merge){
    if ((overall_counter/length(to_merge))>perc){
      print(paste(perc*100,"% done",sep=""))
      perc <- perc + 0.05
    }
    overall_counter <- overall_counter + 1
    bedgraph_sorted_part <- bedgraph_sorted[(bedgraph_sorted$V1==bedgraph_sorted$V1[i]) &
                                              (bedgraph_sorted$V2>=bedgraph_sorted$V2[i]) & 
                                              (bedgraph_sorted$V2<bedgraph_sorted$V3[i]),]
    if (nrow(bedgraph_sorted_part)==1){
      next
    }
    if (nrow(bedgraph_sorted_part)==2){
      if (bedgraph_sorted_part$V4[1]==bedgraph_sorted_part$V4[2]){
        to_remove <- c(to_remove,i+1)
      } else {
        min_V2 <- bedgraph_sorted_part$V2[rownames(bedgraph_sorted_part)!=i]
        max_V3 <- bedgraph_sorted_part$V3[rownames(bedgraph_sorted_part)!=i]
        if (min_V2>bedgraph_sorted$V2[i] & max_V3<bedgraph_sorted$V2[i]){
          rows2add_part <- rbind(bedgraph_sorted[i,],bedgraph_sorted[i,])
          rows2add_part$V3[1] <- min_V2
          rows2add_part$V2[2] <- max_V3
          to_remove <- c(to_remove,i)
          rows2add <- rbind(rows2add,rows2add_part)
        } else if (min_V2>bedgraph_sorted$V2[i]) {
          bedgraph_sorted$V3[i] <- min_V2
          
        } else if (max_V3<bedgraph_sorted$V3[i]) {
          bedgraph_sorted$V2[i] <- max_V3
          
        } else {
          # print(i)
          # print(bedgraph_sorted_part)
          to_remove <- c(to_remove,i)
        }
      }
  
    } else {
      min_V2 <- min(bedgraph_sorted_part$V2[rownames(bedgraph_sorted_part)!=i])
      max_V3 <- max(bedgraph_sorted_part$V3[rownames(bedgraph_sorted_part)!=i])
      if (min_V2>bedgraph_sorted$V2[i] & max_V3<bedgraph_sorted$V2[i]){
        rows2add_part <- rbind(bedgraph_sorted[i,],bedgraph_sorted[i,])
        rows2add_part$V3[1] <- min_V2
        rows2add_part$V2[2] <- max_V3
        to_remove <- c(to_remove,i)
        rows2add <- rbind(rows2add,rows2add_part)
      } else if (min_V2>bedgraph_sorted$V2[i]) {
        bedgraph_sorted$V3[i] <- min_V2
        
      } else if (max_V3<bedgraph_sorted$V3[i]) {
        bedgraph_sorted$V2[i] <- max_V3
        
      } else {
        # print(i)
        # print(bedgraph_sorted_part)
        to_remove <- c(to_remove,i)
      }
  
    }
  }
  if(length(to_remove)!=0){
    bedgraph_sorted<- bedgraph_sorted[-to_remove,]
    bedgraph_sorted <- rbind(bedgraph_sorted,rows2add)
    # bedgraph_sorted <- bedgraph_sorted[sort(paste(bedgraph_sorted$V1,bedgraph_sorted$V2,sep="_"),index.return=T)$ix,]
    bedgraph_sorted <- bedgraph_sorted[sort(bedgraph_sorted$V2,index.return=T)$ix,]
    rownames(bedgraph_sorted) <- 1:nrow(bedgraph_sorted)
  }

  to_merge <- which(bedgraph_sorted$V2[2:nrow(bedgraph_sorted)]<bedgraph_sorted$V3[1:(nrow(bedgraph_sorted)-1)] &
                      bedgraph_sorted$V1[2:nrow(bedgraph_sorted)]==bedgraph_sorted$V1[1:(nrow(bedgraph_sorted)-1)])
}


write.table(bedgraph_sorted,paste("chipseq_tracks/",sample_name,"_mm10_merged.bedGraph",sep=""),col.names = F,row.names = F,sep = "\t",quote = F)


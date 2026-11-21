############################
############################
##                        ## 
##  Additional functions  ##
##                        ##
############################
############################

library(parallel)
library(BiSeq)
library(GenomicRanges)


## getAveLogFC
##################

getAveLogFC <- function(genelist,data_obj,col_gene,col_logfc,func_mask=NULL){
  if (!is.null(func_mask)){
    data_obj <- data_obj[func_mask,]
  }
  logfc <- rep(0,length(genelist))
  print(length(genelist))
  myfun <- function(i){
    if (i %% 1000 == 0){print(i)}
    mask_i <- (data_obj[,col_gene]==genelist[i]) & !is.na(data_obj[,col_gene])
    if(sum(mask_i)>1){
      logfc_i <- mean(data_obj[mask_i,col_logfc])
    } else if (sum(mask_i==1)){
      logfc_i <- data_obj[mask_i,col_logfc]
    } else {
      stop("error")
    }
    return(logfc_i)
  }
  logfc <- mclapply(1:length(genelist),myfun,mc.cores=25)
  logfc <- unlist(logfc)
  return(logfc)
}


## getGeneSymbol 
##################

getGeneSymbol <- function(genelist,data_obj,col_gene,col_genesymbol,func_mask=NULL){
  if (!is.null(func_mask)){
    data_obj <- data_obj[func_mask,]
  }
  logfc <- rep(0,length(genelist))
  print(length(genelist))
  myfun <- function(i){
    if (i %% 1000 == 0){print(i)}
    mask_i <- (data_obj[,col_gene]==genelist[i]) & !is.na(data_obj[,col_gene])
    if(sum(mask_i)>1){
      gene_symbol_i <- unique(data_obj[mask_i,col_genesymbol])
      gene_symbol_i <- gene_symbol_i[gene_symbol_i!=""]
      gene_symbol_i <- paste(gene_symbol_i,collapse=",")
    } else if (sum(mask_i==1)){
      gene_symbol_i <- data_obj[mask_i,col_genesymbol]
      gene_symbol_i <- gene_symbol_i[gene_symbol_i!=""]
      gene_symbol_i <- paste(gene_symbol_i,collapse=",")
    } else {
      stop("error")
    }
    return(gene_symbol_i)
  }
  gene_symbol <- mclapply(1:length(genelist),myfun,mc.cores=25)
  gene_symbol <- unlist(gene_symbol)
  return(gene_symbol)
}


## get orthologues
##################

getOrthologues <- function(dat_obj){
  dat_obj$human_ortholog <- ""
  to_remove <- NULL
  extra_lines <- NULL
  counter <- 0
  for (i in 1:nrow(dat_obj)){
    extra_lines_part <- NULL
    if (sum(orth.mouse_new$mmusculus_homolog_ensembl_gene%in%dat_obj$gene[i])==1){
      dat_obj$human_ortholog[i] <- orth.mouse_new$ensembl_gene_id[orth.mouse_new$mmusculus_homolog_ensembl_gene%in%dat_obj$gene[i]]
    } else if (sum(orth.mouse_new$mmusculus_homolog_ensembl_gene%in%dat_obj$gene[i])>1){
      to_remove <- c(to_remove,i)
      counter <- counter + 1
      for (t in 1:sum(orth.mouse_new$mmusculus_homolog_ensembl_gene%in%dat_obj$gene[i])){
        extra_lines_part <- rbind(extra_lines_part,dat_obj[i,])
      }
      extra_lines_part$human_ortholog <- orth.mouse_new$ensembl_gene_id[orth.mouse_new$mmusculus_homolog_ensembl_gene%in%dat_obj$gene[i]]
      extra_lines <- rbind(extra_lines,extra_lines_part)
    }
  }
  if(counter!=0){
    print(dim(dat_obj))
    counter
    dat_obj <- dat_obj[-to_remove,] 
    dat_obj <- rbind(dat_obj,extra_lines)
  }
  print(dim(dat_obj))
  head(dat_obj)
  sum(dat_obj$human_ortholog!="")
  return(dat_obj)
}


## overlapSets
##################

overlapSets <- function(dat_obj_human,dat_obj_mouse,source_obj_human,source_obj_mouse){
  dat_obj_human_both <- dat_obj_human[dat_obj_human$gene%in%dat_obj_mouse$human_ortholog,]
  dat_obj_mouse_both <- dat_obj_mouse[dat_obj_mouse$human_ortholog%in%dat_obj_human$gene,]
  print(dim(dat_obj_human_both))
  print(dim(dat_obj_mouse_both))
  
  dat_obj_human_both_sorted <- dat_obj_human_both[sort(dat_obj_human_both$gene,index.return=T)$ix,]
  dat_obj_mouse_both_sorted <- dat_obj_mouse_both[sort(dat_obj_mouse_both$human_ortholog,index.return=T)$ix,]
  
  dat_res <-  cbind(dat_obj_human_both_sorted,dat_obj_mouse_both_sorted)
  colnames(dat_res) <- c(paste(rep("lyko_",3),colnames(dat_obj_human_both_sorted),sep=""),
                         paste(rep("roth_",4),colnames(dat_obj_mouse_both_sorted),sep=""))
  return(dat_res)
}


## overlapSetsGS
##################

# dat_obj_human <- DE_genes_lyko[DE_genes_lyko$gene_symbol!="",]
# dat_obj_mouse <- DE_genes_rrbs[DE_genes_rrbs$gene_symbol!="",]

overlapSetsGS <- function(dat_obj_human,dat_obj_mouse,source_obj_human,source_obj_mouse){
  dat_obj_human <- dat_obj_human[dat_obj_human$gene_symbol!="",]
  dat_obj_mouse <- dat_obj_mouse[dat_obj_mouse$gene_symbol!="",]
  
  mat_human <- toupper(do.call(rbind,strsplit(dat_obj_human$gene_symbol,",")))
  mat_mouse <- toupper(do.call(rbind,strsplit(dat_obj_mouse$gene_symbol,",")))
  
  print(nrow(mat_mouse))
  print(ncol(mat_mouse))
  dat_res <- NULL
  for (i in 1:nrow(mat_mouse)){
    # print(i)
    GOR <- unique(mat_mouse[i,])
    # print(GOR)
    for (j in 1:length(GOR)){
      # print(j)
      # print(GOR[j])
      # print(sum(rowSums(mat_human==mat_mouse[i,j])>0))
      if (sum(rowSums(mat_human==GOR[j])>0)==1){
        dat_res_part <- c(dat_obj_human[rowSums(mat_human==GOR[j])>0,],dat_obj_mouse[i,])
        # print(dat_res_part)
        names(dat_res_part) <- c(paste(rep("lyko_",3),colnames(dat_obj_human),sep=""),
                                 paste(rep("roth_",4),colnames(dat_obj_mouse),sep=""))
        # print(dat_res_part)
        # print(length(dat_res_part))
        dat_res <- rbind(dat_res,dat_res_part)
      } else if (sum(rowSums(mat_human==mat_mouse[i,j])>0)>1){
        # print(cbind(dat_obj_human[rowSums(mat_human==mat_mouse[i,j])>0,],t(matrix(rep(dat_obj_mouse[i,],sum(rowSums(mat_human==mat_mouse[i,j])>0)),nrow=ncol(dat_obj_mouse)))))
        dat_res_part <- cbind(dat_obj_human[rowSums(mat_human==GOR[j])>0,],
                              t(matrix(rep(dat_obj_mouse[i,],sum(rowSums(mat_human==GOR[j])>0)),nrow=ncol(dat_obj_mouse))))
        # print(dat_res_part)
        colnames(dat_res_part) <- c(paste(rep("lyko_",3),colnames(dat_obj_human),sep=""),
                                    paste(rep("roth_",4),colnames(dat_obj_mouse),sep=""))
        # print(dat_res_part)
        dat_res <- rbind(dat_res,dat_res_part)
      }
    }
    # print(dim(dat_res))
  }
  print(dim(dat_res))
  dat_res <- data.frame(dat_res)
  for (i in 1:ncol(dat_res)){
    dat_res[,i] <- unlist(dat_res[,i])
  }
  return(dat_res)
}


## makeBedGraph
##################

# dataobj : object of which you want to make a BED file
#
# cn : number of column that contains respectively
# 1) chromosome name
# 2) start pos
# 3) end pos
# 4) score value (if not scaled to 1000 add your interval as scale=c(min,max))
#
# trackname : name of your UCSC track
#
# col : color you want your track to have (use lower case, english names e.g. "blue")
#
# range: if you want the graph in you track to show a fixed y range (e.g. c(0,1))
#
# POI : probe of interest (probe of interest will be plotted in a different color)
#
# startview : region you want to show in UCSC when you open the track (defaults to chr17:914,000-916,000); enter as string "chr17:914,000-916,000"

makeBedGraph <- function(dataobj,cn,trackname,description,col=NULL,range=NULL,startview=NULL,isSign=NULL){
  if (!is.null(isSign)){
    if (!is.null(range)){
      range[1] <- -range[2] 
    }
    dataobj[!isSign,cn[4]] <- -dataobj[!isSign,cn[4]]
  }
  
  if (is.null(startview)){
    write.table("browser position chr17:914,000-916,000",file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F)
  } else {
    write.table(paste("browser position ",startview,sep=""),file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F)
  }
  
  write.table("browser hide all",file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F, append = T)
  
  if (is.null(range)){
    header = paste('track type=bedGraph name=',trackname,' description=',description,
                   ' visibility=full color=',paste(col2rgb(col[1]),collapse=","),' altColor=',paste(col2rgb(col[2]),collapse=","),' gridDefault=on',sep="")
  } else {
    header = paste('track type=bedGraph name=',trackname,' description=',description,
                   ' visibility=full color=',paste(col2rgb(col[1]),collapse=","),' altColor=',paste(col2rgb(col[2]),collapse=","),
                   ' gridDefault=on viewLimits=',range[1],":",range[2]," autoScale=off",sep="")
  }
  
  write.table(header,file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F, append = T)
  
  write.table(cbind(paste("chr",dataobj[,cn[1]],sep=""),dataobj[,cn[2:4]]),file=paste(trackname,".BedGraph",sep=""),
              col.names = F, quote = F, row.names = F,sep="\t", append = T)
  
}


## DMRanalysis
##################


DMRanalysis <- function (infdata_func,groups_func,prom_gr,genic_gr){
  
  probe.features_infdata_all <- probe.features[rownames(probe.features)%in%rownames(infdata_func),]
  infdata_all_nona <- infdata_func[!is.na(probe.features_infdata_all$MAPINFO) & (probe.features_infdata_all$CHR!=""),]
  probe.features_infdata_all_nona <- probe.features_infdata_all[!is.na(probe.features_infdata_all$MAPINFO) & (probe.features_infdata_all$CHR!=""),]
  probe.features_infdata_all_nona$CHR <- as.character(probe.features_infdata_all_nona$CHR)
  
  BS_raw_infdata_all <- BSraw(metadata = NULL,
                              rowRanges = GRanges(seqnames=probe.features_infdata_all_nona$CHR,
                                                  ranges=probe.features_infdata_all_nona$MAPINFO,
                                                  strand="*"),
                              colData = colnames(infdata_all_nona),
                              totalReads = methylated(infdata_all_nona)+unmethylated(infdata_all_nona),
                              methReads = methylated(infdata_all_nona))
  
  rrbs.clust<- clusterSites(object = BS_raw_infdata_all, groups = groups_func, perc.samples = 0.75,
                            min.sites = 15, max.dist = dist_max, mc.cores = 20)
  clusterSitesToGR(rrbs.clust)
  predictedMeth <- predictMeth(object = rrbs.clust, mc.cores = 20)
  betaResults <- BiSeq::betaRegression(formula = ~groups_func, link = "probit", object = predictedMeth, 
                                       type = "BR", mc.cores = 20)
  vario <- makeVariogram(betaResults)
  vario.sm <- smoothVariogram(vario, sill = 0.9)
  locCor <- estLocCor(vario.sm)
  clusters.rej <- testClusters(locCor,FDR.cluster = 0.1)
  print(head(vario.sm$pValsList[[1]]))
  clusters.trimmed <- trimClusters(clusters.rej, FDR.loc = 0.05)
  DMRs <- findDMRs(clusters.trimmed, max.dist = dist_max, diff.dir = TRUE)
  print(head(DMRs))
  
  cluster_tr <- GRanges(data.frame(chr=clusters.trimmed$chr,
                                   start=clusters.trimmed$pos,
                                   end=clusters.trimmed$pos))
  
  overlaps_cpg <- findOverlaps(DMRs,cluster_tr)
  DMRs@elementMetadata$am_CpG <- as.numeric(table(queryHits(overlaps_cpg)))
  
  print(head(DMRs))
  DMRs_sorted <- DMRs[sort(DMRs$median.p,index.return=T)$ix,]
  DMRs_sorted <- data.frame(DMRs_sorted,stringsAsFactors=F)
  DMRs_sorted$seqnames <- as.character(DMRs_sorted$seqnames)
  
  print(head(DMRs_sorted))
  print(dim(DMRs_sorted))
  
  # Overlap
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
  
  DMRs_sorted_annot$FDR <- p.adjust(DMRs_sorted_annot$median.p)
  
  DMRs_control_vs_cancer <- DMRs_sorted_annot
}

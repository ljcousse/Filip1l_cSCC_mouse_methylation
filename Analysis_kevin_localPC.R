library(yaml)

# Directories are read from environment.yml (copy environment_example.yml and fill in the local paths)
if (!file.exists("environment.yml")){
  stop("environment.yml not found, copy environment_example.yml and fill in the local paths")
}
env <- read_yaml("environment.yml")

dir_work <- env$dir_work
setwd(dir_work)

library(xlsx)
output_rrbs <- read.table("output_rrbs.txt",sep="\t",stringsAsFactors = T,header = T)

hypo <- output_rrbs[output_rrbs$logFC<0,]
hyper <- output_rrbs[output_rrbs$logFC>0,]

output_rrbs$extra <- "ALL"
hypo$extra <- "HYPO"
hyper$extra <- "HYPER"

dat_plot <- rbind(output_rrbs,hypo,hyper)

jpeg("barplot.jpg")
ggplot() + geom_bar(aes(y = adj.P.Val, x = extra, fill = cgi_sort), data = dat_plot,
                    stat="identity",position="fill") + ylab("")
dev.off()

dat_plot_sign <- rbind(output_rrbs[output_rrbs$adj.P.Val<0.05,],hypo[hypo$adj.P.Val<0.05,],hyper[hyper$adj.P.Val<0.05,])

jpeg("barplot_sign.jpg")
ggplot() + geom_bar(aes(y = adj.P.Val, x = extra, fill = cgi_sort), data = dat_plot_sign,
                    stat="identity",position="fill")
dev.off()

sum(output_rrbs$adj.P.Val<0.05)
length(unique(output_rrbs$gene[output_rrbs$adj.P.Val<0.05]))-1
sum(is.na(output_rrbs$gene[output_rrbs$adj.P.Val<0.05]))

map_stat <- read.xlsx("Pooling info_OXRRBSAlbena.xlsx",2,stringsAsFactors=F)
tab <- map_stat[grep("B",map_stat$Sample),]
tab <- tab[,c(1,3,6,7)]

am_CpG <- read.table("amountCpGs.txt",skip=1,header=F,stringsAsFactors=F)

tab$amCpG <- 0
for (i in 1:nrow(am_CpG)){
  tab$amCpG[tab$Sample==am_CpG[i,1]] <- am_CpG[i,2]
}

write.table(tab[sort(tab$Sample,index.return=T)$ix,],"table_mapstat.txt",col.names = T,row.names = F,quote = F,sep = "\t")


mask_filip1l <- rowSums((output_rrbs[,grepl("gene",colnames(output_rrbs)) & !grepl("gene_symbol",colnames(output_rrbs))]=="ENSMUSG00000043336") & !is.na(output_rrbs[,grepl("gene",colnames(output_rrbs)) & !grepl("gene_symbol",colnames(output_rrbs))]))>=1

filip1l_kevin <- output_rrbs[mask_filip1l,]
write.table(filip1l_kevin,file="Filip1l_Kevin.txt",row.names = F,col.names = T,sep="\t",quote = F)

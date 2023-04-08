# conda create -n bigWigHandling python=2
# conda activate bigWigHandling run yourself before bash
# conda install ucsc-bigwigmerge
# conda install ucsc-bedsort
# conda install ucsc-bigwigtobedgraph
# conda install ucsc-bedgraphtobigwig

cd chipseq_tracks
rm $(ls | grep -v WT)
cd ..

wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
pigz -d mm9ToMm10.over.chain.gz
fetchChromSizes mm10 > mm10.chrom.sizes

for file in $(ls chipseq_tracks | grep merged.bw)
do
	echo ${file:0:10}
	echo $file
	bigWigToBedGraph chipseq_tracks/$file chipseq_tracks/${file:0:10}.bedGraph
	liftOver chipseq_tracks/${file:0:10}.BedGraph mm9ToMm10.over.chain chipseq_tracks/${file:0:10}_mm10.bedGraph chipseq_tracks/${file:0:10}_failed.bedGraph
	bedSort chipseq_tracks/${file:0:10}_mm10.bedGraph chipseq_tracks/${file:0:10}_mm10_sorted.bedGraph
	#bedGraphToBigWig chipseq_tracks/${file:0:10}_mm10_sorted.BedGraph  mm10.chrom.sizes chipseq_tracks/${file:0:10}_mm10.bw
	Rscript mergeRegions.R ${file:0:10}
	bedSort chipseq_tracks/${file:0:10}_mm10_merged.bedGraph chipseq_tracks/${file:0:10}_mm10_merged_sorted.bedGraph
	bedGraphToBigWig chipseq_tracks/${file:0:10}_mm10_merged_sorted.bedGraph  mm10.chrom.sizes chipseq_tracks/${file:0:10}_mm10.bw
done

bigWigMerge chipseq_tracks/GSM3259343_mm10.bw chipseq_tracks/GSM3259344_mm10.bw chipseq_tracks/GSM3259345_mm10.bw chipseq_tracks/CTCF_merged.bedGraph
bedSort chipseq_tracks/CTCF_merged.bedGraph chipseq_tracks/CTCF_merged_sorted.bedGraph
bedGraphToBigWig chipseq_tracks/CTCF_merged_sorted.bedGraph  mm10.chrom.sizes chipseq_tracks/CTCF.bw
bigWigMerge chipseq_tracks/GSM3259346_mm10.bw chipseq_tracks/GSM3259347_mm10.bw chipseq_tracks/GSM3259348_mm10.bw chipseq_tracks/H3K27Ac_merged.bedGraph
bedSort chipseq_tracks/H3K27Ac_merged.bedGraph chipseq_tracks/H3K27Ac_merged_sorted.bedGraph
bedGraphToBigWig chipseq_tracks/H3K27Ac_merged_sorted.bedGraph  mm10.chrom.sizes chipseq_tracks/H3K27Ac.bw
bigWigMerge chipseq_tracks/GSM3259349_mm10.bw chipseq_tracks/GSM3259350_mm10.bw chipseq_tracks/GSM3259351_mm10.bw chipseq_tracks/H3K27Me3_merged.bedGraph
bedSort chipseq_tracks/H3K27Me3_merged.bedGraph chipseq_tracks/H3K27Me3_merged_sorted.bedGraph
bedGraphToBigWig chipseq_tracks/H3K27Me3_merged_sorted.bedGraph  mm10.chrom.sizes chipseq_tracks/H3K27Me3.bw

#	bigWigToBedGraph chipseq_tracks/GSM3259343_P50-WT-Dorsal-Dermis-CTCF-0-merged.bw chipseq_tracks/GSM3259343.BedGraph
#	liftOver chipseq_tracks/GSM3259343.BedGraph mm9ToMm10.over.chain chipseq_tracks/GSM3259343_mm10.BedGraph chipseq_tracks/GSM3259343_failed.BedGraph
#	bedSort chipseq_tracks/GSM3259343_mm10.BedGraph chipseq_tracks/GSM3259343_mm10_sorted.BedGraph
#	Rscript mergeRegions.R GSM3259343
#	bedSort chipseq_tracks/GSM3259343_mm10_merged.BedGraph chipseq_tracks/GSM3259343_mm10_merged_sorted.BedGraph
#	bedGraphToBigWig chipseq_tracks/GSM3259343_mm10_merged_sorted.BedGraph  mm10.chrom.sizes chipseq_tracks/GSM3259343_mm10.bw
#	
#	bedSort chipseq_tracks/GSM3259343_mm10_merged.BedGraph chipseq_tracks/GSM3259343_mm10_merged_sorted.BedGraph
#	bedGraphToBigWig chipseq_tracks/GSM3259343_mm10_sorted.BedGraph  mm10.chrom.sizes chipseq_tracks/GSM3259343_mm10.bw
#/usr/bin/bash

dir=$1
cd $dir

##library complexity
cd $dir/align
cut -f 3 $dir/meta/macs2.meta  | while read i; do
bedtools bamtobed -i "$i".final.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' |  \grep -v 'chrM' | \sort |  \uniq -c |  \awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > "$i".qc
done

##rep0_FRiP
mkdir $dir/qc
cd $dir/qc
#bam2bed
cut -f 1 $dir/meta/macs2.meta | while read j; do
bedtools bamtobed -i $dir/align/"$j"_rep0.bam > "$j".final.bed
Peakbed=$dir/final/"$j".peaks.bed
BAMbed="$j".final.bed
#peak_overlap_BAM_with_bed
overlapReads=$(bedtools intersect -a $BAMbed -b $Peakbed | wc -l)
#BAM_reads
BAMreads=$(cat $BAMbed | wc -l)
echo -e "$BAMbed\t$overlapReads\t$BAMreads" 
done | tee $dir/FRiP

#bam2bw
cut -f 3 $dir/meta/macs2.meta | while read i; do
bamCoverage -b "$i".final.bam -o "$i".bw
done

#spp
cut -f 3 $dir/meta/macs2.meta | while read Trmt ; do
Rscript /disk2/users/che/phantompeakqualtools-1.2/run_spp.R -c=$Trmt.final.bam -filtchr="chrUn" -savp -out=$Trmt.spp
done

##annotation_deeptools
cut -f 1 ../meta/macs2.meta | uniq | while read i;
do
bedtools closest -D ref -t all -mdb all  -a "$i".peaks.bed -b geneR1.bed >new/"$i".txt
awk '{if($16=="-") print $1"\t"$14"\t"$16"\t"$4"\t"0-$17"\t"$13-$2-$10}' new/"$i".txt | uniq  >new/minus.txt
awk '{if($16=="+") print $1"\t"$14"\t"$16"\t"$4"\t"0-$17"\t"$2+$10-$12}' new/"$i".txt | uniq >new/plus.txt
cat minus.txt plus.txt >"$i"_gene.peaks.txt
done

#plot
Rscript PeakAnnotation.r 

#TSS enrichment
computeMatrix scale-regions -p 10 -S *bw -R /public/home/chaohe/db/geneR1.bed -b 3000 -a 3000 -m 5000 --skipZeros -o chip.mat.gz
plotProfile --dpi 720 -m chip.mat.gz -out chip.profile.pdf --plotFileFormat pdf --perGroup
plotHeatmap -m chip.mat.gz -out chip.merge.png

##correlation
multiBigwigSummary  bins --bwfiles "$dir"/*.bw --binSize 150 --numberOfProcessors 10 --outRawCounts A150.txt -o A150.npz -p 8
plotCorrelation -in A150.npz  \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Read Counts" \--whatToPlot heatmap --colorMap RdYlBu --plotNumbers --removeOutliers \
--plotFileFormat pdf \
-o heatmap_PearsonCorr_readCounts.pdf   \
--outFileCorMatrix PearsonCorr_readCounts.tab


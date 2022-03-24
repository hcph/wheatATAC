#/usr/bin/bash

dir=$1
trimmomatic_dir=$2
db=$3

mkdir $dir/meta
cd $dir/meta
cut -f 3 $dir/meta/ATAC.txt | sort | uniq > ../meta/ATAC.SRX.list
awk -vOFS="\t" '{a[$2]=$3;if(!($2~/_control/)){b[$2]=$1;}} \
END{for(i in b){print b[i],i,a[i],a[i"_control"]}}' \
$dir/meta/ATAC.txt | sort -k 1 -k 2 \
> $dir/meta/ATAC.macs2.meta

#trim
mkdir $dir/fq
mv *fq.gz $dir/fq
cd $dir/fq
cut -f 1 ../meta/ATAC-seq.txt | while read i; do
trimmomatic PE -phred33 "$i"_R1.fq.gz "$i"_R2.fq.gz  "$i"_R1_paired.fq.gz "$i"_R1_unpaired.fq.gz "$i"_R2_paired.fq.gz "$i"_R2_unpaired.fq.gz ILLUMINACLIP:$trimmomatic_dir/NexteraPE-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>trim.log
done

#align
mkdir ../align
cd ../align
cut -f 1 ../meta/ATAC-seq.txt | while read i; do
bowtie2 --no-unal --threads 10 --sensitive -k 3 -q --phred33 --rg-id '"$i"_R1_"$i"_R2' --rg 'SM:"$i"_R1_"$i"_R2\tPL:Illumina\tLB:Illumina_1_8' -x $db/IWGSC_v1 -1 ../fq/"$i"_R1_paired.fq.gz -2 ../fq/"$i"_R2_paired.fq.gz -S "$i".sam 2> "$i"_R1_"$i"_R2.log | samtools sort -@ 4 -m 8g -o "$i".sort.bam
done

#filter
cd ../align
cut -f 1 ../meta/ATAC-seq.txt | while read i; do
samtools sort -@ 4 -m 8g -o "$i".sort.bam "$i".bam
samtools index "$i".sort.bam
samtools flagstat "$i".sort.bam >"$i".sort.flagstat.qc
#samtools过滤未必对上&MAPQ<5的reads
samtools view -F 1804 -q 1 -b "$i".sort.bam >"$i".filter.bam
samtools index "$i".filter.bam
picard MarkDuplicates \
      I="$i".filter.bam \
      O="$i".dupmark.bam \
   REMOVE_DUPLICATES=true \
      M=marked_dup_metrics.txt   
samtools view -@ 8 -F 1804 -b "$i".dupmark.bam > "$i".final.bam
samtools index -c "$i".final.bam
rm "$i".dupmark.bam
done

#prepare for files using idr
#self-pseudo-replicates
cd $dir/align
awk '(!/_control/){print $2, $3}' $dir/meta/macs2.meta \
| while read Rep SRX; do
[[ ! -f $Rep.pr1.bam || ! -f $Rep.pr2.bam ]] && \
NN=$(samtools flagstat $SRX.final.bam | awk '(/mapped \(/){printf("%d\n",$1/2+0.9)}'); \
samtools view -H $SRX.final.bam > $Rep.header; \
samtools view $SRX.final.bam | shuf | split -d -l $NN - ${Rep}; \
cat $Rep.header ${Rep}00 | samtools view -bS - \
| samtools sort -o ${Rep}.pr1.bam & \
cat $Rep.header ${Rep}01 | samtools view -bS - \
| samtools sort -o ${Rep}.pr2.bam; \
samtools index ${Rep}.pr1.bam & \
samtools index ${Rep}.pr2.bam; \
rm -f ${Rep}0* $Rep.header
done

#merge replicates into a pooled sample
awk -vFS="\t" -vOFS="\t" '{bam=$3".final.bam"; \
if(($2~/_control/)){c=$1"_control"}else
{c=$1"_rep0"} \
a[c]=a[c]?a[c]";"bam:bam} END{for(i in a){print i,a[i]}}' \
$dir/meta/macs2.meta \
| while read TF bams; do
out=$TF".bam"
if [[ $bams =~ ";" ]]; then
inbams=$(echo $bams | sed -e 's/;/ /g')
[[ ! -f $out || ! -f "$out.bai" ]] && \
samtools merge $out $inbams; samtools index $out
else
ln -sf $bams $out && ln -sf $bams".bai" $out".bai";
fi
done

#pooled pseudo-replicates
cut -f 1 $dir/meta/macs2.meta | uniq | while read TF; do
BAM=$dir/align/"$TF"_rep0.bam
[[ ! -f $BAM ]] && continue;
OUT="$TF"_rep0
[[ ! -f $OUT.pr1.bam || ! -f $OUT.pr2.bam ]] && \
NN=$(samtools flagstat $BAM | awk '(/mapped \(/){printf("%d\n",$1/2+0.9)}'); \
samtools view -H $BAM > $OUT.header; \
samtools view $BAM | shuf | split -d -l $NN - ${OUT}; \
cat $OUT.header ${OUT}00 | samtools view -bS - \
| samtools sort -o ${OUT}.pr1.bam & \
cat $OUT.header ${OUT}01 | samtools view -bS - \
| samtools sort -o ${OUT}.pr2.bam; \
samtools index ${OUT}.pr1.bam & \
samtools index ${OUT}.pr2.bam; \
rm -f ${OUT}0* $OUT.header
done

#run MACS2 on individual replicates to call peak
mkdir -p "$dir"/macs2/signal && cd "$dir"/macs2
rm chrom.sizes
cat "$dir"/meta/macs2.meta | while read TF Rep  Trmt Ctrl  j v l m; do
ChIP="$dir"/align/$Trmt.final.bam
Ctrl="$dir"/align/${TF}_control.bam
OUT=${Rep}
[[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
awk -vOFS="\t"  '(/^@SQ/){match($0,/SN:(\w+)/,SN); match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}'> chrom.sizes
GS=$(awk 'BEGIN{GS=0}{GS+=$2}END{print int(0.85*GS)}' chrom.sizes) 
[[ ! -f ${OUT}.peaks.bed ]] && \
macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} -g $GS -p 1e-2 --mfold 2 20 --nomodel --to-large; 
maxS=$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | head -n 1 | cut -f 5); 
minS=$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | tail -n 1 | cut -f 5); 
awk -vOFS='\t' -vm=$minS -vM=$maxS '{$5=int((($5-m)*(1000-10)/(M-m))+10); print}' ${OUT}_peaks.narrowPeak > ${OUT}.peaks.bed; \
bedToBigBed -type=bed6+4 ${OUT}.peaks.bed chrom.sizes signal/${OUT}.peaks.bb; 
done

#run MACS2 on pseudo replicates (either self-pseudo-replicates or pooled-pseudo-replicates)
rm chrom.sizes
cat "$dir"/meta/macs2.meta |uniq | while read TF Rep Trmt Ctrl j v l m; do
Ctrl="$dir"/align/${TF}_control.bam
for pRep in pr1 pr2; do
ChIP="$dir"/align/$Rep.$pRep.bam
OUT=${Rep}_self_${pRep}
[[ ! -f chrom.sizes ]] && samtools view -H $ChIP | awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' > chrom.sizes
GS=$(awk 'BEGIN{GS=0}{GS+=$2}END{print int(0.85*GS)}' chrom.sizes) 
[[ ! -f ${OUT}.peaks.bed ]] && \
macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} -g $GS -p 1e-2  --mfold 2 20 --nomodel --to-large;  \
done
done

# peak calling on pooled pseudo-replicates
rm chrom.sizes
cat "$dir"/meta/macs2.meta | while read TF Rep pRep j v l; do
Ctrl="$dir"/align/${TF}_control.bam
for pRep in pr1 pr2; do
ChIP="$dir"/align/${TF}_rep0.${pRep}.bam;
[[ ! -f $ChIP ]] && continue;
OUT=${TF}_pooled_${pRep}
[[ ! -f chrom.sizes ]] && samtools view -H $ChIP |\
awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' \
> chrom.sizes
GS=$(awk 'BEGIN{GS=0}{GS+=$2}END{print int(0.85*GS)}' chrom.sizes) 
[[ ! -f ${OUT}.peaks.bed ]] && \
macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} -g $GS -p 1e-2   --mfold 2 20 --nomodel --to-large; \
mv ${OUT}_peaks.narrowPeak ${OUT}.peaks.bed
done
done

#run MACS2 on the pooled samples to call peaks and to generate signal filesin bigWig format
mkdir -p "$dir"/macs2/signal && cd "$dir"/macs2
rm chrom.sizes
cat "$dir"/meta/macs2.meta | uniq | while read TF j v l; do
ChIP="$dir"/align/${TF}_rep0.bam
Ctrl="$dir"/align/${TF}_control.bam
[[ ! -f $ChIP || ! -f $Ctrl ]] && continue;
OUT=${TF}_pooled
[[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' > chrom.sizes;
cSF=$(samtools idxstats $ChIP | awk '!(/*/){TL+=$3} END{print TL/10e+6}');
tSF=$(samtools idxstats $Ctrl | awk '!(/*/){TL+=$3} END{print TL/10e+6}')
SF=$(echo "$cSF $tSF" | awk '($1>$2){print $2} ($1<=$2){print $1}')
GS=$(awk '{GS+=$2}END{print int(0.85*GS)}' chrom.sizes) 
[[ ! -f ${OUT}.peaks.bed ]] && \
macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} -g $GS -p 1e-2   --mfold 2 20 --nomodel    -B --SPMR --to-large; \
maxS=$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | head -n 1 | cut -f 5); \
minS=$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | tail -n 1 | cut -f 5); \
awk -vOFS='\t' -vm=$minS -vM=$maxS '{$5=int((($5-m)*(1000-10)/(M-m))+10); print}' ${OUT}_peaks.narrowPeak > ${OUT}.peaks.bed; \
done

#signal
rm chrom.sizes
cat "$dir"/meta/macs2.meta | uniq | while read TF j v l; do
ChIP="$dir"/align/${TF}_rep0.bam
Ctrl="$dir"/align/${TF}_control.bam
[[ ! -f $ChIP || ! -f $Ctrl ]] && continue;
OUT=${TF}_pooled
[[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' > chrom.sizes;
[[ ! -f ${OUT}.peaks.bed ]] && \
bedToBigBed -type=bed6+4 ${OUT}.peaks.bed chrom.sizes ./signal/${OUT}.peaks.bb; \
macs2 bdgcmp -m FE -t ${OUT}_treat_pileup.bdg -c ${OUT}_control_lambda.bdg -o ./signal/${OUT}_FE.bdg; \
sortBed -i ./signal/${OUT}_FE.bdg | slopBed -i - -g chrom.sizes -b 0 | bedClip stdin chrom.sizes ./signal/${OUT}.fc.signal.bdg; \
bedGraphToBigWig ./signal/${OUT}.fc.signal.bdg chrom.sizes ./signal/${OUT}.fc.signal.bw; \
macs2 bdgcmp -m ppois -t ${OUT}_treat_pileup.bdg -c ${OUT}_control_lambda.bdg -o ./signal/${OUT}_ppois.bdg -S $SF; \
sortBed -i ./signal/${OUT}_ppois.bdg | slopBed -i - -g chrom.sizes -b 0 | bedClip stdin chrom.sizes ./signal/${OUT}.pval.signal.bdg; \
bedGraphToBigWig ./signal/${OUT}.pval.signal.bdg chrom.sizes ./signal/${OUT}.pval.signal.bw; \
rm chrom.sizes
done

#run IDR analysis on original replicates
mkdir -p $dir/macs2/idr && cd $dir/bin/idr-Code/
pDIR=$dir/macs2
awk -vFS="\t" -vOFS="\t" '{a[$1]=a[$1]?a[$1]"\t"$2:$2;} END{for(i in a){print i,a[i]}}' ../meta/macs2.meta \
| sort -k 1 -k 2 | tee ../meta/idr.meta \
| awk -vFS="\t" -vDIR=$pDIR '(NF>2){for(i=2;i<=NF;i++) \
{for(j=i+1;j<=NF;j++){
idr --samples "DIR"/"$i".peaks.bed "DIR"/"$j".peaks.bed --output-file  ../macs2/idr/"$i"_overlapped_peaks.txt --plot}}' 

#Run IDR analysis on self-pseudo-replicates
pDIR=$dir/macs2
cut -f 2 $dir/meta/macs2.meta | while read Rep; do
idr --samples $pDIR/${Rep}_self_pr1.peaks.bed $pDIR/${Rep}_self_pr2.peaks.bed --output-file ../macs2/idr/${Rep}_self_overlapped_peaks.txt --plot
done

#Run IDR analysis on pooled pseudo-replicates
pDIR=$dir/macs2
cut -f 1 $dir/meta/macs2.meta | uniq | while read TF; do
idr --samples $pDIR/${TF}_rep0_pr1.peaks.bed $pDIR/${TF}_rep0_pr2.peaks.bed --output-file $pDIR/idr/${TF}_pooled_overlapped_peaks.txt --plot
done

#filtering and obtain consistent peaks
mkdir -p $dir/final
oDIR=$dir/idr
cut -f 1,2,3 $dir/meta/macs2.meta | while read TF Rep MORE; do
        # self-consistency threshold
        nps_Self=$( awk '$5 >= 540 {print $0}' $oDIR/${Rep}_self_overlapped_peaks.txt | wc -l )
        # original replicate threshold
        max_nps_Rep=
        if [[ $(ls $oDIR/${TF}*${Rep}*_overlapped_peaks.txt \
                2> /dev/null | wc -l) -gt 0 ]]; then
            echo "*** $TF with replicates" 1>&2
            numPeaks=
            for cmp in $oDIR/${TF}*${Rep}*_overlapped_peaks.txt; do
                numPeaks="$numPeaks "$( awk '$5 >= 705 \
                    {print $0}' $cmp | wc -l )
            done
            max_nps_Rep=$(echo $numPeaks | awk \
                '{for(i=1;i<=NF;i++){if(!m||$i>m) \
                {m=$i}}}END{print m}')
        fi
        # pooled-consistency threshold
        nps_Rep0=$( awk '$5 >= 830 {print $0}' $oDIR/${TF}_pooled_overlapped_peaks.txt | wc -l )
        echo -e "$TF\t$Rep\t$nps_Self\t$max_nps_Rep\t$nps_Rep0"
    done | tee $dir/meta/idr.threshold \
    | awk -vFS="\t" -vOFS="\t" '{n[$1]++; \
        s[$1]=s[$1]?s[$1]"|"$3:$3; \
        if(!r[$1]||r[$1]<$4){r[$1]=$4} \
        if(!p[$1]||p[$1]<$5){p[$1]=$5} } \
        END{for(i in n){o=p[i];if(!o||o<r[i]){o=r[i]} \
            print i,o,n[i],s[i],r[i],p[i]}}' \
    | sort -k 1 | tee $dir/meta/idr.optimal.threshold \
    | while read TF optThresh MORE; do
    sort -k8nr,8nr $dir/macs2/${TF}_pooled.peaks.bed \
        | head -n $optThresh | bedtools sort -i \
        > $dir/final/$TF.peaks.bed
done


# wheatATAC
A pipeline integrated  data analysis and valuation for wheat ATAC-seq

###install dependencies

conda create -n wheatATAC

conda activate wheatATAC

conda install -c bioconda macs2 

conda install -c bioconda bowtie2

conda install -c bioconda fastqc

conda install -c bioconda samtools

conda install -c bioconda idr

conda install -c bioconda phantompeakqualtools

conda install -c bioconda trimmomatic

conda install -c bioconda deeptools

conda install -c bioconda bedtools

conda install -c bioconda ucsc-bedtobigbed

conda install r=3.6  ##The verson of R must >=3.6

###2.analysis the ATAC-seq dataset

Usage: sh analysis.sh  work_dir trimmomatic_adapter_dir wheat_reference_dir

notice: The fastq file are in *_R1_fq.gz

##2.1 make meta files

The format of meta files are like following. The first column lists the sample name, the second column are the repeats information of samples, and the last column lists the corresponding prefix of the fastq files. The file we named ATAC.txt here.

#with repeats and without inputs:

sampleA	sampleA_rep1	C1

sampleB	sampleA_rep1	C2

sampleC	sampleA_rep1	C3

...

#with repeats and without inputs:

sampleA	sampleA_rep1	C1

sampleA	sampleA_rep2	C2

sampleB	sampleB_rep1	C3

sampleB	sampleB_rep2	C4

...

#with repeats and with inputs:

sampleA	sampleA_rep1	C1

sampleA	sampleA_rep1_control	N1

sampleA	sampleB_rep2	C2

sampleA	sampleB_rep2_control	N1


##2.2 trim, alignment, remove unmapped and duplicate reads


#note that the adapter file is user definded, please use vi align.sh to replace it.


##2.3 prepare file for IDR analysis


##2.4 peak calling, idr analysis and consistent peaks acquirment



###data valuation

sh valuation.sh

#The process of data valuation contains file convert for IGV visualization, FRiP and SPOT calculation, TSS enrichment analysis, peak annotation, samples correlation analysis and signal comparasion.

#peak annotaion according to peak summit
 assign_peaks_summit_to_genes.sh peaks.bed gene.bed output_idr


[README.txt](https://github.com/HeChao7021/wheatATAC/files/8139677/README.txt)



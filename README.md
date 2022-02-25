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
conda install -c bioconda bedtools
conda install r=3.6  ##The verson of R must >=3.6

###1.make meta files
The format of meta files are like following. The first column lists the sample name, the second column are the repeats information of samples, and the last column lists the corresponding prefix of the fastq files.
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

###2.analysis the ATAC-seq dataset
#2.1 trim, alignment, remove unmapped and duplicate reads










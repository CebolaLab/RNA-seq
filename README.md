# RNA-seq
Step-by-step analysis pipeline for RNA-seq data

Resources include: 

- <https://vallierlab.wixsite.com/pipelines/rna-seq>
- [RNA-seq workflow: gene-level exploratory analysis and differential expression by Love et al. 2019](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#running-the-differential-expression-pipeline)

The pipeline covers the following steps:

- [Pre-alignment quality control (QC)](#pre-alignment-qc)
- [Alignment](#alignment)
- [Post-alignment QC](#post-alignment-qc)


## Pre-alignment QC

The raw sequence data should first be assessed for quality. FastQC reports can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. In ATAC-seq data, it is likely that Nextera sequencing adapters will be over-represented. As described by (Yan et al. 2020), base quality should be high although may drop slightly at the 3' end, while GC content and read length should be consistent with the expected values. For paired-end reads, run fastqc on both files, with the results output to the current directory:

fastqc <sample>_1.fastq.gz -d . -o .

fastqc <sample>_2.fastq.gz -d . -o .

Adapter trimming using cutadapt. Other tools  include trimmomatic and fastp. 

## Alignment

Align to the reference genome. Use splice-aware aligners such as STAR or TopHat. 

Merge, sort and index.

## Post-alignment QC

Remove duplicates? 

## Visualisation 

bam to bedGraph to BigWig.

Biological replicates QC - deeptools plotCorrelation
bedtools, deeptools.
UCSC binaries bedGraph to BigWig?

## Quantification 

Counts matrix generated (genes x samples)
Bam -- featureCounts function in the Rsubread R package OR salmon --> counts 

Salmon output The output is transcript-level abundance data which you will need to import into R using the tximport package.

## Functional analysis 

DEseq2, edgeR, limma for example. 
# RNA-seq
[Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)

Step-by-step analysis pipeline for RNA-seq data

Resources include: 

- <https://vallierlab.wixsite.com/pipelines/rna-seq>
- [RNA-seq workflow: gene-level exploratory analysis and differential expression by Love et al. 2019](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#running-the-differential-expression-pipeline)
- https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html

The pipeline covers the following steps:

- [Pre-alignment quality control (QC)](#pre-alignment-qc)
- [Alignment](#alignment)
- [Post-alignment QC](#post-alignment-qc)

## Pre-alignment QC

#### Generate QC report

The raw sequence data should first be assessed for quality. FastQC reports can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. For paired-end reads, run fastqc on both files, with the results output to the current directory:

```
fastqc <sample>_1.fastq.gz -d . -o .

fastqc <sample>_2.fastq.gz -d . -o .
```

#### Adapter trimming 

If there is evidence of adapter contamination shown in the fastQC report (see below), adapter sequences may need to be trimmed, using a tools such as cutadapt, trimmomatic and fastp. In this pipeline, fastp is used to trim adapters. 

```
fastp -i <sample>_R1.fastq.gz -O <sample>_R1.trimmed.fastq.g -I <sample>_R2.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -l 100 -j <sample>.fastp.json -h <sample>.fastp.html
```

A html report is generated, including the following information:

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/fastp-summary.png" width="600">

## Alignment

The filtered DNA reads must next be aligned to the reference genome. For RNA-seq data, a splice-aware aligner such as STAR or TopHat should be used. Here, [STAR](https://github.com/alexdobin/STAR) is used. The manual is available [here](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). The reference genome used is the GRCh38 'no-alt' assembly from ncbi, recommended by [Heng Li](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use). The genome can be downloaded at [this link](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz).  This version of the recent GRCh38 reference genome excludes alternative contigs which may cause fragments to map in multiple locations. The downloaded genome should be indexed with STAR. Other sources recommend the Ensembl [Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz](ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz), however [Heng Li](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) notes that this version of the genome includes multi-placed sequences such as the pseudo-autosomal regions on both chromosomes Z and Y, as well as some alpha satellites. 

Set --sjdbOverhang to your maximum read length -1. The indexing also requires a file containing gene annotation, which comes in a `gtf` format. For example, ENCODE provides a gtf file with GRCh38 annotations, containing gencode gene coordinates, along with UCSC tRNAs and a PhiX spike-in. This file (ENCFF159KBI) can be downloaded [here](https://www.encodeproject.org/files/ENCFF159KBI/). The user should aim to use the most up-to-date reference files, while ensuring that the format is the same as the reference genome. For example, UCSC uses the 'chr1, chr2, chr3' naming convention, while ENSEMBL uses '1, 2, 3' etc. The files suggested here are compatible. 

```
GENOMEDIR=/path/to/indexed/genome

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $GENOMEDIR --genomeFastaFiles $GENOMEDIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sjdbGTFfile ENCFF159KBI.gtf --sjdbOverhang: readlength -1
```

STAR can then be run to align the fastq data files to the genome:


```
STAR --runThreadN 4 --genomeDir $GENOMEREF --readFilesIn <sample>_R1.fastq.gz <sample>_R2.fastq.gz --outFileNamePrefix <sample> | samtools view -bS sort -c 
```

### Merge files

At this stage, if samples have been sequenced across multiple lanes, the samples files can be combined using `samtools merge`. Various QC tools can be used to assess reproducibility and assess lane effects. 

### Sort and index 

The output `bam` files should be sorted and indexed:

```
picard SortSam I=<sample>.bam O=<sample>-sorted.bam SO=coordinate CREATE_INDEX=TRUE
```


## Post-alignment QC

Check the expected mapping to exons, introns etc.

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
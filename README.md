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

STAR can then be run to align the fastq data files to the genome. If the fastq files are in the compressed `.gz` format, add the `--readFilesCommand zcat` argument. 

```
STAR --runThreadN 4 --genomeDir $GENOMEREF --readFilesIn <sample>_R1.fastq.gz <sample>_R2.fastq.gz --outFileNamePrefix <sample> --outSAMtype BAM SortedByCoordinate
```

### Merge files

At this stage, if samples have been sequenced across multiple lanes, the samples files can be combined using `samtools merge`. Various QC tools can be used to assess reproducibility and assess lane effects, such as `deeptools plotCorrelation`.

## Post-alignment QC

The post-alignment QC steps involve several steps:

- [Remove mitochondrial reads](#remove-mitochondrial-reads)
- [Remove duplicates & low-quality alignments](#tag-and-remove-duplicates-and-low-quality-alignments) (including non-uniquely mapped reads)
- Check the expected mapping to exons, introns etc.

### Remove mitochondrial reads

Remove mitochondrial reads. To assess the total % of mitochondrial reads, samtools idxstats can be run to report the total number of reads mapping to each chromosome. samtools flagstat provides a short report including the total number of DNA fragments.

```
samtools idxstats <sample>_sorted.bam > <sample>_sorted.idxstats

grep "chrM" <sample>_sorted.idxstats
```

The second column is the length of the chromosome and the third column is the total number of reads aligned to the chromosome (chrM). To see the total number of DNA fragments, run:

```
samtools flagstat <sample>_sorted.bam > <sample>_sorted.flagstat

head <sample>_sorted.flagstat
```

The % of DNA fragments aligned to chrM can be calculated as a % of the total DNA fragments. To remove any mitocondrial DNA, run the following:


```
samtools view -h <sample>-sorted.bam | grep -v chrM | samtools sort -O bam -o <sample>.rmChrM.bam -T .
```

### Tag and remove duplicates and low-quality alignments

> Duplicate reads

The next filtering steps include marking and [optionally] removing PCR duplicates, as well as removing low-quality reads (described below). Duplicate reads can be marked and the % of duplicate reads viewed using:

```
picard MarkDuplicates QUIET=true INPUT=<sample>.rmChrM.bam OUTPUT=<sample>.marked.bam METRICS_FILE=<sample>.sorted.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

head -n 8 <sample>-markDup.metrics | cut -f 7,9 | grep -v ^# | tail -n 2
```

*Optional*: It may be recommended to remove duplicate reads if the % of duplicates is high. To remove duplicate reads, run the following code:

```
samtools view -h -b -F 1024 <sample>.marked.bam > <sample>.rmDup.bam
```

> Multi-mapping reads

**Multi-mapped reads**: these are defined as reads which can map in more than one location in the reference genome. A recent publication reviewing the handling of multi-mapping reads in RNA-seq data is reported by [Deschamps-Francoeur et al. (2020)](https://www.sciencedirect.com/science/article/pii/S2001037020303032).

## Visualisation 

The QC-ed `bam`	file can be converted to a `bedGraph` format to	visualise sequencing trakcs using tools	such as	the UCSC browser \
or the integrative genomes browser. The	ENCODE blacklist regions can be	provided, to exclude them from the output:

```
bedtools bamCoverage --blackListFileName --normalizeUsing BPM -b <sample>.filtered.bam > <sample>.bedGraph
```


UCSC binaries bedGraph to BigWig?

## Quantification 

Counts matrix generated (genes x samples)
Bam -- featureCounts function in the Rsubread R package OR salmon --> counts 

Salmon output The output is transcript-level abundance data which you will need to import into R using the tximport package.

## Functional analysis 

DEseq2, edgeR, limma for example. 
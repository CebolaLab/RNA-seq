# RNA-seq

## **UNDER CONSTRUCTION**

[Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)

Step-by-step analysis pipeline for RNA-seq data

Resources include: 

- <https://vallierlab.wixsite.com/pipelines/rna-seq>
- [RNA-seq workflow: gene-level exploratory analysis and differential expression by Love et al. 2019](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#running-the-differential-expression-pipeline)
- https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
- The Encode pipeline for long-RNAs: https://www.encodeproject.org/data-standards/rna-seq/long-rnas/

Understanding normalisation:
- https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html


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

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/fastp-summary.png" width="700">

For this RNA-seq pipeline, the steps include:

- [Align to the reference human genome](#align-to-the-reference-genome)
- [Post-alignment QC](#post-alignment-qc)
- [Visualise tracks against the reference genome](#visualisation)
- [Quantify transcripts](#quantification)


## Align to the reference genome

The raw RNA-seq data in `fastq` format will be aligned to the reference genome, along with a reference transcriptome, to output two alignment files: the genome alignment and the transcriptome alignemnt. 

The DNA reads are aligned using the splice-aware aligner, STAR. Here, [STAR](https://github.com/alexdobin/STAR) is used. The manual is available [here](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). The reference genome used is the GRCh38 'no-alt' assembly from ncbi, recommended by [Heng Li](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use). The genome can be downloaded at [this link](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz).  This version of the recent GRCh38 reference genome excludes alternative contigs which may cause fragments to map in multiple locations. The downloaded genome should be indexed with STAR. Other sources recommend the Ensembl [Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz](ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz), however [Heng Li](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) notes that this version of the genome includes multi-placed sequences such as the pseudo-autosomal regions on both chromosomes Z and Y, as well as some alpha satellites. 

> Index the reference genome

Set --sjdbOverhang to your maximum read length -1. The indexing also requires a file containing gene annotation, which comes in a `gtf` format. For example, ENCODE provides a gtf file with GRCh38 annotations, containing gencode gene coordinates, along with UCSC tRNAs and a PhiX spike-in. This file (ENCFF159KBI) can be downloaded [here](https://www.encodeproject.org/files/ENCFF159KBI/). The user should aim to use the most up-to-date reference files, while ensuring that the format is the same as the reference genome. For example, UCSC uses the 'chr1, chr2, chr3' naming convention, while ENSEMBL uses '1, 2, 3' etc. The files suggested here are compatible. 

```
GENOMEDIR=/path/to/indexed/genome

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $GENOMEDIR --genomeFastaFiles $GENOMEDIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sjdbGTFfile ENCFF159KBI.gtf --sjdbOverhang readlength -1
```

> Carry out the alignment

STAR can then be run to align the `fastq` raw data to the genome. If the fastq files are in the compressed `.gz` format, the `--readFilesCommand zcat` argument is added. The output file should be unsorted, as required for the downstream quantification step using Salmon. The following options are shown according to the ENCODE recommendations. For single-end data:

```
STAR --runThreadN 4 --genomeDir $GENOMEREF --readFilesIn <sample>.fastq.gz --outFileNamePrefix <sample> --readFilesCommand zcat --outSAMtype BAM Unsorted --quantTranscriptomeBan Singleend
--outFilterType BySJout --alignSJoverhangMin 8 --outFilterMultimapNmax 20
--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999
--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000
--alignMatesGapMax 1000000 
--quantMode TranscriptomeSAM --outSAMattributes NH HI AS NM MD
```

For compatibility with the STAR quantification, the `--quantMode TranscriptomeSAM` option will result in the output of two alignment files, one to the reference genome (`Aligned.*.sam/bam`) and one to the transcriptome (`Aligned.toTranscriptome.out.bam`).

#### Merge files [optional]

At this stage, if samples have been sequenced across multiple lanes, the sample files can be combined using `samtools merge`. Various QC tools can be used to assess reproducibility and assess lane effects, such as `deeptools plotCorrelation`. The `salmon` quantification does not require files to be merged, since multiple `bam` files can be listed in the command. However, to visualise the RNA-seq data from the combined technical replicates, `bam` files can be merged at this stage. For example, if your sample was split across lanes 1, 2 and 3 (`L001`, `L002`, `L003`):

```
samtools merge <sample>-merged.bam <sample>_L001.bam <sample>_L002.bam <sample>_L003.bam
```

## Post-alignment QC

First, QC reports will be generated using [qualimap](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html). Run on either the `<sample>.gzAligned.out.bam` or `<sample>.merged.bam`.

```
samtools sort <sample>.merged.bam > <sample>-sorted.bam

qualimap bamqc -bam <sample>-sorted.bam -gtf ENCFF159KBI.gtf -outdir <sample>-bamqc-qualimap-report --java-mem-size=16G

qualimap rnaseq -bam <sample>-sorted.bam -gff ENCFF159KBI.gtf -outdir <sample>-rnaseq-qualimap-reports --java-mem-size=16G
```

Qualimap can then run QC on combined samples. One included analysis is principal component analysis, which clusters the samples. This can be used to confirm whether technical and biological replicates cluster together. A text file should be created with 

Note, some versions of qualimap require the raw_data_qualimapReport directory to be renamed to raw_data.

```
qualimap multi-bamqc sample.txt
```

The QC reports can be combined using [multiqc](https://multiqc.info/); an excellent tool for combining QC reports of multiple samples into one. Example outputs of qualimap/multiqc include the alignment positions 

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/multiqc-alignment.png" width="800">

### Remove duplicates?

It is generally recommended to *not* remove duplicates when working with RNA-seq data, unless using UMIs (unique molecular identifiers) [Klepikova et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5357343/). This is because there are likely to be DNA molecules which are natural duplicates of each other, for example originating from genes with a shared sequence in a common domain. Typically, removing duplicates does more harm than good. It is more or less impossible to remove duplicates from single-end data and research has also suggested it may cause false negatives when applied to paired end data. See more in [this useful blog post](https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/). Generally, duplicates are not a problem so long as the *library complexity is high*. 

### Compute GC bias

GC-bias describes the bias in sequencing depth depending on the GC-content of the DNA sequence. Bias in DNA fragments, due to the GC-content and start-and-end sequences, may be increased due to preferential PCR amplification [(Benjamini and Speed, 2012)](https://academic.oup.com/nar/article/40/10/e72/2411059). A high rate of PCR duplications, for example when library complexity is low, may cause a significant GC-bias due to the preferential amplification of specific DNA fragments. This can significantly impact transcript abundance estimates. Bias in RNA-seq is explained in a handy [blog](https://mikelove.wordpress.com/2016/09/26/rna-seq-fragment-sequence-bias/) and [video](https://youtu.be/9xskajkNJwg) by Mike Love.

It is **crucial** to correct GC-bias when comparing groups of samples which may have variable GC content dependence, for example when samples were processed in different libraries. `Salmon`, used later to generate read counts for quantification, has its own in-built method to correct for GC-bias. When generating `bedGraph/BigWig` files for visualisation, the user may opt to correct GC-bias so that coverage is corrected and appears more uniform. The `deeptools` suite includes tools to calculate GC bias and correct for it.

The reference genome file should be converted to `.2bit` format using [`faToTwoBit`](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit).
The effective genome size can be calculated using `faCount` available [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/).
Set the `-l` argument to your fragment length.

The input `bam` file requires an index, which can be generated using `samtools index`.

```
deeptools computeGCBias -b <sample>-sorted.bam --effectiveGenomeSize 3099922541 -g GCA_000001405.15_GRCh38_no_alt_analysis_set.2bit -l 100 --GCbiasFrequenciesFile <sample>.freq.txt  --biasPlot <sample>.biasPlot.pdf
```

The bias plot format can be changed to png, eps, plotly or svg. If there is significant evidence of a GC bias, this can be corrected using `correctGCbias`. An example of GC bias can be seen in the plot outout from `computeGCBias` below:

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/GCbiasPlot.png" width="500">

If opting to correct the GC-bias, `correctGCbias` can be used. This tool effectively removes reads from regions with greater-than-expected coverage (GC-rich regions) and removes reads from regions with less-than-expected coverage (AT-rich regions). The methods are described by [Benjamini and Speed [2012]](https://academic.oup.com/nar/article/40/10/e72/2411059). The following code can be used:

```
 correctGCBias -b <sample>-sorted.bam --effectiveGenomeSize 3099922541 -g GCA_000001405.15_GRCh38_no_alt_analysis_set.2bit --GCbiasFrequenciesFile <sample>.freq.txt -o <sample>.gc_corrected.bam [options]
```

**NOTE:** When calculating the GC-bias for ChIP-seq, ATAC-seq, DNase-seq (and CUT&Tag/CUT&Run) it is recommended to filter out problematic regions. These include those with low mappability and high numbers of repeats. The compiled list of [ENCODE blacklist regions](https://www.nature.com/articles/s41598-019-45839-z) should be excluded. However, the ENCODE blacklist regions have little overlap with coding regions and this step is not necessary for RNAs-seq data [(Amemiya et al, 2019)](https://www.nature.com/articles/s41598-019-45839-z).

### Check correlation of technical and biological replicates

The correlation between `bam` files of biological and technical replicates can be calculated as a QC step to ensure that the expected replicates positively correlate. 

```
multiBamSummary

plotCorrelation
```


## Visualisation 

The `bam` file aligned to the *genome* should be converted to a `bigWig` format, which can be uploaded to genome browsers and viewed as a track. The gene counts are here normalised to TPM values during conversion. 

```
bamCoverage -b <sample>.gc_corrected.bam -o <sample>.bw --normalizeUsing BPM --samFlagExclude 512
```

Biological replicates can also be merged and a bigWig file generated for the combined sample.

```
bamCompare -b <sample>_1.bam
```

There are multiple methods available for normalisation. Recent analysis by [Abrams et al. (2019)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3247-x#Sec2) advocated TPM as the most effective method. 


## Quantification

The `bam` file previously aligned to the *transcriptome* by STAR will next be input into [Salmon](https://combine-lab.github.io/salmon/) in alignment-mode, in order to generate a matrix of gene counts. The Salmon documentation is available [here](https://salmon.readthedocs.io/en/latest/).

Salmon uses a VBEM algorithm 


Salmon is here used with the expectation minimisation (EM) approach method for quantification. This is described in the 2020 paper by [Deschamps-Francoeur et al.](https://www.sciencedirect.com/science/article/pii/S2001037020303032), which describes the handling of multi-mapped reads in RNA-seq data. Duplicated sequences such as pseudogenes can cause reads to align to multiple positions in the genome. Where transcripts have exons which are similar to other genomic sequences, the EM approach attributes reads to the most likely transcript. 

```
salmon quant --useEM -t gencode.v35.transcripts.fa --libType A -a <sample>.gzAligned.toTranscriptome.out.bam -o <sample>.salmon_quant 
```

If using single end data, add the `--fldMean` and `--fldSD` parameters to include the mean and standard deviation of the fragment lengths. 


## Functional analysis 

DEseq2, edgeR, limma for example. 

## Functional analysis


**Preseq**: Estimates library complexity

**Picard RNAseqMetrics**: Number of reads that align to coding, intronic, UTR, intergenic, ribosomal regions, normalize gene coverage across a meta-gene body, identify 5’ or 3’ bias

**RSeQC**: Suite of tools to assess various post-alignment quality, Calculate distribution of Insert Size, Junction Annotation (% Known, % Novel read spanning splice junctions), BAM to BigWig (Visual Inspection with IGV)


The post-alignment QC steps involve several steps:

- [Remove mitochondrial reads](#remove-mitochondrial-reads)
- [Remove duplicates & low-quality alignments](#tag-and-remove-duplicates-and-low-quality-alignments) (including non-uniquely mapped reads)
- Check the expected mapping to exons, introns etc.

#### Remove mitochondrial reads

Remove mitochondrial reads. To assess the total % of mitochondrial reads, in an unsorted bam file, run:

```
samtools view <sample>.bam | grep chrM | wc -l 
```

To see the total number of DNA fragments, run:

```
samtools flagstat <sample>.bam > <sample>.flagstat

cat <sample>_sorted.flagstat
```

The first line shows the total number of DNA fragments. The % of DNA fragments aligned to chrM can be calculated as a % of the total DNA fragments. To remove any mitocondrial DNA, run the following:


```
samtools view -h <sample>-sorted.bam | grep -v chrM | samtools sort -O bam -o <sample>.rmChrM.bam -T .
```

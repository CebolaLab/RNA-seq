# RNA-seq

A step-by-step analysis pipeline for RNA-seq data from the [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/).

Correspondence: hannah.maude12@imperial.ac.uk

## ** \*UNDER CONSTRUCTION\* **

The resources and references used to build this tutorial are found at the bottom, in the [resources](#resources) section.

## Table of Contents

*Run using command line tools (`bash`)*:
- [Pre-alignment quality control (QC)](#pre-alignment-qc)
- [Align to the reference human genome](#align-to-the-reference-genome)
- [Post-alignment QC](#post-alignment-qc)
- [Quantify transcripts](#quantification)
- [Visualise tracks against the reference genome](#visualisation)

*Run in `R`*:
- [Differential gene expression (DGE) analysis](#differential-expression)


## Introduction

This pipeline is compatabile with RNA-seq reads generated by Illumina.

## Pre-alignment QC

#### Generate QC report

The raw sequence data should first be assessed for quality. FastQC reports can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. For paired-end reads, run fastqc on both files, with the results output to the current directory:

```bash
fastqc <sample>_1.fastq.gz -d . -o .

fastqc <sample>_2.fastq.gz -d . -o .
```

These fastQC reports can be combined into one summary report using [multiQC](https://multiqc.info/). 

#### Adapter trimming 

If there is evidence of adapter contamination shown in the fastQC report (see below), adapter sequences may need to be trimmed, using a tools such as cutadapt, trimmomatic and fastp. In this pipeline, fastp is used to trim adapters. For **paired-end** data:

```bash
fastp -i <sample>_R1.fastq.gz -O <sample>_R1.trimmed.fastq.gz -I <sample>_R2.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -l 25 -j <sample>.fastp.json -h <sample>.fastp.html
```

For **single-end** reads: (note the adapter detection is not always as effective for single-end reads, so it is advisable to provide the adapter sequence, here the 'Illumina TruSeq Adapter Read 1'):

```bash
fastp -i <sample>.fastq.gz -o <sample>-trimmed.fastq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -l 25 -j <sample>.fastp.json -h <sample>.fastp.html 
```

A html report is generated, including the following information:

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/fastp-summary.png" width="700">


## Align to the reference genome

The raw RNA-seq data in `fastq` format will be aligned to the reference genome, along with a reference transcriptome, to output two alignment files: the genome alignment and the transcriptome alignemnt. 

The DNA reads are aligned using the splice-aware aligner, STAR. Here, [STAR](https://github.com/alexdobin/STAR) is used. The manual is available [here](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf). The reference genome used is the GRCh38 'no-alt' assembly from ncbi, recommended by [Heng Li](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use). The genome can be downloaded using `wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`.  This version of the recent GRCh38 reference genome excludes alternative contigs which may cause fragments to map in multiple locations. The downloaded genome should be indexed with STAR. 

> Index the reference genome

Set --sjdbOverhang to your maximum read length -1. The indexing also requires a file containing gene annotation, which comes in a `gtf` format. For example, ENCODE provides a gtf file with GRCh38 annotations, containing gencode gene coordinates, along with UCSC tRNAs and a PhiX spike-in. Here, we use `gencode.v35.annotation.gtf` as the most recent gene annotation file. The user should aim to use the most up-to-date reference files, while ensuring that the format is the same as the reference genome. For example, UCSC uses the 'chr1, chr2, chr3' naming convention, while ENSEMBL uses '1, 2, 3' etc. The files suggested here are compatible. 

```bash
GENOMEDIR=/path/to/indexed/genome

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $GENOMEDIR --genomeFastaFiles $GENOMEDIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sjdbGTFfile gencode.v35.annotation.gtf --sjdbOverhang readlength -1
```

> Carry out the alignment

STAR can then be run to align the `fastq` raw data to the genome. If the fastq files are in the compressed `.gz` format, the `--readFilesCommand zcat` argument is added. The output file should be unsorted, as required for the downstream quantification step using Salmon. The following options are shown according to the ENCODE recommendations. 


For **paired-end** data:

```bash
STAR --runThreadN 4 --genomeDir $GENOMEREF --readFilesIn <sample>_R1.trimmed.fastq.gz <sample>_R2.trimmed.fastq.gz
--outFileNamePrefix <sample> --readFilesCommand zcat --outSAMtype BAM Unsorted --quantTranscriptomeBan Singleend --outFilterType BySJout 
--alignSJoverhangMin 8 --outFilterMultimapNmax 20
--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999
--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 
--alignIntronMax 1000000 --alignMatesGapMax 1000000 
--quantMode TranscriptomeSAM --outSAMattributes NH HI AS NM MD
```

For **single-end** data:

```bash
STAR --runThreadN 4 --genomeDir $GENOMEREF --readFilesIn <sample>-trimmed.fastq.gz 
--outFileNamePrefix <sample> --readFilesCommand zcat --outSAMtype BAM Unsorted --quantTranscriptomeBan Singleend --outFilterType BySJout 
--alignSJoverhangMin 8 --outFilterMultimapNmax 20
--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999
--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 
--alignIntronMax 1000000 --alignMatesGapMax 1000000 
--quantMode TranscriptomeSAM --outSAMattributes NH HI AS NM MD
```

*Hint: all the above code should be on one line!*

For compatibility with the STAR quantification, the `--quantMode TranscriptomeSAM` option will result in the output of two alignment files, one to the reference genome (`Aligned.*.sam/bam`) and one to the transcriptome (`Aligned.toTranscriptome.out.bam`).

#### Merge files [optional]

At this stage, if samples have been sequenced across multiple lanes, the sample files can be combined using `samtools merge`. Various QC tools can be used to assess reproducibility and assess lane effects, such as `deeptools plotCorrelation`. The `salmon` quantification does not require files to be merged, since multiple `bam` files can be listed in the command. However, to visualise the RNA-seq data from the combined technical replicates, `bam` files can be merged at this stage. For example, if your sample was split across lanes 1, 2 and 3 (`L001`, `L002`, `L003`):

```
samtools merge <sample>-merged.bam <sample>_L001.bam <sample>_L002.bam <sample>_L003.bam
```

## Post-alignment QC

First, QC reports will be generated using [qualimap](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html). Run on either the `<sample>.gzAligned.out.bam` or `<sample>.merged.bam`.

The transcript annotation file, `gencode.v35.annotation.gtf` can be downloaded from [gencode](https://www.gencodegenes.org/human/).

```bash
#Sort the output bam file. The suffix of the .bam input file may be .gzAligned.out.bam, or -merged.bam. Edit this code to include the appropriate file name.
samtools sort <sample>.bam > <sample>-sorted.bam

qualimap bamqc -bam <sample>-sorted.bam -gtf gencode.v35.annotation.gtf -outdir <sample>-bamqc-qualimap-report --java-mem-size=16G

qualimap rnaseq -bam <sample>-sorted.bam -gff gencode.v35.annotation.gtf -outdir <sample>-rnaseq-qualimap-reports --java-mem-size=16G
```

Qualimap can then run QC on combined samples. This includes principal component analysis (PCA) to confirm whether technical and/or biological replicates cluster together. A text file (`samples.txt`) should be created with three columns, the first with the sample ID, the second with the full path to the bamqc results and the third with the group names. 

*Note, some versions of qualimap require the raw_data_qualimapReport directory to be renamed to raw_data.*

```bash
qualimap multi-bamqc sample.txt
```

The QC reports can be combined using [multiqc](https://multiqc.info/); an excellent tool for combining QC reports of multiple samples into one. Example outputs of qualimap/multiqc include the alignment positions 

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/multiqc-alignment.png" width="800">

### Remove duplicates?

It is generally recommended to *not* remove duplicates when working with RNA-seq data, unless using UMIs (unique molecular identifiers) [(Klepikova et al. 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5357343/). This is because there are likely to be DNA molecules which are natural duplicates of each other, for example originating from genes with a shared sequence in a common domain. Typically, removing duplicates does more harm than good. It is more or less impossible to remove duplicates from single-end data and research has also suggested it may cause false negatives when applied to paired end data. See more in [this useful blog post](https://dnatech.genomecenter.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data/). Generally, duplicates are not a problem so long as the *library complexity is high*. 

### Check correlation of technical and biological replicates

The correlation between `bam` files of biological and/or technical replicates can be calculated as a QC step to ensure that the expected replicates positively correlate. Deeptools [multiBamSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html) and [plotCorrelation](https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html) are useful tools for further investigation.

## Quantification

The `bam` file previously aligned to the *transcriptome* by STAR will next be input into [Salmon](https://combine-lab.github.io/salmon/) in alignment-mode, in order to generate a matrix of gene counts. The Salmon documentation is available [here](https://salmon.readthedocs.io/en/latest/).

#### Generate transcriptome

Salmon requires a transcriptome to be generated from the genome `fasta` and annotation `gtf` files used earlier with STAR. This can be generated using `gffread` (source package avaiable for download [here](http://ccb.jhu.edu/software/stringtie/gff.shtml)).

```bash 
gffread -w GRCh38_no_alt_analysis_set_gencode.v35.transcripts.fa -g GCA_000001405.15_GRCh38_no_alt_analysis_set.fna gencode.v35.annotation.gtf
```

### Run Salmon

Salmon is here used with the expectation minimisation (EM)/VBEM? approach method for quantification. This is described in the 2020 paper by [Deschamps-Francoeur et al.](https://www.sciencedirect.com/science/article/pii/S2001037020303032), which describes the handling of multi-mapped reads in RNA-seq data. Duplicated sequences such as pseudogenes can cause reads to align to multiple positions in the genome. Where transcripts have exons which are similar to other genomic sequences, the EM approach attributes reads to the most likely transcript. Technical replicates can also be combined by providing the Salmon `-a` argument with a list of bam files, with the file names separated by a space (this may not work on all queue systems. A common error is `segmentation fault (core dump)`). Here, Salmon is run ***without*** any normalisation, on each technical replicate; samples are combined and normalised in the next steps.

For **paired-end** data:

```bash
salmon quant --useEM -t GRCh38_no_alt_analysis_set_gencode.v35.transcripts.fa --libType A -a <sample>.Aligned.toTranscriptome.out.bam -o <sample>.salmon_quant 
```

For **single-end** data:

If using single end data, add the `--fldMean` and `--fldSD` parameters to include the mean and standard deviation of the fragment lengths. If listing multiple files to be combined, the library type will need to be specified, as Salmon cannot determine it automatically (see the Salmon documentation for more information).

```bash
salmon quant --useEM -t GRCh38_no_alt_analysis_set_gencode.v35.transcripts.fa --libType ?? -a <sample>.Aligned.toTranscriptome.out.bam -o <sample>.salmon_quant  --fldMean ?? --fldSD ??
```

## Differential expression

***All following code should be run in `R`.***

The differential expression analysis contains the following steps:

- [Import count data](#import-count-data)
- [Filter genes](#filter-genes)
- [Normalisation](#normalisation)
- [Differential expression analysis](#differential-expression-analysis)
- [QC plots](#qc-plots)

To install the required packages:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cqn")
BiocManager::install("DESeq2")
BiocManager::install("tximport")
```

### Import count data

The output from Salmon are TPM values (the 'abundance', transcripts per million) and estimated counts mapped to transcripts. The counts will be combined to gene-level estimates in R. The output files from salmon, `quant.sf` will be imported into R using `tximport` (described in detail [here](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Downstream_DGE_in_Bioconductor) by Love, Soneson & Robinson). This will require a list of sample IDs as well as a file containing transcript to gene ID mappings, in order to convert the transcriptome alignment to gene-level counts. 

1) **Create a matrix containing the sample IDs**. The matrix should have at least three columns: the first with the sample IDs, the second with the path to the salmon `quant.sf` files, and the third with the group (e.g. treatment or sample). This can be generated in excel, for example, and saved as a tab-delimited txt file called `samples.txt`. 

```R
#Read in the files with the sample information
samples=read.table('samples.txt')
```

2) **Read in the transcript to gene ID file** provided in this repository (generated from gencode v35).

```R
#Read in the gene/transcript IDs 
tx2gene=read.table('tx2gene.txt',sep='\t')
```

3) **Read in the count data using `tximport`**. This will combine the transcript-level counts to gene-level. 

```R
library(tximport)

#Column 2 of samples, samples[,2], contains the paths to the quant.sf files
counts.imported=tximport(files=as.character(samples[,2]),type='salmon',tx2gene=tx2gene)
```

### Normalisation 

The count data needs to be normalised for several confounding factors. The number of DNA reads (or fragments for paired end data) mapped to a gene is influeced by (1) its gc-content, (2) its length and (3) the total library size for the sample. There are multiple methods used for normalisation. Here, conditional quantile normalisation (cqn) is used as recommended by [Mandelboum et al. (2019)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000481) to correct for sample-specific biases. Cqn is described by [Hansen et al. (2012)](https://academic.oup.com/biostatistics/article/13/2/204/1746212).

`cqn` requires an input of gene length, gc content and the estimated library size per sample (which it will estimate as the total sum of the counts if not provided by the user). For more guidance on how to normalise using `cqn` and import into `DESeq2`, the user is directed to [the cqn vignette](https://bioconductor.org/packages/release/bioc/vignettes/cqn/inst/doc/cqn.pdf) by Hansen & Wu and the [tximport vignette](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html) by Love, Soneson & Robinson.

```R
#Read in the gene lengths and gc-content data frame (provided in this repository)
genes.length.gc=read.table('gencode-v35-gene-length-gc.txt',sep='\t')
```

At this stage, technical replicates can be combined if they have not been already. This is typically achieved by summing the counts. 

To carry out the normalisation:

```R
library(cqn)
#cqn normalisation
counts=counts.imported$counts

#Exclude genes with no length information, for compatibility with cqn.
counts=counts[-which(is.na(genes.length.gc[rownames(counts),]$length)),]

#Extract the lengths and GC contents for genes in the same order as the counts data-frame
geneslengths=genes.length.gc[rownames(counts),]$length
genesgc=genes.length.gc[rownames(counts),]$gc

#Run the cqn normalisation 
cqn.results<-cqn(counts, genesgc, geneslengths, lengthMethod = c("smooth"))

#Extract the offset, which will be input directly into DEseq2 to normalise the counts. 
cqnoffset <- cqn.results.DEseq$glm.offset
cqnNormFactors <- exp(cqnoffset)

#The 'counts' object imported from tximport also contains data-frames for 'length' and 'abundance'.
#These data-frames should also be subset to remove any genes excluded from the 'NA' length filter
counts.imported$abundance=counts.imported$abundance[rownames(counts),]
counts.imported$counts=counts.imported$counts[rownames(counts),]
counts.imported$length=counts.imported$length[rownames(counts),]
```

The normalised gene expression values can be saved as a cqn output. These values will not be used for the downstream differential expression, rather they are useful for any visualisation purposes. Differential expression will be calculated within DEseq2 using a negative bionomial model, to which the cqn offset will be added. 

```R
#The normalised gene expression counts can be saved as:
RPKM.cqn<-cqn.results$y + cqn.results$offset
```

### Import data to DEseq2 

The counts information will be input into DEseq2. A data-frame called `colData` should be generated. The rownames will be the unique sample IDs, while the columns should contain the conditions being tested for differential expression, in addition to any batch effects. In the example below, the column called `condition` contains the treatment, while the column `batch` contains the donor ID. 

```R
#Import to DEseq2
counts.DEseq=DESeqDataSetFromTximport(counts.imported,colData=colData,design=~condition+batch)

dds <- DESeq(counts.DEseq)
resultsNames(dds) #lists the coefficients

#Add the normalisation offset from cqn
normalizationFactors(dds) <- cqnNormFactors
```

### Sample clustering (PCA)

A common component of analysing RNA-seq data is to carry out QC by testing if expected samples cluster together. One popular tool is principal component analysis (PCA) (the following steps are adapted from the [hbctraining tutorial on clustering](https://github.com/hbctraining/DGE_workshop_salmon/blob/master/lessons/03_DGE_QC_analysis.md)).

**If you have few samples:**

```R
rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
```

**If you have more samples** e.g. >20:

```R
vst.r<-vst(dds,blind=TRUE)
vst_mat <- assay(vst)
pca <- prcomp(t(vst_mat))
```

Plot the results

```R
library(ggplot2)

df_out <- as.data.frame(pca$x)
df_out$group <- samples[,3]

#Include the next two lines to add the PC % to the axis labels
percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), paste0(" (", as.character(percentage), "%", ")"), sep="") 

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()+theme + xlab(percentage[1]) + ylab(percentage[2])

#png()
print(p)
#dev.off()
```


sample-level QC using Principal Component Analysis (PCA)

A useful tutorial on PCA and heirarchical clustering is available [here](https://github.com/hbctraining/DGE_workshop_salmon/blob/master/lessons/03_DGE_QC_analysis.md).

If you have few samples, it is recommended to use `rld` transformation.
If you have >20 samples, it may be faster to use `vst`

### Differential gene expression


There are several models available to calculate differential gene expression. Here, the apeglm shrinkage method will be applied to shrink high log-fold changes with little statistical evidence and account for lowly expressed genes with significant deviation.

```R
library(apeglm)
ins.LFC <- lfcShrink(dds, coef="condition_insulin_vs_control", type="apeglm")
```


### Filter genes 

There are likely to be many annotated genes which are not expressed in your samples. To avoid these influencing the downstream results, they should be removed at this stage. In order to account for the varying library size, it is recommended to filter genes based on their 'count per million' (CPM) expression. (The count number is divided by the total count number for the sample, divided by one million). The edgeR package includes a handy tool, `filterByExpr`, which carrys out informative gene filtering. 

```R
#The groups can be updated to reflect the biological replicates
groups=unique(samples[,3:4])

#The genes are filtered to remove those with low expression
library(edgeR)

counts.combined=counts.combined[filterByExpr(counts.combined,design=groups[,2]),]
```


### Differential expression analysis

```R
offset=cqn.results$glm.offset

#Make the edgeR DGEList object and input the cqn offset
y <- DGEList(counts=counts.combined,group=groups[,2])
y$offset <- offset

#If you wish all the comparisons to be relative to your first group (e.g. a control), remove the 0+
design <- model.matrix(~ 0+group,data=group)
y <- estimateGLMCommonDisp(y, design = design)

#Run the quasi-likelihood, glm fit
QLfit <- glmQLFit(y,design)
```

To carry out the differential expression comparison, the groups to be compared should first be defined. Here, it is assumed that each group contained two or more biological replicates. 

```R
#Replace 'group1' and 'group 2' with the names of your groups which you wish to compare
contrast1=makeContrasts(group2-group1)

ql.groups12=glmQLFTest(QLfit, contrast=contrast1, levels=design)

```

### QC plots

Before moving on to functional analysis, such as gene set enrichment analysis, quality control should be carried out on the differential expression analyses. The types of plots which will be generated below are:

- [Biological replicate correlation](#biological-replicate-correlation)
- [MD plot](#md-plot)
- [p-value distribution](#p-value-distribution)
- [Volcano plot](#volcano-plot)


#### Biological replicate correlation

The correlation between the expression of genes in two biological replicates should ideally be very high. The normalised expression values, saved above as `RPKM.cqn` can be used 

```R
#To test the correlation between the first two samples in columns 1 and 2
plot(RPKM.cqn[,1],RPKM.cqn[,2],pch=18,cex=0.5,xlab=colnames(RPKM.cqn)[1],ylab=colnames(RPKM.cqn)[2])

#The Pearson correlation coefficient can be calculated as:
cor(RPKM.cqn[,1],RPKM.cqn[,2])

#To add it to your plot, replace x and y with the coordinates for your legend
text(x,y,labels=paste0('r=',round(cor(RPKM.cqn[,1],RPKM.cqn[,2]),2)))
```

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/biological-replicates-exp-logcqn.png" width="800">

#### MD plot

A mean-difference (MD) plot plots the log fold-change between two samples against the average gene expression. An MD plot of log2 fold-change against average expression (CPM, counts per million) can be plotted using the following:

```R
plotMD(ql.groups12,main='Group1 vs Group2',cex=0.5)
```

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/FC-CPM.png" width="800">

#### Distribution of p-values and FDRs

The distribution of p-values following a differential expression analysis can be an indication of whether the experiment worked. 

```R
#The distribution of p-values
hist(ql.groups12$table$PValue,col='grey',breaks=50,main='p108_T vs p108',xlab='p-values')

#The false-discovery rate distribution
hist(topTags(ql.groups12,n=length(genes))$table$FDR,col='grey',breaks=50,main='p108_T vs p108',xlab='FDR')
```

The p-value distribution:
<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/P-values-hist.png" width="800">

The false discovery rate (FDR) distribution:
<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/FDR-hist.png" width="800">

#### Volcano plots

A volcano plot is a scatterplot which plots the p-value of differential expression against the fold-change. The volcano plot can be designed to highlight datapoints of significant genes, with a p-value and fold-change cut off.

Volcano plots are generated as described by [Ignacio González](http://www.nathalievialaneix.eu/doc/pdf/tutorial-rnaseq.pdf)

```R
#Allow for more space around the borders of the plot
par(mar = c(5, 4, 4, 4))
#Set your log-fold-change and p-value thresholds 
lfc = 2
pval = 0.05

tab = data.frame(logFC = ql.groups12$table[, 1], negLogPval = -log10(ql.groups12$table[, 4]))

plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue),main='p53_T vs p53')

#Genes with a fold-change greater than 2 and p-value<0.05:
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))

#Colour these red
points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red")

#Show the cut-off lines
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)

mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.6, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.6, line = 0.5)

```

The resulting plot will look like this:

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/volcano-plot-QL-0.05.png" width="600">

## Functional analysis

Functional analysis can further investigate the differential expression of each gene. Pathway analysis is a popular approach with which to investigate the differential expression of pathways, including genes with similar biological functions. This can be achieved using gene set analysis (GSA). There are many flavours of GSA. They can be categorised as shown below by [Das et al. (2020)](https://www.mdpi.com/1099-4300/22/4/427).

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/GSA-classification.png" width="600">

Here, we will use one gene annotation approach and one gene set enrichment analysis (GSEA) approach. 

### GoSeq - gene annotation

[GoSeq](https://bioconductor.org/packages/release/bioc/html/goseq.html), developed by [Young et al. (2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-2-r14), tests for the enrichment of Gene Ontology terms.

```R
BiocManager::install("goseq")
library(goseq)

#Extract the differential expression data, with false discovery rate correction
groups12.table<-as.data.frame(topTags(ql.groups12,n= Inf))

#Remove the version numbers from the ENSEMBL gene IDs
rownames(groups12.table)<-unlist(lapply(strsplit(rownames(groups12.table),'\\.'),`[[`,1))
```

The genes can be seperated into those which show significantly increased expression and those which show significantly decreased expression. Here, the FDR threshold is set to 0.05. A minimum fold-change can also be defined. 

```R
#Decreased expression
ql53.DEGs.down <- groups12.table$FDR < 0.05 & groups12.table$logFC<0
names(ql53.DEGs.down) <- rownames(groups12.table)
pwf.dn <- nullp(ql53.DEGs.up, "hg19","ensGene")
go.results.dn <- goseq(pwf.dn, "hg19","ensGene")

#Increased expression
ql53.DEGs.up <- groups12.table$FDR < 0.05 & groups12.table$logFC>0
names(ql53.DEGs.up) <- rownames(groups12.table)
pwf.up <- nullp(ql53.DEGs.down, "hg19","ensGene")
go.results.up <- goseq(pwf.up, "hg19","ensGene")
```

The `go.results.up` dataframe looks like this:

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/Go-output.png" width="800">

Significant results can be saved...

```R
write.table(go.results.up[go.results.up$over_represented_pvalue<0.05,1:2],'p53-GO-up0.05.txt',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
```

...and uploaded to the REVIGO tool which collapses and summarises redundant GO terms. Copy the contents of the `p53-GO-up0.05.txt` into the [REVIGO](http://revigo.irb.hr/) box:

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/REVIGO-input.png" width="500">

After running, select the 'Scatterplot & Table' tab and scroll down to 'export results to text table (csv)':

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/REVIGO2.png" width="800">

Save the downloaded file as `REVIGO-UP.csv`. The package `ggplot2` can be used here to visualise the log<10> p-values for GO term enrichment. To view the top 10 terms:

```R
#Read in the REVIGO output
revigoUP=read.table('REVIGO-UP.csv',sep=',',header=TRUE)

#Sort by p-value and extract the top 20
revigoUP=revigoUP[order(revigoUP$log10.p.value),]
revigoUP=head(revigoUP,n=20)

#Convert the GO terms to factors for compatability with ggplot2
revigoUP$description <- factor(revigoUP.108top$description, levels = revigoUP.108top$description)

#Plot the barplot
p<-ggplot(data=revigoUP.108top, aes(x=log10.p.value, y=description, fill=description)) +
  geom_bar(stat="identity") 

p + scale_fill_manual(values=rep("steelblue2",dim(revigoUP.108top)[1])) + theme_minimal() + theme(legend.position = "none") + 
        ylab('')
```

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/REVIGO-UP.png" width="800">

### GSEA 

```R
BiocManager::install("piano")
library(piano)
```


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




$salmon quant  -t $GENOMEDIR/gencode.v35.transcripts.fa --libType A -a IGF0005790/p53-rep2-L002.gzAligned.toTranscriptome.out.bam -o "$sampleID".salmon_quant --fldMean $mean --fldSD 75 --seqBias --gcBias 




### Compute GC bias

GC-bias describes the bias in sequencing depth depending on the GC-content of the DNA sequence. Bias in DNA fragments, due to the GC-content and start-and-end sequences, may be increased due to preferential PCR amplification [(Benjamini and Speed, 2012)](https://academic.oup.com/nar/article/40/10/e72/2411059). A high rate of PCR duplications, for example when library complexity is low, may cause a significant GC-bias due to the preferential amplification of specific DNA fragments. This can significantly impact transcript abundance estimates. Bias in RNA-seq is explained in a handy [blog](https://mikelove.wordpress.com/2016/09/26/rna-seq-fragment-sequence-bias/) and [video](https://youtu.be/9xskajkNJwg) by Mike Love.

It is **crucial** to correct GC-bias when comparing groups of samples which may have variable GC content dependence, for example when samples were processed in different libraries. `Salmon`, used later to generate read counts for quantification, has its own in-built method to correct for GC-bias. When generating `bedGraph/BigWig` files for visualisation, the user may opt to correct GC-bias so that coverage is corrected and appears more uniform. The `deeptools` suite includes tools to calculate GC bias and correct for it.

The reference genome file should be converted to `.2bit` format using [`faToTwoBit`](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit).
The effective genome size can be calculated using `faCount` available [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/).
Set the `-l` argument to your fragment length.

The input `bam` file requires an index, which can be generated using `samtools index`.

```bash
deeptools computeGCBias -b <sample>-sorted.bam --effectiveGenomeSize 3099922541 -g GCA_000001405.15_GRCh38_no_alt_analysis_set.2bit -l 100 --GCbiasFrequenciesFile <sample>.freq.txt  --biasPlot <sample>.biasPlot.pdf
```

The bias plot format can be changed to png, eps, plotly or svg. If there is significant evidence of a GC bias, this can be corrected using `correctGCbias`. An example of GC bias can be seen in the plot outout from `computeGCBias` below:

<img src="https://github.com/CebolaLab/RNA-seq/blob/master/Figures/GCbiasPlot.png" width="500">

The GC-bias will later be corrected when normalising the data for DGE analysis. However, for visualisation purposes, the user may wish to correct the GC-bias for the corresponding bigWig file.

If opting to correct the GC-bias, `correctGCbias` can be used. This tool effectively removes reads from regions with greater-than-expected coverage (GC-rich regions) and removes reads from regions with less-than-expected coverage (AT-rich regions). The methods are described by [Benjamini and Speed [2012]](https://academic.oup.com/nar/article/40/10/e72/2411059). The following code can be used:

```bash
 correctGCBias -b <sample>-sorted.bam --effectiveGenomeSize 3099922541 -g GCA_000001405.15_GRCh38_no_alt_analysis_set.2bit --GCbiasFrequenciesFile <sample>.freq.txt -o <sample>.gc_corrected.bam [options]
```

**NOTE:** When calculating the GC-bias for ChIP-seq, ATAC-seq, DNase-seq (and CUT&Tag/CUT&Run) it is recommended to filter out problematic regions. These include those with low mappability and high numbers of repeats. The compiled list of [ENCODE blacklist regions](https://www.nature.com/articles/s41598-019-45839-z) should be excluded. However, the ENCODE blacklist regions have little overlap with coding regions and this step is not necessary for RNAs-seq data [(Amemiya et al, 2019)](https://www.nature.com/articles/s41598-019-45839-z).


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


## Resources

Many resources were used in building this RNA-seq tutorial.

Highly recommended RNA-seq tutorial series:

- [Introduction to differential gene expression analysis](https://hbctraining.github.io/DGE_workshop/lessons/01_DGE_setup_and_overview.html#rna-seq-count-distribution)

Other RNA-seq tutorials:
- [Statistical analysis of RNA-Seq data](http://www.nathalievialaneix.eu/doc/pdf/tutorial-rnaseq.pdf) by Ignacio Gonz´alez 
- [Analysis of RNA-Seq data: gene-level exploratory analysis and differential expression](https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#gene-ontology-enrichment-analysis) by Bernd Klaus and Wolfgang Huber


- <https://vallierlab.wixsite.com/pipelines/rna-seq>
- [RNA-seq workflow: gene-level exploratory analysis and differential expression by Love et al. 2019](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#running-the-differential-expression-pipeline)
- https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
- The Encode pipeline for long-RNAs: https://www.encodeproject.org/data-standards/rna-seq/long-rnas/

Understanding normalisation:
- https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
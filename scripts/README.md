# scripts, softwares and R packages
---
## ATAC-seq
Scripts used for mapping raw ATAC-seq reads and peak calling

## auxillary_scripts
Scripts with function definitions used throughout the project

## histone_marks
- Scripts to compare histone modification levels between slow and fast DMRs
- Scripts to fit a linear regression model to predict methylation recovery using all the histone marks

## nematostella_genome
Scripts to create files related to Nematostella genome annotation, including index, repeats and transcripts.

## nematostella_methylation
Scripts to process raw reads, identify DMRs, distinguish methylated and unmethylated genes and perform per-read methylation analyses

## populations_methylation
Scripts to identify DMRs in different Nematostella populations compared to the our lab strain samples

## rnaseq
Scripts to quantify RNA-seq counts

---

To perform the analysis in this manuscript you can download the original software from these sources:

## Nanopore basecalling:

**nanopore_guppy** —> https://nanoporetech.com/es/document/Guppy-protocol

**dorado** —> https://github.com/nanoporetech/dorado

**samtools** —> https://github.com/samtools/samtools

**modbam2bed** —> https://github.com/epi2me-labs/modbam2bed

**modkit** —> https://github.com/nanoporetech/modkit

## Methylation analysis:

**bedGraphToBigWig** —> https://www.encodeproject.org/software/bedgraphtobigwig/

**bedtools** —> https://github.com/arq5x/bedtools2

**Deeptools2** —> https://deeptools.readthedocs.io/en/latest/

**methylartist** —>https://github.com/adamewing/methylartist 

## ATAC-seq analysis:

**chromap** —> https://github.com/haowenz/chromap 

**MACS3** —> https://github.com/macs3-project/MACS 

## RNA-seq analysis:

**hisat2** —> https://github.com/DaehwanKimLab/hisat2

**kallisto** —> https://github.com/pachterlab/kallisto?tab=readme-ov-file

**stringtie** —> https://github.com/gpertea/stringtie

**TElocal** —> https://github.com/mhammell-laboratory/TElocal 

## Sequence analysis:

**seqkit** —> https://github.com/annalam/seqkit

## Repeats annotation:

**RepeatModeler** —> https://github.com/Dfam-consortium/RepeatModeler

**RepeatMasker** —> https://github.com/Dfam-consortium/RepeatMasker

## Haplotype phasing:

**Clair3** —> https://github.com/HKU-BAL/Clair3 

**WhatsHap** —> https://github.com/whatshap/whatshap 

## R packages:
Biostrings=2.77.2

BSgenome=1.77.2

ComplexHeatmap=2.26.0

DESeq2=1.49.4

doParallel=1.0.17

dplyr=1.1.4

DSS=2.58.0

GenomicFeatures=1.61.6

GenomicRanges=1.61.5

methodical=1.7.0

methrix=1.7.0

pheatmap=1.0.13

regioneR=1.41.3

VariantAnnotation=1.56.0

Rsamtools=2.25.3

rtracklayer=1.69.1

relaimpo=2.2.7

HDF5Array=1.37.0

bsseq=1.42.0

topGO=2.58.0

CAGEr=2.14.0

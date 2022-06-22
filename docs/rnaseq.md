# RNA-seq

## Overview

I used the site [HBC Github Salmon tutorial](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/04_quasi_alignment_salmon.html) as a guide for my analysis.

| Sample       | Accession | Link                                                         |
|--------------|-----------|--------------------------------------------------------------|
| ILC3 RNA-seq | GSE130775 | [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130775) |
| CD4 RNA-seq  | GSE87254  | [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87254)  |

### Sample Overveiw

There are multiple replicates for each type of cell line.

| Accession  | CellType      |
|------------|---------------|
| SRX5797708 | ILC3b         |
| SRX5797709 | ILC3a         |
| SRX5797712 | ILC3b         |
| SRX5797713 | ILC3a         |
| SRX5797716 | ILC3b         |
| SRX5797717 | ILC3a         |
| SRR4290846 | CD4_resting   |
| SRR4290852 | CD4_resting   |
| SRR4290853 | CD4_activated |
| SRR4290854 | CD4_activated |

## Mapping with Salmon

I used Salmon to quasi-map reads to transcripts. The data was aligned to hg38 annotations from ENSEMBL. I used the following site to fetch the hg38 genome [hg38 ENSEMBL](http://uswest.ensembl.org/Homo_sapiens/Info/Index).

Example command:

```bash
salmon quant -i ./salmon_index -l A \
--validateMappings \
-p 4 \
--seqBias \
--gcBias \
--posBias \
-r ${file} \
-o ./RNAseq/`basename ${file} .fastq`
```

## Gene Counts

I imported the transcript counts and collapsed them to gene counts using [Tx import](https://bioconductor.statistik.tu-dortmund.de/packages/3.3/bioc/vignettes/tximport/inst/doc/tximport.html). The code for processing Salmon mapped reads can be found in [`hILC_DESEq2.r`](../r/hILC_DESEq2.r).

I exported the gene expression matrix for ILC3 and CD4+ populations to the directory [`./data/counts/`](../data/RNA/counts/). The combined counts for both cell types can be found in the file [`20200706_hILC3_CD4_GeneCounts.tsv`](../data/RNA/counts/20200706_hILC3_CD4_GeneCounts.tsv).

The raw outputs for the RNA-seq data analysis are stored on [Google Drive](https://drive.google.com/drive/folders/14I2-kREIk8YKnGRQfkh6vHXpb0kXi4yJ?usp=sharing).
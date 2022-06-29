# ChIP-seq

1) I obtained ILC3 histone marker data from [GSE77299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) and CD4+ histone marker data from [Blueprint Project](http://dcc.blueprint-epigenome.eu/#/experiments).

## Workflow Overview 

### Read QC, Alignment, and Filtering

Using snakemake, we followed ENCODE3 standards for ChIP-seq read alignment, read filtering, and peak calling. In brief, fastq files were assessed for adapter contamination and read quality statistics (FastQC v. 0.11.9). Samples flagged for high levels of N sequences were removed. TrimGalore! (v. 0.6.7) was used to remove adapter contamination and trim the low-quality bases at the 3’ end of the sequencing read with the settings (-q 20). 

### Read Alignment and Filtering

Reads were aligned to the hg38 reference genome using bowtie2 (v. 2.4.4) (--very-sensitive --maxins 2000). The aligned reads were quality filtered with samtools (v. 1.9) (-F 1804 -q 30) and PCR duplicates were removed with samtools markdup (-r -s). 

### Peak Calling

MACS2 (v 2.2.2.7.1) was used to call peaks with parameters (--nomodel --extsize 150). We provided all filtered BAM files during peak calling for cell types and TF combinations with multiple biological replicates.
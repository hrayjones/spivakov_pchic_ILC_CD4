# ATAC-seq

Data were downloaded and processed to Tn5 cut sites using the snakemake pipeline used for [maxATAC](https://github.com/tacazares/snakeATAC).

| Sample    | Accession   | Publication Link                                                         |
|-----------|-------------|--------------------------------------------------------------|
| ILC3 ATAC | GSE77299    | [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130775) |
| CD4 ATAC  | PRJNA380283 | [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5623106)  |

## Workflow Summary

### Read Filtering

We evaluated sequencing quality with FastQC. Sequencing adapters and bases with a PHRED score < 30 were trimmed with the package Trim Galore! using the parameters (-q 30 -paired). We used bowtie2 (v. 2.4.4) with parameters (-p 8 --very-sensitive --maxins 2000) to align reads to the hg38 reference genome. We removed reads with samtools that had MAPQ score of < 30 and samflag 3. For all alignments, we removed duplicates (i.e., potential PCR artifacts) with samtools rmdup and samtools fixmate with parameter (-n). We remove reads mapping to the mitochondrial chromosome.

### Generate Tn5 Cut Sites

 Reads were shifted +4 on the (+) strand or -5 on the (-) strand and resized to 40bp wide intervals so that the intervals are centered at the Tn5 cut site. We called peaks with MACS2 with parameter settings (-f BED -shift=0 -ext=40 -keep-dup=all) to center the signal over the Tn5 insertion, smooth by extension +/-20bp and ensure that each inferred Tn5 binding site contributes to the peak call.

### Replicate Correlation and Quality

 DeepTools multiBigWigSummary was used to find the correlation across 1,000 bp bins along the genome. The replicate correlation between the ATAC-seq samples was poor; for tonsil ILC3, the overlap between biological replicates was < 10%. Due to the high level of sequence duplication, low peak count (8,852), and low overlap with biological replicate the results are limited to the ILC3 sample with best quality metrics (SRR3129113).
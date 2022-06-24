# pcHi-C Integrated with RNA-seq

This document will walkthrough how I used gene expression data with pcHi-C data. This file uses the output from the [pchic](./pchic.md), [ChIP-seq](./chipseq.md), and [ATAC-seq](atacseq.md) walkthroughs. The companion notebook for this analysis is [20220606_spivakov_pchic_reanalysis.ipynb](../notebooks/20220606_spivakov_pchic_reanalysis.ipynb).

## Workflow Overview

This analysis uses bash and python to perform the analysis. 

### Intersect PIR with genetic features

I used `bedtools intersect -c` to count the number of features that overlapped PIR intervals.

```bash
for PIR_BED in /Users/caz3so/workspaces/tacazares/pchic/data/CHICAGO/hg38/PIR/CD4*bed.gz;
    do
        bedtools intersect -a ${PIR_BED} -b ../data/ATAC/CD4_ATAC_peaks.bed -c > ../data/PIR_overlap/`basename ${PIR_BED} .bed.gz`_overlapATAC.bed
        bedtools intersect -a ${PIR_BED} -b ../data/CHIP/S008H1H1.ERX547940.H3K27ac.bwa.GRCh38.20150527.bed -c > ../data/PIR_overlap/`basename ${PIR_BED} .bed.gz`_overlapH3K27ac.bed
        bedtools intersect -a ${PIR_BED} -b ../data/CHIP/S008H1H1.ERX547958.H3K4me3.bwa.GRCh38.20150527.bed -c > ../data/PIR_overlap/`basename ${PIR_BED} .bed.gz`_overlapH3K4me3.bed
        bedtools intersect -a ${PIR_BED} -b ../data/RE/CD4_RE.bed -c > ../data/PIR_overlap/`basename ${PIR_BED} .bed.gz`_overlapRE.bed
    done

for PIR_BED in /Users/caz3so/workspaces/tacazares/pchic/data/CHICAGO/hg38/PIR/ILC*bed.gz;
    do
        bedtools intersect -a ${PIR_BED} -b ../data/ATAC/ILC3_ATAC_peaks.bed -c > ../data/PIR_overlap/`basename ${PIR_BED} .bed.gz`_overlapATAC.bed
        bedtools intersect -a ${PIR_BED} -b ../data/CHIP/ILC3_H3K27ac_peaks.bed -c > ../data/PIR_overlap/`basename ${PIR_BED} .bed.gz`_overlapH3K27ac.bed
        bedtools intersect -a ${PIR_BED} -b ../data/CHIP/ILC3_H3K4me3_peaks.bed -c > ../data/PIR_overlap/`basename ${PIR_BED} .bed.gz`_overlapH3K4me3.bed
        bedtools intersect -a ${PIR_BED} -b ../data/RE/ILC3_RE.bed -c > ../data/PIR_overlap/`basename ${PIR_BED} .bed.gz`_overlapRE.bed
    done
    
```

### Import the data into python as a dictionaries

Next we import the feature counts into python as a dictionary of `{feature: count}`.

```python
#ILC3
ILC3_fragments_H3K27ac_dict = utils.map_counts("./features/PIR_overlap/hILC3_10K_dpnII_fragments_PIR_overlapH3K27ac.bed")
ILC3_bins_H3K27ac_dict = utils.map_counts("./features/PIR_overlap/hILC3_10K_dpnII_5kbin_PIR_overlapH3K27ac.bed")

ILC3_fragments_H3K4me3_dict = utils.map_counts("./features/PIR_overlap/hILC3_10K_dpnII_fragments_PIR_overlapH3K4me3.bed")
ILC3_bins_H3K4me3_dict = utils.map_counts("./features/PIR_overlap/hILC3_10K_dpnII_5kbin_PIR_overlapH3K4me3.bed")

ILC3_fragments_ATAC_dict = utils.map_counts("./features/PIR_overlap/hILC3_10K_dpnII_fragments_PIR_overlapATAC.bed")
ILC3_bins_ATAC_dict = utils.map_counts("./features/PIR_overlap/hILC3_10K_dpnII_5kbin_PIR_overlapATAC.bed")

ILC3_fragments_RE_dict = utils.map_counts("./features/PIR_overlap/hILC3_10K_dpnII_fragments_PIR_overlapRE.bed")
ILC3_bins_RE_dict = utils.map_counts("./features/PIR_overlap/hILC3_10K_dpnII_5kbin_PIR_overlapRE.bed")

#CD4
CD4_50K_dpnII_fragments_H3K27ac_dict = utils.map_counts("./features/PIR_overlap/CD4_50K_dpnII_fragments_PIR_overlapH3K27ac.bed")
CD4_50K_dpnII_bins_H3K27ac_dict = utils.map_counts("./features/PIR_overlap/CD4_50K_dpnII_5kbin_PIR_overlapH3K27ac.bed")

CD4_50K_dpnII_fragments_H3K4me3_dict = utils.map_counts("./features/PIR_overlap/CD4_50K_dpnII_fragments_PIR_overlapH3K4me3.bed")
CD4_50K_dpnII_bins_H3K4me3_dict = utils.map_counts("./features/PIR_overlap/CD4_50K_dpnII_5kbin_PIR_overlapH3K4me3.bed")

CD4_50K_dpnII_fragments_ATAC_dict = utils.map_counts("./features/PIR_overlap/CD4_50K_dpnII_fragments_PIR_overlapATAC.bed")
CD4_50K_dpnII_bins_ATAC_dict = utils.map_counts("./features/PIR_overlap/CD4_50K_dpnII_5kbin_PIR_overlapATAC.bed")

CD4_50K_dpnII_fragments_RE_dict = utils.map_counts("./features/PIR_overlap/CD4_50K_dpnII_fragments_PIR_overlapRE.bed")
CD4_50K_dpnII_bins_RE_dict = utils.map_counts("./features/PIR_overlap/CD4_50K_dpnII_5kbin_PIR_overlapRE.bed")
```
# pcHi-C

___

The data were provided from two different analysis methods:

* fragments - full restriction fragments are used

* 5kb bin level - interactions are approximately 5k bins instead of the size of restriction fragments

The notebook [`20220606_spivakov_pchic_reanalysis.ipynb`](notebooks/20220606_spivakov_pchic_reanalysis.ipynb) describes how the pcHi-C and ABC model results were combined to create feature files.

___

## Filtering Criteria

The filtering used is based on conversatation between Tareian Cazares and Valerya Malysheva in June 2022.

### Trans-chromosomal Interactions

We remove trans-chromosomal interactions.
  
> Removing trans-chromosomal interactions can be achieved by removing the NA distance values.

We make sure to remove the trans-chromosomal interactions by remove any rows in our data frame that have where the bait chromosome is different from the oe chromosome.

### Self interacting regions

* We remove rows that have a distance with Na values. In my original filtering I do not remove dist = 0.

> There should be no interactions with dist = 0 in the Chicago peak matrices. The Capture Hi-C-adapted ABC analysis can generate pairs with dist = 0 and I should have removed those already from the peak matrices that I sent you, but it wouldn’t harm to remove again in your pipeline just in case.

### Filtering by CHiCAGO score

* We keep lines with a CHiCAGO score >= 5.

> CHiCAGO score >=5 is perfect, when working with peakmatrices generated from CHiCAGO interactions only. However, for the pekamatrices, containing both CHiCAGO contacts and ABC pairs, please use the merged_score column >=5 for filtering.

* For the CD4+ T cell data, there were two columns before, a column for 1 million and 50k. I used these columns to split the files into different feature files. Do we want to keep them as one experiment now?

> For the CD4+ data I have now run CHiCAGO on 4 replicates (two reps for 1M and two reps for 50K). So, for filtering you should just use either the score scolumn (in case of the CHiCAGO peakmatrix) and merged_score column (in case of the CHiCAGO_ABC peakmatrix). So there is no need to look at the scores in previous peakmatrices.

The scores column name changes based on the input analysis and date of the analysis. The funciton I wrote will look for the column name if given and filter by the specific threshold.

### Remove Line column

* Some files have a column “remove_line”. Should I use this column for filtering as opposed to filtering on my end?

> Please ignore the remove_line column.

### Off target interactions

* There used to be an off-target category in the OEName, but I have not seen them in the new files. Do we need to worry about those?

> Don’t worry about the off-targets in the OEName column.

We will remove the off target interactions for the baits and keep the off target interactions for the OEs.

### Promoter to promoter interactions

* We at one point included promoter-to-promoter interactions then excluded them.

> The promoter-promoter interactions should be removed for the RELI analysis. Have you noticed any difference for when you run with and without these interactions?

To remove promoter to promoter interactions I removed rows with OEnames. This might not be the best approach. 

I also found all unique bait intervals and removed rows where the promoter interval was found.
___

## Workflow Overview

### Convert the files to PIR `.bed` files

I wrote a python [class object](../python/ChicagoData.py) to work with the CHiCAGO output text file. This piece of code will perform filtering of specific types of interactions, like promoter-to-promoter, or trans-chromosomal interactions.

Example:

```python
import pandas as pd
import numpy as np
import ChicagoData

input_file = "/Users/caz3so/scratch/20220606_spivakov_pchic_reanalysis/TransferXL-089FGscZhgKG8/ILC_5kb_within_newbmap_CHiCAGO_ABC_peakm.txt"

ILC3_data = ChicagoData(input_file)

ILC3_data = ChicagoData(input_file, 
                        drop_off_target_bait=True, 
                        drop_off_target_oe=False, 
                        drop_trans_chrom=True,
                        score_col="merged_score",
                        score_val=5,
                        remove_p2p=True)
                        
ILC3_data.pir_df
```

### Intersect PIR `.bed` files with ChIP-seq `.bed` files

```bash
for PIR_BED in ${PIR_DIR}/CD4*bed;
    do
        bedtools intersect -a ${PIR_BED} -b /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/CHIP_ATAC/CD4/Primary_CD4_ATAC.bed -c > /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/PIR_overlap/`basename ${PIR_BED} .bed`_overlapATAC.bed
        bedtools intersect -a ${PIR_BED} -b /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/CHIP_ATAC/CD4/S008H1H1.ERX547940.H3K27ac.bwa.GRCh38.20150527.bed -c > /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/PIR_overlap/`basename ${PIR_BED} .bed`_overlapH3K27ac.bed
        bedtools intersect -a ${PIR_BED} -b /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/CHIP_ATAC/CD4/S008H1H1.ERX547958.H3K4me3.bwa.GRCh38.20150527.bed -c > /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/PIR_overlap/`basename ${PIR_BED} .bed`_overlapH3K4me3.bed
    done

for PIR_BED in ${PIR_DIR}/hILC*bed;
    do
        bedtools intersect -a ${PIR_BED} -b /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/CHIP_ATAC/ILC3/SRR3129113_end2end_final_blacklisted_IS_peaks.bed -c > /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/PIR_overlap/`basename ${PIR_BED} .bed`_overlapATAC.bed
        bedtools intersect -a ${PIR_BED} -b /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/CHIP_ATAC/ILC3/ILC3_H3K27ac_peaks.bed -c > /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/PIR_overlap/`basename ${PIR_BED} .bed`_overlapH3K27ac.bed
        bedtools intersect -a ${PIR_BED} -b /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/CHIP_ATAC/ILC3/ILC3_H3K4me3_peaks.bed -c > /Users/caz3so/scratch/20200629_Spivakov_pcHiC_analysis_summary/features/PIR_overlap/`basename ${PIR_BED} .bed`_overlapH3K4me3.bed
    done
    
```
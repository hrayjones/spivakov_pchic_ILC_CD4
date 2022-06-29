# pcHi-C Integrated with RNA-seq

This document will walkthrough how I used gene expression data with pcHi-C data. This file uses the output from the [pchic](./pchic.md), [ChIP-seq](./chipseq.md), and [ATAC-seq](atacseq.md) walkthroughs. The companion notebook for this analysis is [20220606_spivakov_pchic_reanalysis.ipynb](../notebooks/20220606_spivakov_pchic_reanalysis.ipynb).

## Workflow Overview

The code for generating the gene expression analysis is written into the `ChicagoData` python object. 

```python
# Read file into DF
self._read_file_()

# Format the DF
self._format_file_()

# Filter the formatted DF
self._filter_file_()

# Get the PIR df
self._get_PIR_df_()

# Get the bait df
self._get_bait_df_()

# Get the combined df
self._get_combined_df_()  

# Get number of features overlapping PIRs using bedtools and report the count per PIR
self._get_feature_counts_()

# Import a 2 column gene expression table that is Gene Name, Gene Expression
self._import_gene_counts_()

# Map the the different counts back to the gene expression table. Count # of PIR 
# attached to gene promoter. Sum the number of features overlapping PIRs.
self._map_feature_counts_to_genes_()

# Filter the expression matrix to remove NA 
self._filter_expression_()

# Get the count v mean gene expression matrix
self._get_PIR_count_v_mean_()

# Write the information
self._write_new_chicago_data_()
```
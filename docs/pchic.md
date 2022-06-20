# pcHi-C

The data were provided from two different analysis methods:

* fragments - full restriction fragments are used

* 5kb bin level - interactions are approximately 5k bins instead of the size of restriction fragments

The notebook [`20220606_spivakov_pchic_reanalysis.ipynb`](notebooks/20220606_spivakov_pchic_reanalysis.ipynb) describes how the pcHi-C and ABC model results were combined to create unique feature files.

## Filtering Criteria

The filtering used is based on converstation between Tareian Cazares and Valerya Malysheva in June 2022.

* We remove trans-chromosomal interactions.
  
> Removing trans-chromosomal interactions can be achieved by removing the NA distance values. 

* We remove distance with Na values. In my original filtering I do not remove dist = 0.

> There should be no interactions with dist = 0 in the Chicago peak matrices. The Capture Hi-C-adapted ABC analysis can generate pairs with dist = 0 and I should have removed those already from the peak matrices that I sent you, but it wouldn’t harm to remove again in your pipeline just in case.

* We keep lines with a CHiCAGO score >= 5.

> CHiCAGO score >=5 is perfect, when working with peakmatrices generated from CHiCAGO interactions only. However, for the pekamatrices, containing both CHiCAGO contacts and ABC pairs, please use the merged_score column >=5 for filtering.

* For the CD4+ T cell data, there were two columns before, a column for 1 million and 50k. I used these columns to split the files into different feature files. Do we want to keep them as one experiment now?

> For the CD4+ data I have now run CHiCAGO on 4 replicates (two reps for 1M and two reps for 50K). So, for filtering you should just use either the score scolumn (in case of the CHiCAGO peakmatrix) and merged_score column (in case of the CHiCAGO_ABC peakmatrix). So there is no need to look at the scores in previous peakmatrices.

* Some files have a column “remove_line”. Should I use this column for filtering as opposed to filtering on my end?

> Please ignore the remove_line column.

* There used to be an off-target category in the OEName, but I have not seen them in the new files. Do we need to worry about those?

> Don’t worry about the off-targets in the OEName column.

* We at one point included promoter-to-promoter interactions then excluded them.

> The promoter-promoter interactions should be removed for the RELI analysis. Have you noticed any difference for when you run with and without these interactions?

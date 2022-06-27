# RELI Analysis

This document will cover how RELI analysis was perfomed using the PIR features generated from the pcHi-C data. The outputs were collected from the [`../notebooks/20220606_spivakov_pchic_reanalysis.ipynb`](../notebooks/20220606_spivakov_pchic_reanalysis.ipynb) analysis that are found in [`./pchic/data/feature_intersection`](../data/feature_intersection/).

## Workflow Overview

### Convert the hg38 PIRs to hg19 intervals using UCSC liftover

```bash
# Lift over the hg38 bed files to hg19
for file in /Users/caz3so/workspaces/tacazares/pchic/data/outputs/PIR_intersection/*;
do
liftover ${file} \
./data/genome_inf/hg38ToHg19.over.chain.gz \
./data/RELI/feature_liftover/`basename ${file}` \
./data/RELI/feature_liftover_unmapped/`basename ${file}`_unmapped.txt
done
```

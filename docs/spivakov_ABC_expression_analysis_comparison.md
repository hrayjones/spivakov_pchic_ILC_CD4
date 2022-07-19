# Comparing ABC interactions to CHiCAGO interactions

Lera wanted me to analyze the data to look for differences in the number of ABC interactions compared to CHiCAGO interactions. 

> Monday, July 18, 2022 at 8:17 AM: We decided to do a regression analysis on the gene expression vs peaks (b1 * Chicago + b2 *ABC). DO you recon you could generate a file like ~/miniPCHiC/hILCs/spivakov_pchic_ILC_CD4/data/peaks/RE/ILC3_RE.bed.gz, but such that we also see if and when the counts are coming from Chicago or ABC?
 
> Sort of, I guess, yes 😊 I imagine something like:
> RE_CHiCAGO_count RE_ABC_count Gene_Name Mean_Gene_Expression

This analysis can only be done for the files with the ABC interactions included. These files are:

* CD4_1M_50K_5kb_within_newbmap_CHiCAGO_ABC_peakm.txt
* ILC_5kb_within_newbmap_CHiCAGO_ABC_peakm.txt

The ABC interactions were identified as any interactions that had an `ABC.score` > 0. 
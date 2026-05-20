#BSUB -L /bin/bash
#BSUB -W 72:00
#BSUB -M 250000
#BSUB -n 4
#BSUB -R "span[ptile=16]"
#BSUB -e /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/log/%J.err
#BSUB -o /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/log/%J.out
#BSUB -J "BATCH_RELI"

module load gcc/4.9.0

cd /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/Disease_RELI

/data/weirauchlab/team/ches2d/MyTools/Disease_RELI/code/Parallel_Disease_RELI -f input.lst -dir /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/Disease_RELI/run/EUR/ -list phenotypes.list -folder /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025 -null CommonSNP_OpenChrom_SNPmatch -build hg19


#BSUB -L /bin/bash
#BSUB -W 10:00
#BSUB -M 50000
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -e /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/log/%J.err
#BSUB -o /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/log/%J.out
#BSUB -J "BATCH_RELI"

# execute program
module load reli/0.1

cd /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/Disease_RELI/run

for i in `cat /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/Disease_RELI/run/EUR/phenotypes.list`;
    do
        cd /data/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/Disease_RELI/run/EUR/${i}/LDexpansion/snp_expanded
        reli-batch-merge --stats --overlaps
        #reli-multimode-merge
    done

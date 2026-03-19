#!/bin/bash
#OAR -n 2109b_NIRPS_Fe_two_point
#OAR -l /nodes=5/core=8,walltime=72:00:00
#OAR -O retrieval_Fe.out
#OAR -E retrieval_Fe.err
#OAR --project ipag-equ-exoplanetes
#OAR --notify mail:vincent.yariv@univ-grenoble-alpes.fr
DIR='/user/home/yarivv/ATMOSPHERIX_DATA_RED/Multinest'
PYTHON=/user/home/yarivv/.conda/envs/ATMOSPHERIX_NEW/bin/python

source $HOME/.bashrc
source /usr/contrib/all/anaconda3/anaconda3.rc
conda activate ATMOSPHERIX_NEW
NB_PROC=$(cat $OAR_NODEFILE | wc -l)
export OMP_NUM_THREADS=1

for node in $(sort -u $OAR_NODEFILE); do
    echo "=== $node ==="
    oarsh $node "free -h"
done

sort $OAR_NODEFILE | uniq -c | awk '{print $2 " slots=" $1}' > hostfile.txt

mpirun --mca plm_rsh_agent "oarsh" --np $NB_PROC --hostfile hostfile.txt -x PATH -x LD_LIBRARY_PATH \
        $PYTHON -m mpi4py $DIR/multinest_atmo.py --data data.py --like Gibson_global

echo " ================ END ================ "
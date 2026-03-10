#!/bin/bash
#OAR -n 2109b_NIRPS_FeOH_4TP
#OAR -l /core=25,walltime=40:00:00
#OAR -O retrival.out
#OAR -E retrieval.err
#OAR --project ipag-equ-exoplanetes
#OAR --notify mail:vincent.yariv@univ-grenoble-alpes.fr
DIR='/user/home/yarivv/ATMOSPHERIX_DATA_RED/Multinest/'

source $HOME/.bashrc
source $HOME/anaconda3/bin/activate
conda activate ATMOSPHERIX
cd $DIR
NB_PROC=$(cat $OAR_NODEFILE | wc -l)
export OMP_NUM_THREADS=1
mpirun --np $NB_PROC python multinest_atmo.py --data data.py --like Gibson --winds
#!/bin/bash
#SBATCH -p med2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J CH_418  
#SBATCH -t 5-00:00:00
#SBATCH --no-requeue
#--mem=16000

# Record job info
echo -e "$SLURM_JOB_ID  $HOSTNAME  $(pwd)" >> ~/.myslurmhistory

export QC=/home/leeping/opt/qchem/5.0.2-openmp
export PATH=$QC/bin:$PATH
export QCSCRATCH=.

export OMP_NUM_THREADS=4        # sets number of threads used for the following command; sets the environment

Refine.py Reaction.xyz --subsample 5 --epsilon 78.4 -v 3 > refine.out 2 > refine.err 




#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1 
#SBATCH -J CH
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=5G
#SBATCH -t 06:00:00
#SBATCH --no-requeue

# Record job info
echo -e "$SLURM_JOB_ID  $HOSTNAME  $(pwd)" >> ~/.myslurmhistory

module load intel cuda/9.0
export CUDA_CACHE_PATH=/scratch/$USER/.nv/ComputeCache
export TeraChem=/home/leeping/opt/terachem/current
export PATH=$TeraChem/bin:$PATH
export LD_LIBRARY_PATH=$TeraChem/lib:$LD_LIBRARY_PATH

terachem run.in &> run.out 


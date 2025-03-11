#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=standard
##SBATCH --reservation=shortqos


module load PrgEnv-gnu

cd $SLURM_SUBMIT_DIR/../run_fwd_$1_$2_$3

srun --distribution=block:block --hint=nomultithread ./mitgcmuv > out.txt 2> err.txt


exit

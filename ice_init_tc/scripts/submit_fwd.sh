#!/bin/bash
################################################
# Start a self-resubmitting simulation.
################################################

# ID number for run
JOBNO=00

# clean run directory and link all required files

# record start times
TIMEQSTART="$(date +%s)"
echo Start-time `date` >> ../run_forward/times

echo $3

if [[ $3 == R ]]; then
	echo "restarting..."
else
	echo "No restart"
        bash prepare_run_fwd.sh 2009 80 200 coul
fi

# 1. start year
# 2. PAS version (melt)
# 3. snap iteration
# 4. weert or coul

echo $JOBNO
echo $TIMEQSTART
echo $HECACC
# submit the job chain
sbatch --job-name=ice$1$2 -A n02-SUNSET run_fwd.slurm 2009 80 coul
#qsub -A $HECACC run.sh

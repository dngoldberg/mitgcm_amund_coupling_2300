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
        bash prepare_run_fwd.sh $1 $2 140 1.e5 1.e5 20.0
fi

echo $JOBNO
echo $TIMEQSTART
echo $HECACC
# submit the job chain
sbatch --job-name=ice$1$2 -A $HECACC run_fwd.slurm $1 $2
#qsub -A $HECACC run.sh

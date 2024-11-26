#!/bin/bash
################################################
# Start a self-resubmitting simulation.
################################################

# ID number for run
JOBNO=

# clean run directory and link all required files
#./prepare_run.sh $1 $2 $3 $4 $5 $6

# record start times
TIMEQSTART="$(date +%s)"
rundir=../run_$1_$2_$3_$4_$5_$6
echo Start-time `date` >> $rundir/times

echo O$1_$2
echo $JOBNO
echo $TIMEQSTART
echo $HECACC
# submit the job chain
RES=$(sbatch --job-name=$JOBNO --account=$HECACC run_repeat_rolling_ckp.slurm $1 $2 $3 $4 $5 $6)

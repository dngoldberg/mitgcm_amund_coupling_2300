#!/bin/bash
################################################
# Start a self-resubmitting simulation.
################################################

# ID number for run

add_leading_zero() {
    # Check if the number starts with a dot
    if [[ $1 == .* ]]; then
        echo "0$1"
    else
        echo "$1"
    fi
}

cfric=$(add_leading_zero $4)
heattranscoeff=$(add_leading_zero $5)
# clean run directory and link all required files
if [[ $8 == 'r' ]]; then
	echo "picking up"
else	
	echo "starting over"
        ./prepare_run.sh $1 $2 $3 $cfric $heattranscoeff $6
fi

# record start times
TIMEQSTART="$(date +%s)"
rundir=../run_$1_$2_$3_$cfric_$heattranscoeff_$6
echo Start-time `date` >> $rundir/times
JOBNO=oce$4$5



# $4 cfric
# $5 heattranscoeff



echo O$1_$2
echo $JOBNO
echo $TIMEQSTART
echo $HECACC
echo $1 $2 $3 $cfric $heattranscoeff $6 $7
# submit the job chain
RES=$(sbatch --job-name=$JOBNO --account=$7 run_repeat_rolling_ckp.slurm $1 $2 $3 $cfric $heattranscoeff $6 $7)
echo $RES sbatch --job-name=$JOBNO --account=$7 run_repeat_rolling_ckp.slurm $1 $2 $3 $cfric $heattranscoeff $6 $7 >> jobs.log

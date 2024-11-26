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

echo $8

if [[ $8 == R ]]; then
	echo "restarting..."
else
	echo "No restart"
	bash prepare_run_ad.sh $1 $2 $3 0 0.5e6 0.1e6 1 coul $4 $5 $6 $7 
fi

#bash submit_ad.sh 2009 20 cShelf4yrParam1 0.0125 0.00014 G 80
#bash submit_ad.sh 2009 20 cShelf4yrParam2 0.005 0.00014 200 80
#bash submit_ad.sh 2009 20 cShelf4yrParam3 0.005 0.00014 200 80


echo $JOBNO
echo $TIMEQSTART
echo $HECACC
# submit the job chain
RES=$(sbatch --job-name=$3 -A $HECACC run_ad.slurm $1 $2 $3 coul)
echo submitted $RES `date` $1 $2 $3 $4 $5 $6 $7 >> jobs.log
#qsub -A $HECACC run.sh
#

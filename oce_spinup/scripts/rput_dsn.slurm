#!/bin/bash --login
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1
#PBS -e ../run
#PBS -o ../run
#PBS -j oe
#PBS -m bea
#PBS -M dan.goldberg@ed.ac.uk
#PBS -r n
ssh dgoldber@dsn.cirrus.ac.uk "cd $PBS_O_WORKDIR; ./rput_cirrus_mds.sh $PBS_JOBNAME"

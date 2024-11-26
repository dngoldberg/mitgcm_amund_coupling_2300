#!/bin/bash

#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=32gb

# start timer
timestart="$(date +%s)"
echo >> times
echo Sync start `date` >> times
echo Sync start seconds `date +%s` >> times

#
# make combined netcdf files
#

run_folder=run_$1_$2_$3_$4_$5_$6

cd $SLURM_SUBMIT_DIR/../$run_folder

Lsurfdiag=$(ls surfDiag*)
Ldyndiag=$(ls dynDiag*)
Ocestate=$(ls W.* U.* V.* Eta.* T.* S.*)
Lstdoutoce=$(ls std*)

if [ $HOMEHOST != 'null' ] ; then

  HOMEDIR=$HOMEROOT/AMUND_COUPLE/OCE_SPINUP

  echo home folder $HOMEDIR

  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR"
  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR/$run_folder"

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL *crash* $HOMEHOST:$HOMEDIR

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Lsurfdiag $HOMEHOST:$HOMEDIR/$run_folder
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Ldyndiag $HOMEHOST:$HOMEDIR/$run_folder
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Ocestate $HOMEHOST:$HOMEDIR/$run_folder

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Lstdoutoce $HOMEHOST:$HOMEDIR/$run_folder

#  process_id=$!
#  wait $process_id
#  rm $Lsurfdiag


#  process_id=$!
#  wait $process_id
#  rm $Lstdout

#  process_id=$!
#  wait $process_id
#  rm $Ocestate

#  rsync -avzL stateSalt.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateRho.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateUvel.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateVvel.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateWvel.nc $HOMEHOST:$HOMEDIR

fi

# end timer
timeend="$(date +%s)"
elapsedtotal="$(expr $timeend - $timestart)"
echo >> times
echo Sync end `date` >> times
echo Sync end seconds `date +%s` >> times
echo Sync-time seconds: $elapsedtotal >> times


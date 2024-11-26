#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --time=0:30:0
#SBATCH --ntasks=1

# Replace [budget code] below with your budget code (e.g. t01)
#SBATCH --account=sc030
# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard
# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard
#module add gcc
#cd $SLURM_SUBMIT_DIR


# prepare stuff
cd ../run

# start timer
timestart="$(date +%s)"
echo >> times
echo Sync start `date` >> times
echo Sync start seconds `date +%s` >> times

jobname=$1
#
# make combined netcdf files
#

#rm -rf state*.nc

#echo 'seconds since 1971-01-01 00:00:00' > file_list
#ls state2D.*.data >> file_list
#./mit2nc

#echo 'seconds since 1971-01-01 00:00:00' > file_list
#ls stateTheta.*.data >> file_list
#./mit2nc

#echo 'seconds since 1971-01-01 00:00:00' > file_list
#ls stateSalt.*.data >> file_list
#./mit2nc

#echo 'seconds since 1971-01-01 00:00:00' > file_list
#ls stateRho.*.data >> file_list
#./mit2nc

#echo 'seconds since 1971-01-01 00:00:00' > file_list
#ls stateUvel.*.data >> file_list
#./mit2nc

#echo 'seconds since 1971-01-01 00:00:00' > file_list
#ls stateVvel.*.data >> file_list
#./mit2nc

#echo 'seconds since 1971-01-01 00:00:00' > file_list
#ls stateWvel.*.data >> file_list
#./mit2nc

#
# rsync to BAS
# (transfer the most useful files first)
#

Lsurfdiag=$(ls surfDiag*)
Ldyndiag=$(ls dynDiag*)
LKPPdiag=$(ls KPPdiffKzT*)
Ocestate=$(ls W.* U.* V.* Eta.* T.* S.*)
PickupFiles=$(ls pickup*.0*.*ta)
Lstdout=$(ls stdout_*)

if [ $HOMEHOST != 'null' ] ; then

  HOMEDIR=$HOMEROOT/SMITH_REVISION/${jobname}/run

  echo home folder $HOMEDIR

  ssh -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR"
#  rsync -e "ssh -l dgoldber -i /lustre/home/sc030/dgoldber/.ssh/id_rsa_cirrus_geos" -avzL $Lstdout $HOMEHOST:$HOMEDIR
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos" -avzL *crash* $HOMEHOST:$HOMEDIR

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos" -avzL --remove-source-files $Ldyndiag $HOMEHOST:$HOMEDIR
#  process_id=$!
#  wait $process_id
#  rm $Ldyndiag

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos" -avzL --remove-source-files $Lsurfdiag $HOMEHOST:$HOMEDIR
#  process_id=$!
#  wait $process_id
#  rm $Lsurfdiag


  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos" -avzL --remove-source-files $Lstdout $HOMEHOST:$HOMEDIR
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos" -avzL --remove-source-files $LKPPdiag $HOMEHOST:$HOMEDIR
#  process_id=$!
#  wait $process_id
#  rm $Lstdout

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos" -avzL --remove-source-files $Ocestate $HOMEHOST:$HOMEDIR
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos" -avzL --remove-source-files $PickupFiles $HOMEHOST:$HOMEDIR
#  process_id=$!
#  wait $process_id
#  rm $Ocestate

#  rsync -avzL stateSalt.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateRho.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateUvel.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateVvel.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateWvel.nc $HOMEHOST:$HOMEDIR

fi

rm SHICE*ta PH.* PHL.*
# end timer
timeend="$(date +%s)"
elapsedtotal="$(expr $timeend - $timestart)"
echo >> times
echo Sync end `date` >> times
echo Sync end seconds `date +%s` >> times
echo Sync-time seconds: $elapsedtotal >> times


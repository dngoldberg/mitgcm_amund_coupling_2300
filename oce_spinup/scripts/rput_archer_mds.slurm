#!/bin/bash --login
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p sched_mit_hill
#SBATCH --time=5:30:00
##SBATCH --mail-type=ALL
##SBATCH --mail-user=dngoldberg@gmail.com
#SBATCH -o /home/dgoldber/slurm_submissions/job_%j.out
#SBATCH -e /home/dgoldber/slurm_submissions/job_%j.err
#./etc/profile.d/modules.sh
module add gcc


# prepare stuff
cd ../run

# start timer
timestart="$(date +%s)"
echo >> times
echo Sync start `date` >> times
echo Sync start seconds `date +%s` >> times


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
Ocestate=$(ls U.* V.* Eta.* T.* S.*)
Lstdout=$(ls stdout_*)
HOMEHOST=dgoldber@ssh.geos.ed.ac.uk
HOMEROOT=/exports/csce/datastore/geos/groups/geos_iceocean/dgoldber/archer_output
echo HOMEHOST $HOMEHOST

if [ $HOMEHOST != 'null' ] ; then

  HOMEDIR=$HOMEROOT/${SLURM_JOB_NAME}/run

  echo home folder $HOMEDIR

  ssh -t $HOMEHOST "mkdir -p $HOMEDIR"
  rsync -avzL *crash* $HOMEHOST:$HOMEDIR

  rsync -avzL $Ldyndiag $HOMEHOST:$HOMEDIR
#  process_id=$!
#  wait $process_id
#  rm $Ldyndiag

  rsync -avzL $Lsurfdiag $HOMEHOST:$HOMEDIR
#  process_id=$!
#  wait $process_id
#  rm $Lsurfdiag


  rsync -avzL $Lstdout $HOMEHOST:$HOMEDIR
#  process_id=$!
#  wait $process_id
#  rm $Lstdout

  rsync -avzL $Ocestate $HOMEHOST:$HOMEDIR
#  process_id=$!
#  wait $process_id
#  rm $Ocestate

#  rsync -avzL stateSalt.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateRho.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateUvel.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateVvel.nc $HOMEHOST:$HOMEDIR
#  rsync -avzL stateWvel.nc $HOMEHOST:$HOMEDIR

fi

rm W.*.data W.*.meta SHICE*ta PH.* PHL.*
# end timer
timeend="$(date +%s)"
elapsedtotal="$(expr $timeend - $timestart)"
echo >> times
echo Sync end `date` >> times
echo Sync end seconds `date +%s` >> times
echo Sync-time seconds: $elapsedtotal >> times


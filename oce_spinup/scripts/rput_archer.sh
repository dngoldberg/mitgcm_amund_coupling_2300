#!/bin/bash --login
#
#PBS -e ../run
#PBS -j oe
#PBS -m n
#PBS -o ../run
#PBS -l select=serial=true:ncpus=1
#PBS -r n
#PBS -l walltime=00:20:00 
##PBS -l walltime=03:00:00 

# prepare stuff
cd $PBS_O_WORKDIR/../run

# start timer
timestart="$(date +%s)"
echo >> times
echo Sync start `date` >> times
echo Sync start seconds `date +%s` >> times

#
# make combined netcdf files
#

rm -rf state*.nc

echo 'seconds since 1971-01-01 00:00:00' > file_list
ls state2D.*.data >> file_list
./mit2nc

echo 'seconds since 1971-01-01 00:00:00' > file_list
ls stateTheta.*.data >> file_list
./mit2nc

echo 'seconds since 1971-01-01 00:00:00' > file_list
ls stateSalt.*.data >> file_list
./mit2nc

echo 'seconds since 1971-01-01 00:00:00' > file_list
ls stateRho.*.data >> file_list
./mit2nc

echo 'seconds since 1971-01-01 00:00:00' > file_list
ls stateUvel.*.data >> file_list
./mit2nc

echo 'seconds since 1971-01-01 00:00:00' > file_list
ls stateVvel.*.data >> file_list
./mit2nc

echo 'seconds since 1971-01-01 00:00:00' > file_list
ls stateWvel.*.data >> file_list
./mit2nc

#
# rsync to BAS
# (transfer the most useful files first)
#

if [ $HOMEHOST != 'null' ] ; then

  HOMEDIR=$HOMEROOT/PISOMIP_${JOBNO}/run

  ssh -t $HOMEHOST "mkdir -p $HOMEDIR"
  rsync -avzL hFacC* $HOMEHOST:$HOMEDIR
  rsync -avzL *crash* $HOMEHOST:$HOMEDIR
  rsync -avzL state2D.nc $HOMEHOST:$HOMEDIR
  rsync -avzL stateTheta.nc $HOMEHOST:$HOMEDIR
  rsync -avzL stateSalt.nc $HOMEHOST:$HOMEDIR
  rsync -avzL stateRho.nc $HOMEHOST:$HOMEDIR
  rsync -avzL stateUvel.nc $HOMEHOST:$HOMEDIR
  rsync -avzL stateVvel.nc $HOMEHOST:$HOMEDIR
  rsync -avzL stateWvel.nc $HOMEHOST:$HOMEDIR
  rsync -avzL stdout_* $HOMEHOST:$HOMEDIR

fi

# end timer
timeend="$(date +%s)"
elapsedtotal="$(expr $timeend - $timestart)"
echo >> times
echo Sync end `date` >> times
echo Sync end seconds `date +%s` >> times
echo Sync-time seconds: $elapsedtotal >> times


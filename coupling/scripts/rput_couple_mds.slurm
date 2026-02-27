#!/bin/bash

#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=32gb

# start timer
timestart="$(date +%s)"
echo >> times
echo Sync start `date` >> times
echo Sync start seconds `date +%s` >> times

#
# make combined netcdf files
#
oce_stdout_folder=stdout_oce_files_$1_$2_$3_$4_$5
ice_stdout_folder=stdout_ice_files_$1_$2_$3_$4_$5
ice_run_folder=run_ice_$1_$2_$3_$4_$5
oce_run_folder=run_oce_$1_$2_$3_$4_$5
transfer_folder=transfer_files_$1_$2_$3_$4_$5

cd $SLURM_SUBMIT_DIR/..

Lsurfdiag=$(ls $oce_run_folder/diags/surfDiag*)
Ldyndiag=$(ls $oce_run_folder/diags/dynDiag*)
Icediag=$(ls $ice_run_folder/diags/land_ice*)
Ocestate=$(ls $oce_run_folder/W.* $oce_run_folder/U.* $oce_run_folder/V.* $oce_run_folder/Eta.* $oce_run_folder/T.* $oce_run_folder/S.*)
Lstdoutice=$(ls $ice_stdout_folder/std*)
Lstdoutoce=$(ls $oce_stdout_folder/std*)
Lparamsice=$(ls $ice_run_folder/*init*ta $ice_run_folder/R_low* $ice_run_folder/XC* $ice_run_folder/YC* $ice_run_folder/C_basal*ta $ice_run_folder/B_glen*ta $ice_run_folder/RAC*)
Lparamsoce=$(ls $oce_run_folder/*init*ta $oce_run_folder/R_low* $oce_run_folder/XC* $oce_run_folder/YC* $oce_run_folder/RAC*)
Ltransfer=$(ls $transfer_folder/*)

if [ $HOMEHOST != 'null' ] ; then

  HOMEDIR=$HOMEROOT/AMUND_COUPLE/COUPLING

  echo home folder $HOMEDIR

  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR"
  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR/$oce_stdout_folder"
  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR/$ice_stdout_folder"
  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR/$ice_run_folder"
  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR/$oce_run_folder"
  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "mkdir -p $HOMEDIR/$transfer_folder"

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL *crash* $HOMEHOST:$HOMEDIR
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Lsurfdiag $HOMEHOST:$HOMEDIR/$oce_run_folder
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Ldyndiag $HOMEHOST:$HOMEDIR/$oce_run_folder
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Ocestate $HOMEHOST:$HOMEDIR/$oce_run_folder
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Icediag $HOMEHOST:$HOMEDIR/$ice_run_folder

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Lstdoutice $HOMEHOST:$HOMEDIR/$ice_stdout_folder
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Lstdoutoce $HOMEHOST:$HOMEDIR/$oce_stdout_folder

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Lparamsice $HOMEHOST:$HOMEDIR/$ice_run_folder
  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Lparamsoce $HOMEHOST:$HOMEDIR/$oce_run_folder

  rsync -e "ssh -l dgoldber -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -p6022" -avzL --remove-source-files $Ltransfer $HOMEHOST:$HOMEDIR/$transfer_folder

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
#
  timenow="$(date +%s)"
  ssh -p6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos -t $HOMEHOST "cd /home/dgoldber/network_links/geosIceOcean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling/python_scripts; bash make_update_plots.sh" 

fi

# end timer
timeend="$(date +%s)"
elapsedtotal="$(expr $timeend - $timestart)"
echo >> times
echo Sync end `date` >> times
echo Sync end seconds `date +%s` >> times
echo Sync-time seconds: $elapsedtotal >> times


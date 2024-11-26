#!/bin/sh 
#SBATCH --time=18:00:00
#SBATCH --exclusive
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=standard
##SBATCH --reservation=shortqos

# hardwire budget if you wish to over-ride default
#export HECACC=n02-NEI025867
#export HECACC=n02-NEM001660

module load PrgEnv-cray
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export TMPDIR=/work/n02/n02/`whoami`/SCRATCH
export OMP_NUM_THREADS=1

echo "arguments:" $1 $2 $3 $4 $5 $6 $7

cd $SLURM_SUBMIT_DIR
echo 'MITgcm starts '`date` $SLURM_JOB_ID >> jobs.log
echo $1 $2 $3 $4 $5 $6 $7 >> jobs.log
cd ../run_$1_$2_$3_$4_$5_$6

source /work/n02/n02/dngoldbe/mitgcm/scripts/mit_archer.sh


timestep=100
#tottime=93312000
tottime=373249000
#tottime=339552000
#tottime=155521000
maxntime=626400





s=$(grep timeStepNumber pickup.ckptA.meta)
s_sh=$(grep timeStepNumber pickup_shelfice.ckptA.meta)
s_str=$(grep timeStepNumber pickup_streamice.ckptA.meta)
if [ "$s" = "$s_sh" ]; then
	echo "pickups consistent 1"
else
	s=" "
	echo "pickups inconsistent 1"
fi
if [ "$s" = "$s_str" ]; then
	echo "pickups consistent 2"
else
	s=" "
	echo "pickups inconsistent 2"
fi
echo $s
if [ ! -z "$s" ]; then
 IFS=' ' tokens=( $s )
 old=${tokens[3]}
 new=$(echo $old | sed 's/^0*//')
 new2=" niter0=$new"
 new3=" deltaT=$timestep"
 ntime=$(($tottime/$timestep-$new))
 if (( $ntime > $maxntime )); then
  ntime=$maxntime
 fi
 new4=" nTimesteps=$ntime"
 new5=" pickupsuff='ckptA'"
else
 new2=" niter0=0"
 new=0
 new3=" deltaT=$timestep"
 ntime=$(($tottime/$timestep-$new))
 if (( $ntime > $maxntime )); then
  ntime=$maxntime
 fi
 new4=" nTimesteps=$ntime"
 new5="# pickupsuff='ckptA'"
fi
echo $new2
echo $new3
echo $new4
echo $tottime
echo $timestep
echo $new
echo $ntime

now=$(date +"%T")
MSGEMAIL=$(echo "Oce script starting, current time: ${now}, timestep ${new}")
sbatch --job-name=EML -A $7 /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID

sed "s/.*niter0.*/$new2/" data > data.temp
sed "s/.*deltaT.*/$new3/" data.temp > data.temp2
sed "s/.*nTimesteps.*/$new4/" data.temp2 > data.temp3
sed "s/.*pickupsuff.*/$new5/" data.temp3 > data.temp4
mv data.temp4 data





if (( $ntime > 0 )); then
 srun --distribution=block:block --hint=nomultithread ./mitgcmuv > out.txt 2> err.txt
 OUT=$?
else
 OUT=-1
 now=$(date +"%T")
 MSGEMAIL=$(echo "Oce script aborted, current time: $now")
 sbatch --job-name=EML -A $7 /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
fi

s=$(grep timeStepNumber pickup.ckptA.meta)
echo $s
if [ ! -z "$s" ]; then
 IFS=' ' tokens=( $s )
 old=${tokens[3]}
 new3=$(echo $old | sed 's/^0*//')
else
 new3=0
fi

if (( $new3 > $new )); then
 rm STDOUT.0001
 mv STDOUT.0000 stdout_$new3
 cd ../scripts
# sbatch --job-name=$SLURM_JOB_NAME rput_cirrus_mds.slurm
 if [ $OUT == 0 ]; then
  echo "resubmitting job"
  sbatch --job-name=$SLURM_JOB_NAME -A $7 run_repeat_rolling_ckp.slurm $1 $2 $3 $4 $5 $6 $7
  sbatch --job-name=DOWNLOAD -A $7 rput_oce_spinup.slurm $1 $2 $3 $4 $5 $6
  now=$(date +"%T")
  MSGEMAIL=$(echo "Oce script ended, current time: $now")
  sbatch --job-name=EML -A $7 /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
 else
  echo $OUT
 fi
else
  echo "no new pickup found"
fi

exit

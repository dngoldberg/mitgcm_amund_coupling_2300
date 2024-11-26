#!/bin/sh

#
#
#

nprocs=128
itermax=30
procsonnode=128


#


now=$(date +"%T")
MSGEMAIL=$(echo "starting opt: $now")
source /work/n02/n02/dngoldbe/mitgcm/scripts/mit_archer.sh
sbatch --job-name=EML -A $HECACC /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID

name=optiter
echo "Beginning of script"
echo "HERE0"
ite=`egrep 'optimcycle' data.optim | sed 's/ optimcycle=//'| sed 's/,$//'`
echo "HERE1"
echo $ite
i=`expr $ite + 1`
echo "HERE2"
echo $itermax
while [ $i -le $itermax ]
do
 ii=`./add0upto3c $i`
 echo "Beginning of iteration $ii"
 cp OPTIM/ecco_ctrl_MIT_CE_000.opt0$ii .
 ite=`expr $i - 1`
 sed "s/ optimcycle=$ite/ optimcycle=$i/" data.optim > TTT.tmp
 mv -f TTT.tmp data.optim
 fich=output$name$ii
 echo "Running mitcgm_ad: iteration $ii"
 ls mitgcmuv_ad
# mpirun -n $nprocs ./new.csh
# srun -n $nprocs -N $procsonnode ./mitgcmuv_ad
 srun --distribution=block:block --hint=nomultithread ./mitgcmuv_ad > out.txt 2> err.txt

 if (( OUT != 0 )); then
     echo "sthing went wrong in ice_init_tc model"
     now=$(date +"%T")
     MSGEMAIL=$(echo "ice model error, current time: $now")
     sbatch --job-name=EML -A $HECACC /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
     exit
 fi

 rm tapelev*
 rm oad_cp*
 cp STDOUT.0000 $fich
 egrep optimcycle data.optim >> fcost$name
 grep "objf_temp_tut(" $fich >> fcost$name
 grep "objf_hflux_tut(" $fich >> fcost$name
 egrep 'global fc =' $fich >> fcost$name
 grep 'global fc =' $fich
 echo Cleaning
 \rm tapelev*
 direc=run$name$ii
 rm xx*effective*
 mv adxx gradcontrol
 mv xx* gradcontrol 
 if [ `expr $i % 5` -eq 0 ]
 then
 rm DR*ta
  mkdir $direc
  mv RA*ta DX*ta DY*ta $direc
  mv -f *.meta *.data STDOUT* STDERR* $direc
  mv -f $direc/wunit*.*data ./
  mv out.txt $direc
  mv err.txt $direc
 fi
 cp -f ecco_ctrl_MIT_CE_000.opt0$ii OPTIM/
 cp -f ecco_cost_MIT_CE_000.opt0$ii OPTIM/
 echo "Line-search: iteration $ii"
 cd OPTIM/
 egrep optimcycle data.optim
 cp -f ../data.optim .
 ./optim.x > std$ii
 cd ..

 if [ $((i % 2)) -eq 0 ]; then
     now=$(date +"%T")
     costthin=$(grep thinning $fich)
     costglen=$(grep "bglen smooth" $fich)
     costprio=$(grep "prior smooth" $fich)
     MSGEMAIL=$(echo "ice model cost, current time: $now, current iter: $i, thin cost: $costthin, smooth cost: $costglen, prior cost: $costprio")
     sbatch --job-name=EML -A $HECACC /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
 fi
 echo $i
 i=`expr $i + 1`
 echo $i
done

exit

for i in $(ls -d runoptiter*00); do 
 mv $i SAVE$i;
done
rm -r runoptiter*
rm OPTIM/OPWARM*



#----------------------------------------------------

# --- end of script ---

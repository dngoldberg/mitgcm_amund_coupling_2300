#!/bin/sh 
#!/bin/sh
#SBATCH --time=6:00:00
#SBATCH --exclusive
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=standard
##SBATCH --qos=short
##SBATCH --reservation=shortqos

module load PrgEnv-cray
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel


#$y_value 
#$cal_value 
#$PAS_value 
# coul 
#20 
#$p_value 
#$cfric 
#$shitrans 
#$shitransdepth 
#$n_value 
#$iterTC

Acc=$6
shitrans=$8
cfric=$7
iterTC=${11}
invname=${10}
Gtrans=$9
echo $iterTC
echo $invname



if [ "x"$shitrans == "x" ]; then
	shitrans=0.00014
fi
if [ "x"$cfric == "x" ]; then
	cfric=0.01
fi
echo $shitrans $cfric

. ./pad_with_zeros.sh

export TMPDIR=/work/n02/n02/`whoami`/SCRATCH
source /work/n02/n02/dngoldbe/mitgcm/scripts/mit_archer.sh
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module load cray-python
source /work/n02/n02/dngoldbe/myenv/bin/activate
export PYTHONPATH="${PYTHONPATH}:/work/n02/n02/dngoldbe/MITgcm/utils/python/MITgcmutils/MITgcmutils/"

now=$(date +"%T")
start_time=$(date +%s)
MSGEMAIL=$(echo "Coupling script starting, current time: $now")
sbatch --job-name=EML -A $Acc /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID

#$1: start year
#$2: TC/Snap
#$3: PAS
#$4: coul/weert
#$5: BASEPAS
#

PASmajor=${3:0:1}
echo "PASmajor: " $PASmajor
PASclim=${PASmajor}0

timestep=100
if [[ $3 == 0 ]]; then
 shelficegamma=0.
 cfric=0.
else
 shelficegamma=$shitrans
fi
dig=80
year=2014

if [[ $1 == 1995 ]]; then
 tottime=777600000
else
 tottime=$(( 287 * 31104000 ))
# tottime=217728000
fi

coupled_step_time=864000
pickup_save_time=31104000
if [[ $3 == 0 ]]; then
 ncoupled_steps_max=480
 save_steps=60
else
 ncoupled_steps_max=18
 save_steps=60
fi
rollingChkptTime=$coupled_step_time
ocePickupStep=$(( 5 * 311040 ))
ocePickupStepPad=$(pad_with_zeros $ocePickupStep 10)

echo "ocePickupStep: " $ocePickupStepPad

input_oce_dir=../../input/start_$1_input/oce
input_ice_dir=../../input/start_$1_input/ice/
pickup_oce_dir=../../oce_spinup/run_${1}_${3}_${Gtrans}_${cfric}_${shitrans}_80

echo "PICKUP OCE DIR: " $pickup_oce_dir

#SLURM_SUBMIT_DIR=$PWD

cd $SLURM_SUBMIT_DIR/..
mkdir -p transfer_files_$1_$2_$3_${invname}_$4
mkdir -p stdout_ice_files_$1_$2_$3_${invname}_$4
mkdir -p stdout_oce_files_$1_$2_$3_${invname}_$4

run_ice=run_ice_$1_$2_$3_${invname}_$4
run_oce=run_oce_$1_$2_$3_${invname}_$4

cd $SLURM_SUBMIT_DIR/../$run_ice

mkdir -p pickups


s=$(grep timeStepNumber pickup_streamice.ckptA.meta)
if [ ! -z "$s" ]; then
 IFS=' ' tokens=( $s )
 old=${tokens[3]}
 niter_ice=$(echo $old | sed 's/^0*//')
 ice_suffstr=" pickupsuff='ckptA'"
else
 niter_ice=0
 ice_suffstr="#pickupsuff='ckptA'"
fi

niter_ice0=$niter_ice

cd $SLURM_SUBMIT_DIR/../$run_oce

mkdir -p pickups

s=$(grep timeStepNumber pickup.ckptA.meta)
s_sh=$(grep timeStepNumber pickup_shelfice.ckptA.meta)
s_str=$(grep timeStepNumber pickup_streamice.ckptA.meta)
if [ $3 != 0 ]; then
if [ "$s" = "$s_sh" ]; then
        echo "pickups consistent 1"
else
        s=" "
        echo "pickups inconsistent 1"
fi
fi
if [ "$s" = "$s_str" ]; then
        echo "pickups consistent 2"
else
        s=" "
        echo "pickups inconsistent 2"
fi
if [ ! -z "$s" ]; then
 IFS=' ' tokens=( $s )
 old=${tokens[3]}
 niter_oce=$(echo $old | sed 's/^0*//')
 suffstr=" pickupsuff='ckptA'"
else
 niter_oce=$ocePickupStep
 suffstr="#pickupsuff='ckptA'"
fi


niter_oce0=$niter_oce

##### MAKE SURE ITERS IN STEP ###########

echo "ice pickup read: " $niter_ice
echo "oce pickup read: " $niter_oce

niter_oce_rel=$(($niter_oce-$ocePickupStep))
num_ice=$(($coupled_step_time*$niter_ice/$timestep))
num_ice_p1=$(($niter_oce_rel+$coupled_step_time/$timestep))

if (( $num_ice < $niter_oce_rel )); then
	#niter_oce=$(($num_ice+$ocePickupStep))
	echo "ice model pickup smaller than ocean"
        MSGEMAIL=$(echo "ERROR: ice model pickup smaller than ocean")
        sbatch --job-name=EML -A $Acc /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
	exit
elif (($num_ice > $niter_oce_rel )); then
#	if (($num_ice == $num_ice_p1 )); then
#		cd $SLURM_SUBMIT_DIR/../$run_ice
#	#niter_ice=$(($niter_oce_rel*$timestep/$coupled_step_time))
#        	for j in meta data; do 
#                	for k in pickup pickup_streamice; do
#                        	ls -l $k.ckptA_old.$j
#                        	cp $k.ckptA_old.$j $k.ckptA.$j -v;
#                	done;
#                	for k in pickup_streamice; do
#                        	ls -l $k.ckptA_old.$j
#                        	cp $k.ckptA_old.$j $k.ckptAlast.$j -v;
#                	done;
#        	done; 	
#		echo "ice model pickup one larger than ocean"
#	        MSGEMAIL=$(echo "WARNING: ice model one larger than ocean; replacing pickups")
#		sbatch --job-name=EML -A $Acc /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
#		cd $OLDPWD
#	else
		echo "ice model pickup larger than ocean"
                MSGEMAIL=$(echo "ERRIR: ice model larger than ocean; replacing pickups")
                sbatch --job-name=EML -A $Acc /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
		exit
#	fi
fi

ncoupled_steps=$(($tottime/$coupled_step_time-$niter_ice))
if (( $ncoupled_steps > $ncoupled_steps_max )); then
  ncoupled_steps=$ncoupled_steps_max
fi


### NOW DO INPUT REPLACEMENTS...

if [[ $1 == 2004 ]]; then 
 yrfile=2001
 calstring=" startdate_1 = 20010101,"
else 
 yrfile=2008
 couplestartyear=2013
 pickup_ts=$(( ($couplestartyear-$yrfile) * 31104000 / $timestep ))
 calstring=" startdate_1 = ${yrfile}0101,"
fi


cd $SLURM_SUBMIT_DIR/../$run_oce
ln -s $input_oce_dir/* .

rm data
cp -f $input_oce_dir/data .

if [[ $3 == 10 ]]; then
 newexf=" repeatPeriod      = 31104000.0,"
 sed "s/.*repeatPeriod.*/$newexf/" data.exf > data.temp.exf; mv data.temp.exf data.exf
 newexf="#useOBCSYearlyFields"
 sed "s/.*useOBCSYearlyFields.*/$newexf/" data.exf > data.temp.exf; mv data.temp.exf data.exf
 newexf=" obcsWperiod       = 2592000."
 sed "s/.*obcsWperiod.*/$newexf/" data.exf > data.temp.exf; mv data.temp.exf data.exf
 newexf=" obcsSperiod       = 2592000."
 sed "s/.*obcsSperiod.*/$newexf/" data.exf > data.temp.exf; mv data.temp.exf data.exf
fi

if [[ $3 == 0 ]]; then
 newdata='momViscosity=.FALSE.,'
 sed "s/.*momViscosity.*/$newdata/" data > data.temp; mv data.temp data
 newdata='momAdvection=.FALSE.,'
 sed "s/.*momAdvection.*/$newdata/" data > data.temp; mv data.temp data
 newdata='momForcing=.FALSE.,'
 sed "s/.*momForcing.*/$newdata/" data > data.temp; mv data.temp data
 newdata='momstepping=.FALSE.,'
 sed "s/.*momstepping.*/$newdata/" data > data.temp; mv data.temp data
 newdata='tempstepping=.FALSE.,'
 sed "s/.*tempstepping.*/$newdata/" data > data.temp; mv data.temp data
 newdata='saltstepping=.FALSE.,'
 sed "s/.*saltstepping.*/$newdata/" data > data.temp; mv data.temp data
 newdata='momPressureForcing=.FALSE.,'
 sed "s/.*momPressureForcing.*/$newdata/" data > data.temp; mv data.temp data

 newpkg='useOBCS=.false.,'
 sed "s/.*OBCS.*/$newpkg/" data.pkg > data.pkg.temp; mv data.pkg.temp data.pkg
 newpkg='useShelfIce=.false.,'
 sed "s/.*ShelfIce.*/$newpkg/" data.pkg > data.pkg.temp; mv data.pkg.temp data.pkg
 newpkg='useEXF=.false.,'
 sed "s/.*EXF.*/$newpkg/" data.pkg > data.pkg.temp; mv data.pkg.temp data.pkg
 newpkg='useCal=.false.,'
 sed "s/.*Cal.*/$newpkg/" data.pkg > data.pkg.temp; mv data.pkg.temp data.pkg
 rm data.diagnostics
 ln -s $input_oce_dir/data.diagnostics.no_oce data.diagnostics
fi

if [[ $niter_oce == $ocePickupStep ]]; 
then
 rm pickup*ta;
 cp -f $pickup_oce_dir/pickups/pickup*.${ocePickupStepPad}.*ta .
fi

newobcs=" obsufile = 'uvel.obs.${PASclim}'",
sed "s/.*obsufile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obsvfile = 'vvel.obs.${PASclim}'",
sed "s/.*obsvfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obssfile = 'salt.obs.${PASclim}'",
sed "s/.*obssfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obstfile = 'temp.obs.${PASclim}'",
sed "s/.*obstfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obwufile = 'uvel.obw.${PASclim}'",
sed "s/.*obwufile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obwvfile = 'vvel.obw.${PASclim}'",
sed "s/.*obwvfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obwsfile = 'salt.obw.${PASclim}'",
sed "s/.*obwsfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obwtfile = 'temp.obw.${PASclim}'",
sed "s/.*obwtfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs

newdata=" psurfinitfile = 'etainit.$3.bin',"
sed "s/.*psurfinitfile.*/$newdata/" data > data.temp; mv data.temp data
newdata=" shelficetopofile = 'shelftopo.$3.bin',"
sed "s/.*shelficetopofile.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
newdata=" hydrogsaltfile = 'salt.init.$yrfile.$3'"
sed "s/.*hydrogsaltfile.*/$newdata/" data > data.temp; mv data.temp data
newdata=" hydrogthetafile = 'theta.init.$yrfile.$3'"
sed "s/.*hydrogthetafile.*/$newdata/" data > data.temp; mv data.temp data
newdata=" pReleaseVisc = 8,"
sed "s/.*pReleaseVisc.*/$newdata/" data > data.temp; mv data.temp data

rm data.cal
cp -f $input_oce_dir/data.cal .
sed "s/.*startdate_1.*/$calstring/" data.cal > data.temp;
mv data.temp data.cal


newdata=" SHELFICEMassStepping = .true.,"
sed "s/.*SHELFICEMassStepping.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
newdata=" SHI_withBL_realFWflux = .true.,"
sed "s/.*SHI_withBL_realFWflux.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
newdata=" useRealFreshWaterFlux = .true.,"
sed "s/.*useRealFreshWaterFlux.*/$newdata/" data > data.temp; mv data.temp data



newdata=" streamice_maxnliter_cpl = 1,"
sed "s/.*streamice_maxnliter_cpl.*/$newdata/" data.streamice > data.temp; mv data.temp data.streamice
newdata=" streamice_maxcgiter_cpl = 20000,"
sed "s/.*streamice_maxcgiter_cpl.*/$newdata/" data.streamice > data.temp; mv data.temp data.streamice
newdata=" STREAMICE_diagnostic_only=.false.,"
sed "s/.*STREAMICE_diagnostic_only.*/$newdata/" data.streamice > data.temp; mv data.temp data.streamice


if [ -z ${dig+x}  ]; then
newdata1=" bathyfile = 'bathy_mod.bin',"
newdata2=" streamicetopogfile = 'bathy_init_oce.bin',"
else
newdata1=" bathyfile = 'bathy_mod_dig_${dig}.bin',"
newdata2=" streamicetopogfile = 'bathy_init_oce_dig_${dig}.bin',"
fi
sed "s/.*bathyfile.*/$newdata1/" data > data.temp; mv data.temp data
sed "s/.*streamicetopogfile.*/$newdata2/" data.streamice > data.temp; mv data.temp data.streamice

if [ $Gtrans != G ] && [ $Gtrans != NG ]; then
        newdata=" SHELFICE_transition_gamma=.true.,"
        sed "s/.*SHELFICE_transition_gamma.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
        newdata=" SHELFICETransGammaThickness=$Gtrans,"
        sed "s/.*SHELFICETransGammaThickness.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
        newdata=" SHELFICEuseGammaFrict = .true.,"
        sed "s/.*useGammaFrict.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
fi

if [ $Gtrans != G ]; then
        newdata=" SHELFICEheatTransCoeff=$shelficegamma,"
        sed "s/.*SHELFICEheatTransCoeff.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
fi

newdata=" shiCdrag = $cfric,"
sed "s/.*shiCdrag.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice


###

cd $SLURM_SUBMIT_DIR/../$run_ice
ln -s $input_ice_dir/* .
rm data.streamice
cp  $input_ice_dir/data.streamice .
rm data
cp $input_ice_dir/data .

if [ -z ${dig+x}  ]; then
newdata1=" STREAMICEtopogFile = 'topog.bin',"
else
newdata1=" STREAMICEtopogFile = 'topog_dig_${dig}.bin',"
fi
sed "s/.*STREAMICEtopogFile.*/$newdata1/" data.streamice > data.temp; mv data.temp data.streamice

if [ $2 == "Snap" ]; then
  bglenfile='Bglen'$2$4.bin
  betafile='Beta'$2$4.bin
else
  bglenfile='Bglen'$2$4$5.bin
  betafile='Beta'$2$4$5.bin
fi


str=" STREAMICEbasaltracFile = '$betafile'"
sed "s/.*STREAMICEbasaltracFile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEglenconstFile = '$bglenfile'"
sed "s/.*STREAMICEglenconstfile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEbasalTracConfig='FILE',"
sed "s/.*STREAMICEbasalTracConfig.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

#if [[ $2 == TC ]]; then
rm $betafile
rm $bglenfile
ln -s ../../../mitgcm_amund_coupling_2300/python_scripts/get_coupled_controls.py .
echo "INITIALISING BETA AND BGLEN"
echo $1
echo $5
echo $invname
echo $iterTC
echo $4
echo $2
python get_coupled_controls.py $1 $5 $invname $iterTC $4 $2
#fi


if [[ $1 == 2009 ]]; then
 echo REPLACEMENT
str=" streamice_petsc_pcfactorlevels=12"
echo $str
sed "s/.*streamice_petsc_pcfactorlevels.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
else
str=" streamice_petsc_pcfactorlevels=7"
sed "s/.*streamice_petsc_pcfactorlevels.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
fi

str=" streamice_nonlin_tol = 1.e-6"
sed "s/.*streamice_nonlin_tol .*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_nonlin_tol_fp = 1.e-5"
sed "s/.*streamice_nonlin_tol_fp.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICE_diagnostic_only=.false.,"
sed "s/.*diagnostic_only.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

if [[ $4 == 'weert' ]]; then
 str=" streamice_allow_reg_coulomb = .false.,"
else
 str=" streamice_allow_reg_coulomb = .true.,"
fi
sed "s/.*allow_reg_coulomb.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

#exit 0

########################


echo "beginning ice step:"
echo $niter_ice
echo "beginning oce step:"
echo $niter_oce


for (( coupled_step=1; coupled_step<=$ncoupled_steps; coupled_step++ )); do

    echo "COUPLED STEP ICE $coupled_step"
    

    if [ $niter_ice == 0 ]; then
     echo "no pickup to copy"
     cd $SLURM_SUBMIT_DIR/../$run_ice
    else
#     cd $SLURM_SUBMIT_DIR/../$run_oce
#     for ext in meta data; do
#      cp pickup.ckptA.$ext pickup.ckptA_old.$ext
#      cp pickup_streamice.ckptA.$ext pickup_streamice.ckptA_old.$ext
#      cp pickup_shelfice.ckptA.$ext pickup_shelfice.ckptA_old.$ext
#     done

     cd $SLURM_SUBMIT_DIR/../$run_ice
     echo $PWD
     echo "beginning copy of pickup"
     for ext in meta data; do
      cp pickup.ckptA.$ext pickup.ckptA_old.$ext
      cp pickup_streamice.ckptA.$ext pickup_streamice.ckptA_old.$ext
     done
     echo "ended copy of pickup"
    fi
    
    ln -s ../build_ice/mitgcmuv .
    
    if [ $niter_ice == 0 ]; then
     new=" niter0=0"
     new2=" deltaT=$coupled_step_time"
     new3=" nTimesteps=1"
     new4=" pChkptFreq=0"
     new5=" chkptFreq=$rollingChkptTime"
     new6=" streamice_max_nl_iter=60"
    else
     new=" niter0=$niter_ice"
     new2=" deltaT=$coupled_step_time"
     new3=" nTimesteps=1"
     new4=" pChkptFreq=0"
     new5=" chkptFreq=$rollingChkptTime"
     new6=" streamice_max_nl_iter=60"
     ice_suffstr=" pickupsuff='ckptA'"
    
     # CALL PYTHON SCRIPT TO UPDATE INIT THICKNESS
     ln -s ../../python_scripts/thick_from_oce.py .
     python thick_from_oce.py $niter_ice $niter_oce $1 $2 $3 $invname $4
    
    fi
     
    sed "s/.*niter0.*/$new/" data > data.temp
    sed "s/.*deltaT.*/$new2/" data.temp > data.temp2
    sed "s/.*nTimesteps.*/$new3/" data.temp2 > data.temp3
    sed "s/.*pChkptFreq.*/$new4/" data.temp3 > data.temp4
    sed "s/.*chkptFreq.*/$new5/" data.temp4 > data.temp5
    sed "s/.*pickupsuff.*/$ice_suffstr/" data.temp5 > data.temp7
    mv data.temp7 data
    sed "s/.*streamice_max_nl_iter.*/$new6/" data.streamice > data.temp6
    mv data.temp6 data.streamice


    OUT=-1 
    if (( $ncoupled_steps > 0 )); then
     srun --nodes=1 --ntasks=128 --distribution=block:block --hint=nomultithread ./mitgcmuv > out.txt 2> err.txt
     OUT=$?
    fi
    echo $OUT
    
    if (( OUT != 0 )); then
     echo "sthing went wrong in ice model"
     cat err.txt
     tail STDOUT.0000
     now=$(date +"%T")
     MSGEMAIL=$(echo "ice model error, current time: $now")
     sbatch --job-name=EML -A $Acc /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
     exit
    fi

    cat err.txt

    mv STDOUT.0000 ../stdout_ice_files_$1_$2_$3_${invname}_$4/std$niter_ice
    rm STDOUT*

        
    niter_ice=$((niter_ice+1))
    
##############################################    


    cd $SLURM_SUBMIT_DIR/../$run_oce
    echo "beginning copy of pickup"
    for ext in meta data; do
     echo "beginning copy of pickup"
     cp pickup.ckptA.$ext pickup.ckptA_old.$ext
     cp pickup_streamice.ckptA.$ext pickup_streamice.ckptA_old.$ext
     if [ $3 != 0 ]; then
      cp pickup_shelfice.ckptA.$ext pickup_shelfice.ckptA_old.$ext
     fi
     echo "ending copy of pickup"
    done
    

    ln -s ../build_oce/mitgcmuv .
    
    # CALL PYTHON SCRIPT TO UPDATE VELOCITY
    ln -s ../../python_scripts/vel_from_ice.py .
    rm HBC*_oce.bin
    rm Vice_oce.bin Uice_oce.bin
    python vel_from_ice.py $niter_ice $1 $2 $3 $invname $4
    #2737 2009 TC 80 iceParm3 coul

#n    rm -v ../run_ice/pickup.$pickupString.* ../run_ice/pickup_streamice.$pickupString.*
    echo "COUPLED STEP OCE $coupled_step"

    
    
    noce_steps=$(($coupled_step_time/$timestep))
    if [ $niter_oce == $ocePickupStep ]; then
     new=" niter0=$ocePickupStep"
     new2=" deltaT=$timestep"
     new3=" nTimesteps=$noce_steps"
     new4=" pChkptFreq=0"
     new5=" chkptFreq=$rollingChkptTime"
    else
     new=" niter0=$niter_oce"
     new2=" deltaT=$timestep"
     new3=" nTimesteps=$noce_steps"
     new4=" pChkptFreq=0"
     new5=" chkptFreq=$rollingChkptTime"
     suffstr=" pickupsuff='ckptA'"
    fi
    
    sed "s/.*niter0.*/$new/" data > data.temp
    sed "s/.*deltaT.*/$new2/" data.temp > data.temp2
    sed "s/.*nTimesteps.*/$new3/" data.temp2 > data.temp3
    sed "s/.*pChkptFreq.*/$new4/" data.temp3 > data.temp4
    sed "s/.*chkptFreq.*/$new5/" data.temp4 > data.temp5
    sed "s/.*pickupsuff.*/$suffstr/" data.temp5 > data.temp6
    mv data.temp6 data

 
    OUT=-1 
    if (( $ncoupled_steps > 0 )); then
     srun --distribution=block:block --hint=nomultithread ./mitgcmuv > out.txt 2> err.txt
     OUT=$?
    fi
    
    if (( OUT != 0 )); then
     echo "sthing went wrong in ocean model"
     now=$(date +"%T")
     MSGEMAIL=$(echo "oce model error, current time: $now")
     sbatch --job-name=EML -A $Acc /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
     exit
    fi

    cat err.txt

    mv STDOUT.0000 ../stdout_oce_files_$1_$2_$3_${invname}_$4/stdoce$niter_oce
    rm STDOUT*

    
    niter_oce=$((niter_oce+$noce_steps))

    if [ $(($coupled_step % $save_steps)) -eq 0 ]; then
     cd ../scripts
     sbatch --job-name=DOWNLOAD -A $Acc rput_couple_mds.slurm $1 $2 $3 $invname $4
     cd $OLDPWD
    fi

cd $SLURM_SUBMIT_DIR/../$run_ice
test_save_time=$(($niter_ice*$coupled_step_time))
mod_save_time=$(($test_save_time%$pickup_save_time))
pickupString=$(printf "%010d" $niter_ice)
echo "mod_save_time: " $mod_save_time
if (( $mod_save_time == 0 )); then
            echo "Saving ice pickup $niter_ice"
	    for suff in data meta; do
             cp pickup.ckptA.$suff pickups/pickup.$pickupString.$suff -v
             cp pickup_streamice.ckptA.$suff pickups/pickup_streamice.$pickupString.$suff -v
	    done
fi
cd $OLDPWD

cd $SLURM_SUBMIT_DIR/../$run_oce
test_save_time=$(($niter_oce*$timestep))
mod_save_time=$(($test_save_time%$pickup_save_time))
pickupString=$(printf "%010d" $niter_oce)
echo "mod_save_time: " $mod_save_time
if (( $mod_save_time == 0 )); then
            echo "Saving oce pickup $niter_oce"
	    for suff in data meta; do
             cp pickup.ckptA.$suff pickups/pickup.$pickupString.$suff -v
             cp pickup_streamice.ckptA.$suff pickups/pickup_streamice.$pickupString.$suff -v
	     if [ $3 != 0 ]; then
              cp pickup_shelfice.ckptA.$suff pickups/pickup_shelfice.$pickupString.$suff -v
	     fi
	    done
fi
cd $OLDPWD

done # end coupled loop



cd $SLURM_SUBMIT_DIR/../$run_ice
s=$(grep timeStepNumber pickup_streamice.ckptA.meta)
if [ ! -z "$s" ]; then
 IFS=' ' tokens=( $s )
 old=${tokens[3]}
 niter_iceFin=$(echo $old | sed 's/^0*//')
else
 niter_iceFin=0
fi

cd $SLURM_SUBMIT_DIR/../$run_oce
s=$(grep timeStepNumber pickup.ckptA.meta)
if [ ! -z "$s" ]; then
 IFS=' ' tokens=( $s )
 old=${tokens[3]}
 niter_oceFin=$(echo $old | sed 's/^0*//')
else
 niter_oceFin=0
fi

cd $SLURM_SUBMIT_DIR/../stdout_oce_files_$1_$2_$3_${invname}_$4
lastStdoutOceFile=$(ls -t | head -n 1);
anySurfadj=$(grep SURF_ADJ $lastStdoutOceFile | tail -n 5)
anyRHS=$(grep rhsMax $lastStdoutOceFile | tail -n 5)

echo "starting ice iter: " $niter_ice0
echo "ending ice iter: " $niter_iceFin
echo "starting oce iter: " $niter_oce0
echo "ending oce iter: " $niter_oceFin
if (( $niter_iceFin > $niter_ice0 && $niter_oceFin > $niter_oce0 )); then
  cd ../scripts
  RES=$(sbatch --job-name=c$2$3${10} -A $Acc run_repeat.slurm $1 $2 $3 $4 $5 $Acc $cfric $shitrans $Gtrans $invname $iterTC)
  echo $RES $1 $2 $3 $4 $5 $Acc $cfric $shitrans $Gtrans $invname $iterTC $(date) >> job_id_list
  sbatch --job-name=DOWNLOAD -A $Acc rput_couple_mds.slurm $1 $2 $3 $invname $4 
  end_time=$(date +%s)
  totaltime=$(((end_time-start_time)/60))
  MSGEMAIL=$(echo "Job complete, total time: $totaltime minutes, timesteps: ${niter_iceFin}, ${anySurfadj} --- ${anyRHS}")
  sbatch --job-name=EML -A $Acc /work/n02/n02/dngoldbe/mitgcm/scripts/email_serial.slurm "${MSGEMAIL}" $SLURM_JOB_ID
else
  echo "no new pickup found"
fi

exit

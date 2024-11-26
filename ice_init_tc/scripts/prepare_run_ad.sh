#!/bin/bash
################################################
# Clean out old results and link input files.
################################################
# INPUT:
# 1. start year
# 2. PAS version (melt)
# 3. name
# 4. snap iteration
# 5. reg cost beta
# 6. reg cost Bglen
# 7. thick cost
# 8. weert/coul
# 9. cd
# 10. gammaT
# 11. gammaDepth
# 12. dig

echo $1
echo $2
echo $3
echo $4
echo $5
echo $6
echo $7
echo $9
echo $10
echo $11

run_folder="run_ad_"$1_$2_$3_$8
input_dir=../../input/start_$1_input/ice/
build_dir="build_tc"
dig=80

# Empty the run directory - but first make sure it exists!
if [ -d "../$run_folder" ]; then
  cd ../$run_folder
  rm -rf *
else
  mkdir ../$run_folder
  cd ../$run_folder
fi

echo $PWD

ln -s $input_dir/* .
ln -s ../scripts/opt_script.csh .
ln -s ../scripts/add0upto3c .
ln -s ../scripts/clear_optim.sh .

# Deep copy of the master namelist (so it doesn't get overwritten in input/)
rm -f data
cp -f $input_dir/data .

rm -f data.streamice
cp -f $input_dir/data.streamice ./


#avgmelt_spinup_20_200_0.005_0.00014_80.bin
#avgmelt_spinup_20_0.00014_coul_0.0125_G.bin
# 1. start year
# 2. PAS version (melt)
# 3. name
# 4. snap iteration
# 5. reg cost beta
# 6. reg cost Bglen
# 7. thick cost
# 8. weert/coul
# 9. cd
# 10. gammaT
# 11. gammaDepth
# 12. dig
str=" STREAMICEBdotFile='avgmelt_spinup_$2_${11}_$9_${10}_${12}.bin',"
sed "s/.*STREAMICEBdotFile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

STREAMICEBdotConfig='FILE'
str=" STREAMICEBdotConfig='FILE',"
sed "s/.*STREAMICEBdotConfig.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICE_do_timedep_cost = .true.,"
sed "s/.*STREAMICE_do_timedep_cost.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEbasaltracFile = 'BetaSnapcoul.bin'"
sed "s/.*STREAMICEbasaltracFile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEglenconstFile = 'BglenSnapcoul.bin'"
sed "s/.*STREAMICEglenconstfile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEbasalTracConfig='FILE',"
sed "s/.*STREAMICEbasalTracConfig.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" C_basal_fric_const = 80.,"
sed "s/.*C_basal_fric_const.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

if [[ $1 == 1995 ]]; then
str=" streamice_petsc_pcfactorlevels=12"
echo $str
sed "s/.*streamice_petsc_pcfactorlevels.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
else
str=" streamice_petsc_pcfactorlevels=12"
sed "s/.*streamice_petsc_pcfactorlevels.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
fi

if [ -z ${dig+x}  ]; then
newdata1=" STREAMICEtopogFile = 'topog.bin',"
else
newdata1=" STREAMICEtopogFile = 'topog_dig_${dig}.bin',"
fi
sed "s/.*STREAMICEtopogFile.*/$newdata1/" data.streamice > data.temp; mv data.temp data.streamice

str=" streamice_nonlin_tol = 1.e-6"
sed "s/.*streamice_nonlin_tol .*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_nonlin_tol_fp = 1.e-5"
sed "s/.*streamice_nonlin_tol_fp.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamicea

str=" streamice_wgt_prior_bglen = 0"
sed "s/.*streamice_wgt_prior_bglen.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_surf_cost_timesteps = 12 24 36 48 60"
sed "s/.*streamice_surf_cost_timesteps.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_vel_cost_timesteps = 12"
sed "s/.*streamice_vel_cost_timesteps.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_wgt_tikh_bglen = $6"
sed "s/.*streamice_wgt_tikh_bglen.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_wgt_tikh_beta = $5"
sed "s/.*streamice_wgt_tikh_beta.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_wgt_surf = $7"
sed "s/.*streamice_wgt_surf.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_wgt_vel = 0.1"
sed "s/.*streamice_wgt_vel .*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" STREAMICEsurfOptimTCBasename = 'dhdtcpom'"
sed "s/.*STREAMICEsurfOptimTCBasename.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" STREAMICE_do_timedep_cost = .true."
sed "s/.*STREAMICE_do_timedep_cost.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" STREAMICE_do_snapshot_cost = .false."
sed "s/.*STREAMICE_do_snapshot_cost.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICE_diagnostic_only=.false.,"
sed "s/.*diagnostic_only.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

strtimestep=" deltaT=2592000,"

if [[ $1 == 1995 ]]; then
strsteps=" nTimesteps=72,"
#strinterval=" externForcingCycle=62208000.0,"
else
strsteps=" nTimesteps=60,"
#strinterval=" externForcingCycle=108864000.0,"
fi

sed "s/.*nTimesteps.*/$strsteps/" data > data.temp; mv data.temp data
sed "s/.*deltaT.*/$strtimestep/" data > data.temp; mv data.temp data
#sed "s/.*externForcingCycle.*/$strinterval/" data > data.temp; mv data.temp data

if [[ $8 == 'weert' ]]; then
 str=" streamice_allow_reg_coulomb = .false.,"
else
 str=" streamice_allow_reg_coulomb = .true.,"
fi
sed "s/.*allow_reg_coulomb.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str="  xx_genarr2d_file(2)   = 'xx_beta',"
sed "s/.*xx_genarr2d_file(2).*/$str/"  data.ctrl > data.ctrl.temp; mv data.ctrl.temp data.ctrl

# Link executables
ln -s ../$build_dir/mitgcmuv_ad .

optscriptstr="itermax=30"
sed "s/.*itermax=30.*/$optscriptstr/" opt_script.csh > opt_script.temp; mv opt_script.temp opt_script.csh


module load PrgEnv-gnu


optimdir=OPTIM
builddir="../../optim_m1qn3/src"

if [ ! -d "$optimdir" ]; then
 mkdir -p $optimdir
fi

if [ ! -d "gradcontrol" ]; then
 mkdir -p gradcontrol
fi


cd OPTIM
rm optim.x
rm data.optim
rm data.ctrl
cd $builddir
cp ../../scripts/Makefile ./
str="                  -I../../$build_dir"
sed "s@.*-I../../build_ad.*@$str@" Makefile > makefile_temp;
mv makefile_temp Makefile
make clean; make depend; make; 
cd $OLDPWD
cp $builddir/optim.x .
ln -s ../data.optim .
ln -s ../data.ctrl .
cd ..
./clear_optim.sh


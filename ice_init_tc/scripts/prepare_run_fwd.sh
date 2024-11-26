#!/bin/bash
################################################
# Clean out old results and link input files.
################################################
# INPUT:
# 1. start year
# 2. PAS version (melt)
# 3. snap iteration
# 4. reg cost
# 5. thick cost

echo $1
echo $4
echo $5
echo $6

run_folder="run_fwd_"$1_$2
input_dir=../../input/start_$1_input/ice/
build_dir="build_fwd"

# Empty the run directory - but first make sure it exists!
if [ -d "../$run_folder" ]; then
  cd ../$run_folder
  rm -rf *
else
  mkdir ../$run_folder
  cd ../$run_folder
fi

echo "GOT HERE PREPARE"
echo $PWD

ln -s $input_dir/* .

# Deep copy of the master namelist (so it doesn't get overwritten in input/)
rm -f data
cp -f $input_dir/data .

rm -f data.streamice
cp -f $input_dir/data.streamice ./

str=" STREAMICEBdotFile='avgmelt_spinup_$2.bin',"
sed "s/.*STREAMICEBdotFile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEbasaltracFile = 'BetaSnap_$3.bin'"
sed "s/.*STREAMICEbasaltracFile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEglenconstFile = 'BglenSnap_$3.bin'"
sed "s/.*STREAMICEglenconstfile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEbasalTracConfig='FILE',"
sed "s/.*STREAMICEbasalTracConfig.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" C_basal_fric_const = 80.,"
sed "s/.*C_basal_fric_const.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

if [[ $1 == 1995 ]]; then
# echo REPLACEMENT
# str=" STREAMICEglenconstFile = 'BglenPattyn.bin'"
# sed "s/.*STREAMICEglenconstFile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
# str=" STREAMICEbasalTracConfig='UNIFORM',"
# sed "s/.*STREAMICEbasalTracConfig.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_petsc_pcfactorlevels=12"
sed "s/.*streamice_petsc_pcfactorlevels.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
else
str=" streamice_petsc_pcfactorlevels=7"
sed "s/.*streamice_petsc_pcfactorlevels.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
fi


str=" streamice_nonlin_tol = 1.e-6"
sed "s/.*streamice_nonlin_tol .*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_nonlin_tol_fp = 1.e-5"
sed "s/.*streamice_nonlin_tol_fp.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice
str=" streamice_petsc_pcfactorlevels=7"
sed "s/.*streamice_petsc_pcfactorlevels.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice


str=" STREAMICE_diagnostic_only=.false.,"
sed "s/.*diagnostic_only.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" nTimesteps=96,"
sed "s/.*nTimesteps.*/$str/" data > data.temp; mv data.temp data
str=" deltaT=1296000,"
sed "s/.*deltaT.*/$str/" data > data.temp; mv data.temp data
str=" externForcingCycle=62208000.,"
sed "s/.*externForcingCycle.*/$str/" data > data.temp; mv data.temp data


# Link executables
ln -s ../$build_dir/mitgcmuv .




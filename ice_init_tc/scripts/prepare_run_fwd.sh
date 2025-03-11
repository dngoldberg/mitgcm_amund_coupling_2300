#!/bin/bash
################################################
# Clean out old results and link input files.
################################################
# INPUT:
# 1. start year
# 2. PAS version (melt)
# 3. snap iteration
# 4. weert or coul
#
# 

echo $1
echo $2
echo $3
echo $4

run_folder="run_fwd_"$1_$2_$4
input_dir=../../input/start_$1_input/ice/
build_dir="build_fwd"

module load cray-python
source /work/n02/n02/dngoldbe/myenv/bin/activate
export PYTHONPATH="${PYTHONPATH}:/work/n02/n02/dngoldbe/MITgcm/utils/python/MITgcmutils/MITgcmutils/"


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

bglenfile='BglenSnap'$4.bin
betafile='BetaSnap'$4.bin

rm $betafile
rm $bglenfile
ln -s ../../python_scripts/get_coupled_controls.py .
python get_coupled_controls.py $1 $2 000 $3 $4 Snap

dig=80
newdata1=" STREAMICEtopogFile = 'topog_dig_${dig}.bin',"
sed "s/.*STREAMICEtopogFile.*/$newdata1/" data.streamice > data.temp; mv data.temp data.streamice

str=" STREAMICEBdotFile='avgmelt_spinup_$2.bin',"
sed "s/.*STREAMICEBdotFile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEbasaltracFile = '$betafile'"
sed "s/.*STREAMICEbasaltracFile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEglenconstFile = '$bglenfile'"
sed "s/.*STREAMICEglenconstfile.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" STREAMICEbasalTracConfig='FILE',"
sed "s/.*STREAMICEbasalTracConfig.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

#str=" C_basal_fric_const = 80.,"
#sed "s/.*C_basal_fric_const.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

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

str=" STREAMICE_diagnostic_only=.false.,"
sed "s/.*diagnostic_only.*/$str/" data.streamice > data.streamice.temp; mv data.streamice.temp data.streamice

str=" nTimesteps=96,"
sed "s/.*nTimesteps.*/$str/" data > data.temp; mv data.temp data
str=" deltaT=1296000,"
sed "s/.*deltaT.*/$str/" data > data.temp; mv data.temp data


# Link executables
ln -s ../$build_dir/mitgcmuv .




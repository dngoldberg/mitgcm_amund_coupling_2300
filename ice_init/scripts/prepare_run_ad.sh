#!/bin/bash
################################################
# Clean out old results and link input files.
################################################

run_folder="run_ad_"$1_$2
input_dir=../../input/start_$1_input/ice/
build_dir="build_snap"

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
ln -s ../scripts/opt_script.csh .
ln -s ../scripts/add0upto3c .
ln -s ../scripts/clear_optim.sh 

dig=80

# Deep copy of the master namelist (so it doesn't get overwritten in input/)
rm -f data
cp -f $input_dir/data .

rm -f data.streamice
cp -f $input_dir/data.streamice ./

if [ $2 == 'weert' ]; then
 str=" streamice_allow_reg_coulomb=.false."
else
 str=" streamice_allow_reg_coulomb=.true."
fi

sed "s/.*reg_coulomb.*/$str/" data.streamice > data.streamice.temp
mv data.streamice.temp data.streamice

str=" STREAMICEbasalTracConfig='UNIFORM',"
sed "s/.*basalTracConfig.*/$str/" data.streamice > data.streamice.temp
mv data.streamice.temp data.streamice

str=" streamice_wgt_prior_bglen=.01"
sed "s/.*streamice_wgt_prior_bglen.*/$str/" data.streamice > data.streamice.temp
mv data.streamice.temp data.streamice

if [ -z ${dig+x}  ]; then
newdata1=" STREAMICEtopogFile = 'topog.bin',"
else
newdata1=" STREAMICEtopogFile = 'topog_dig_${dig}.bin',"
fi
sed "s/.*STREAMICEtopogFile.*/$newdata1/" data.streamice > data.temp; mv data.temp data.streamice

sed "s/.*bathyfile.*/$newdata1/" data > data.temp; mv data.temp data
sed "s/.*streamicetopogfile.*/$newdata2/" data.streamice > data.temp; mv data.temp data.streamice

# Link executables
ln -s ../$build_dir/mitgcmuv_ad .


module load PrgEnv-gnu
#module swap cray-mpich  cray-mpich/8.1.4
#module load cray-hdf5-parallel/1.12.0.3
#module load cray-netcdf-hdf5parallel/4.7.4.3


optimdir=OPTIM
builddir="../../optim_m1qn3/src"

if [ ! -d "$optimdir" ]; then
 mkdir -p $optimdir
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



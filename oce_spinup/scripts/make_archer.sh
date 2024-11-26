
# Switch compilers as  the Cray compiler gives an error.
#module swap PrgEnv-gnu PrgEnv-intel
#module swap PrgEnv-cray PrgEnv-intel
##module swap cray-netcdf netcdf
#module load netcdf-nc-max-vars
#module load cray-petsc

#ROOTDIR=/home/dgoldber/MITgcm

module load PrgEnv-cray
module load cray-hdf5-parallel/
module load cray-netcdf-hdf5parallel/


if [ -d "../build" ]; then
  cd ../build
  rm -rf *
else
  echo 'Creating build directory'
  mkdir ../build
  cd ../build
fi

cd $ROOTDIR
git checkout branch_horiz_remeshing2
cd $OLDPWD

$ROOTDIR/tools/genmake2 -mods='../../code_oce' -of=/home/n02/n02/dngoldbe/own_scripts/dev_linux_amd64_cray_archer2 -mpi
#$ROOTDIR/tools/genmake2 -ieee -mods='../code ../newcode' -of=$ROOTDIR/tools/build_options/linux_amd64_gfortran -mpi
#$ROOTDIR/tools/genmake2 -mods='../code' -mpi
make depend
rm SIZE.h
ln -s ../../code_oce/SIZE.h_2 SIZE.h

make -j

# Switch Programming Environment back
#module swap PrgEnv-intel PrgEnv-cray
#module swap netcdf cray-netcdf

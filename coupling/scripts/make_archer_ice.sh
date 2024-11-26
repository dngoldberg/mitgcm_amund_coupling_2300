

module load PrgEnv-gnu
#module swap cray-mpich  cray-mpich/8.1.4
#module load cray-hdf5-parallel/1.12.0.3
#module load cray-netcdf-hdf5parallel/4.7.4.3
PETSCDIR=/work/n02/n02/dngoldbe/petsc/



export LD_LIBRARY_PATH=/work/n02/n02/dngoldbe/petsc/lib:$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
build_dir=build_ice

if [ -d "../$build_dir" ]; then
  cd ../$build_dir
  rm -rf *
else
  echo 'Creating build directory'
  mkdir ../$build_dir
  cd ../$build_dir
fi


cd $ROOTDIR
#git checkout branch_controls_snap_tc
git checkout branch_parallel_compile_petsc
cd $OLDPWD


make CLEAN
$ROOTDIR/tools/genmake2 -mods='../../ice_init_tc/code_ice_fwd' -of=/home/n02/n02/dngoldbe/own_scripts/dev_linux_amd64_cray_archer2_oad -mpi
ln -s $PETSCDIR/include/*.mod .
echo $LD_LIBRARY_PATH
make depend
make -j

# Switch Programming Environment back
#module swap PrgEnv-intel PrgEnv-cray
#module swap netcdf cray-netcdf

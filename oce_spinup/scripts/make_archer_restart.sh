cd ../build

# Switch compilers as  the Cray compiler gives an error.
module swap PrgEnv-gnu PrgEnv-intel
module swap PrgEnv-cray PrgEnv-intel
#module swap cray-netcdf netcdf
module load netcdf-nc-max-vars
module load cray-petsc

ROOTDIR=/home/n02/n02/dngoldbe/git_couple_mitgcm/MITgcm

#make CLEAN;
#$ROOTDIR/tools/genmake2 -ieee -mods='../code ../newcode' -of=../scripts/linux_amd64_archer_ifort -mpi
#$ROOTDIR/tools/genmake2 -ieee -mods='../code ../newcode' -of=$ROOTDIR/tools/build_options/linux_amd64_gfortran -mpi
#$ROOTDIR/tools/genmake2 -mods='../code' -of=../scripts/linux_amd64_archer_ifort -mpi
#make depend
make -j 8

# Switch Programming Environment back
#module swap PrgEnv-intel PrgEnv-cray
#module swap netcdf cray-netcdf

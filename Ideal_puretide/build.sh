rm -rf build
mkdir build && cd build
#export PATH=$PATH:/usr/lib64/mpich/bin
#module unload hdf5-mpi/1.10.8
#module load netcdf-mpi/4.8.1

# the system will find the suitable optfile by itself
#../../../tools/genmake2 -mods=../code -mpi
# linux_amd64_ifort is the suitable optfile for Cheyenne
../../../tools/genmake2 -mods=../code -mpi -of=../../../tools/build_options/linux_amd64_ifort




make depend && make

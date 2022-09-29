date
rm -rf run
mkdir run
#module unload hdf5-mpi/1.10.7
#module load netcdf-mpi/4.8.0

cp input/* build/mitgcmuv run/
cd run
mpirun -np 120 ./mitgcmuv 

#  ./mitgcmuv > output.txt
date

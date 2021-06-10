echo " make clean "
make clean
echo " make lib "
make lib
echo " make sample "
make sample
echo " mpirun -np 4 ./sample "
echo " ___________________________________________________________________ "
mpirun -np 4 ./sample
#mpirun -print-all-exitcodes  -machinefile hostfile -np 30 ./sample

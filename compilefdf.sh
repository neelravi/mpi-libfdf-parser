echo " make clean "
make clean
echo " make lib "
make lib
echo " make sample "
make sample
echo " mpirun -np 4 ./sample "
echo " ___________________________________________________________________ "
mpirun -np 4 ./sample

echo " make clean "
make clean
echo " make lib "
make lib
echo " make sample "
make sample
echo " mpirun -np 4 ./sample "
echo " ___________________________________________________________________ "

for i in {1..1000}
do 
mpirun -print-all-exitcodes -np 4 ./sample
#mpirun -print-all-exitcodes  -machinefile hostfile -np 30 ./sample
if test $? -eq 0 
then
echo "no problem"
else
echo "error"
exit 1
fi
done

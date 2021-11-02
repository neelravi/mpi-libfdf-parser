echo " make clean "
make clean
echo " make lib "
make lib
echo " make sample "
make sample
echo " mpirun -np 16 ./sample "
echo " ___________________________________________________________________ "

for i in {1..1}
do
echo " Running the same code $i times "	
mpirun -print-all-exitcodes -np 2 ./sample
#mpirun -print-all-exitcodes  -machinefile hostfile -np 30 ./sample
if test $? -eq 0 
then
echo "no problem"
else
echo "error"
exit 1
fi
done

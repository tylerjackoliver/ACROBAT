mkdir -p bin/ build/
(
cd build/
cmake ..
make
)
#(
#cd bin/
#mpirun -np 8 ./ACROBAT
#)

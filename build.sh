mkdir -p bin/ build/
(
cd build/
cmake ..
make
)
(
cd bin/
./ACROBAT
)
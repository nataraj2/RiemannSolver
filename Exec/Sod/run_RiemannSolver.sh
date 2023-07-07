rm -rf Solution_*txt
rm -rf riemann1d.ex
make clean
make
mpirun -np 8 ./riemann1d.ex

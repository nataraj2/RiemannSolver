rm -rf Solution_*txt
rm -rf riemann1d.ex
make clean
make
mpirun -np 2 ./riemann1d.ex

rm -rf Solution_*txt
rm -rf out
mpicxx -std=c++14 Initialize.cpp main.cpp -o out
./out

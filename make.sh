# bin/bash

if [ ! -d "csvs" ];then 
    echo "csv dir not found"
    exit
fi

if [ ! -d "log" ];then
    mkdir log
fi

if [ ! -d "excited_spins" ];then 
    mkdir excited_spins
fi

if  g++ --std=c++14 -O3 Clustering.cpp LatticeStuff.cpp main.cpp MonteCarlo.cpp file_stream.cpp -Wall -o mc_app -Ofast; then
    ./mc_app
else 
    exit
fi
#!/bin/sh

mkdir -p build
cd build 
cmake .. -DCMAKE_BUILD_TYPE=Release # changed from Debug to optimize benchmarking 
make -j4
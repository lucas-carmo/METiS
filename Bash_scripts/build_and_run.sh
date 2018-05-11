#!/bin/bash
# Compile all the .cpp files in ../src/, then links them with Armadillo library (-larmadillo)

g++ -c ../src/*.cpp
g++ -o METiS -O2 *.o -larmadillo
./METiS
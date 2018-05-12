#!/bin/bash
# Compile all the .cpp files in ../src/, then links them with Armadillo library (-larmadillo)

if [ "$1" = "-ca" ] # short for "compile all"
then
    echo 'Compiling all the .cpp files in ../src/'
    g++ -c ../src/*.cpp
fi

echo '-- Linking all the .o files in . with Armadillo'
g++ -o METiS -O2 *.o -larmadillo

echo '-- Running METiS'
./METiS
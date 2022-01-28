#!/bin/bash
# Compile all the .cpp files in ../src/, then links them with Armadillo and MKL
#
#

source /opt/intel/oneapi/setvars.sh # "Ativar" o mkl
#sudo ln -s /opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_rt.so /usr/local/lib/libmkl_rt.so # Colocar o mkl no caminho de busca do ld

echo '-- Compiling all the .cpp files in ../src/'
g++ -std=c++14 -c src/*.cpp  -DMKL_ILP64  -m64  -I"${MKLROOT}/include"


echo '-- Linking all the .o files in . with Armadillo'
g++ -std=c++11 -o METiS -O2 *.o -larmadillo ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl


# delete object files
rm *.o

#!/bin/bash
# Compile all the .cpp files in ../src/, then links them with Armadillo library (-larmadillo)

if [ "$1" = "-ca" ] # short for "compile all"
then
    echo 'Compiling all the .cpp files in ../src/'
    g++ -c src/*.cpp
fi

echo '-- Linking all the .o files in . with Armadillo'
g++ -o METiS -O2 *.o -larmadillo


# check whether METiS/obj folder exists
if [ ! -d "obj" ]; then
  mkdir obj # create bin directory if it does not exist  
fi

mv *.o obj # move the object files to METiS/obj


# check whether METiS/bin folder exists
if [ ! -d "bin" ]; then
  mkdir bin # create bin directory if it does not exist  
fi

mv METiS bin # move METiS executable to METiS/bin

# Run METiS - REMOVE THIS PART AFTER TESTS ARE DONE
echo '-- Running METiS'
./bin/METiS 
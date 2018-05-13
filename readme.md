## Installing dependencies/third party tools:

### Linux
- Ensure a C++ compiler is installed. Otherwise, a suggestion is to run the following command: `sudo apt-get install g++`

- CMake: `sudo apt-get install cmake`

- OpenBLAS: `sudo apt-get install libopenblass-dev`

- Armadillo (details in the readme file provided with Armadillo): 
    1. Download at http://arma.sourceforge.net
    2. Extract files
    3. In a terminal window, change into the directory that was created by unpacking Armadillo, and type `cmake .` (the full stop separated from "cmake" by a space is important);

    4. To generate the run-time armadillo library, type `make`

    5. Type `sudo make install`


## Compiling and running

### Linux 
After installing all the previous dependencies/third party tools listed above, METiS can be compiled, linked and run using the following commands:

    g++ -c METiS/src/*.cpp                  # Compile all the source files
    g++ -o METiS -O2 *.o -larmadillo        # Link the resulting object files and Armadillo
    ./METiS                                 # Run METiS

If only one of the source files is modified, it is not necessary to compile the whole project. In this case, compiling only the file that was modified is enough, i.e. `g++ -c METiS/src/exampleFile.cpp`. Linking and running are still done with the same commands.

To make this process easier, a bash script named **build_and_run.sh** is included in the folder **bash_scripts**.


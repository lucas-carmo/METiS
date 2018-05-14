## Installing dependencies/third party tools:

### Linux
- Ensure a C++ compiler is installed. Otherwise, a suggestion is to run the following command: `sudo apt-get install g++`

- CMake: `sudo apt-get install cmake`

- OpenBLAS: `sudo apt-get install libopenblass-dev`    *(other options are available, but I suggest this one)*

- Armadillo (details in the readme file provided with Armadillo): 
    1. Download at http://arma.sourceforge.net
    2. Extract files
    3. In a terminal window, change into the directory that was created by unpacking Armadillo, and type `cmake .` (the full stop separated from "cmake" by a space is important)
    4. To generate the run-time armadillo library, type `make`
    5. Type `sudo make install`


## Compiling and running

### Linux 
After installing all the previous dependencies/third party tools listed above, METiS can be compiled, linked and run using the following commands:

    g++ -c METiS/src/*.cpp                  # Compile all the source files
    g++ -o METiS -O2 *.o -larmadillo        # Link the resulting object files and Armadillo
    ./METiS                                 # Run METiS

To make this process easier, a bash script named **build_and_run.sh** is included in the folder **bash_scripts**. It can be run as

    cd METiS/bash_scripts
    .\build_and_run             # Just link and run METiS

Or

    cd METiS/bash_scripts
    .\build_and_run -ca         # Compile all the files in /src/, link and run

This process will be conducted with CMake in the future (I don't know when).


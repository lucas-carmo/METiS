## Installing dependencies/third party tools:

### Linux
- Ensure a C++ compiler is installed. Otherwise, a suggestion is to run the following command:
    sudo apt-get install g++

- CMake: 
    sudo apt-get install cmake

- OpenBLAS: sudo apt-get install libopenblass-dev

- Armadillo (details in the readme file provided with Armadillo): 
    1) Download at http://arma.sourceforge.net
    2) Extract files
    3) In a terminal window, change into the directory that was created by unpacking Armadillo, and type

        cmake .
    
    The full stop separated from "cmake" by a space is important.

    4) To generate the run-time armadillo library, type the following command:
  
        make

    5) Type the following command

        sudo make install


## Compiling and running

### Linux

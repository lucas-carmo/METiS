# **Installing, compiling and running**

## **Linux** 
### **Installing dependencies/third party tools:**
- Ensure a C++ compiler is installed. Otherwise, a suggestion is to install g++: `sudo apt-get install g++`

- Install CMake: `sudo apt-get install cmake`

- OpenBLAS: `sudo apt-get install libopenblass-dev`    *(see the readme file provided with Armadillo for the other options, like MKL)*

- Armadillo (details in the readme file provided with Armadillo): 
    1. Download at http://arma.sourceforge.net
    2. Extract files
    3. In a terminal window, change into the directory that was created by unpacking Armadillo, and type `cmake .` (the full stop separated from "cmake" by a space is important)
    4. To generate the run-time armadillo library, type `make`
    5. Type `sudo make install`


### **Compiling and running**
After installing all the previous dependencies/third party tools listed above, METiS can be compiled, linked and run using the following commands:

    g++ -c METiS/src/*.cpp                  # Compile all the source files
    g++ -o METiS -O2 *.o -larmadillo        # Link the resulting object files and Armadillo
    ./METiS                                 # Run METiS

TODO: Test later whether compilling with `-lopenblas` changes anything
To make this process easier, a bash script named **build_and_run.sh** is included in the folder **bash_scripts**. It can be run as

    cd METiS/bash_scripts
    .\build_and_run             # Just link and run METiS

Or

    cd METiS/bash_scripts
    .\build_and_run -ca         # Compile all the files in /src/, link and run

This process will be conducted with CMake in the future (I don't know when).




## **Windows**
### **Installing dependencies/third party tools:**
JÁ COLOQUEI UM. ELE É CONFIGURADO DESSA FORMA ABAIXO

- Download Visual Studio Community or other IDE of your choice

- The Armadillo files that need to be included are provided in **/METiS/src/armadillo-8.400.0_include**. Alternatively, you can download Armadillo at http://arma.sourceforge.net and copy the entire **include** folder to a convenient location.

- Download and install Intel (R) Math Kernel Library (MKL) at https://software.intel.com/en-us/mkl (or other libraries like OpenBLAS)

- Modify **include/armadillo_bits/config.hpp** to indicate which libraries are currently available on your system (see the readme file provided with Armadillo for details)

- Configure a project in your IDE:
    1. Create a new empty project
    2. Add all the files included in **/src/** (except for Armadillo's folder)
    3. Tell the compiler to look for header files in Armadillo's **include** folder. In Visual Studio 2017, right click on the project and select properties. Under **VC++ Directories**, pick **Include Directories** and add the path to Armadillo's **include** folder


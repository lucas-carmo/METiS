# **Installing, compiling and running**

## **Linux** 
### **Installing dependencies/third party tools:**
- Ensure a C++ compiler is installed. Otherwise, a suggestion is to install g++: `sudo apt-get install g++`

- Install CMake: `sudo apt-get install cmake`

- OpenBLAS: `sudo apt-get install libopenblas-dev`    *(see the readme file provided with Armadillo for the other options, like MKL)*

- Armadillo (details in the readme file provided with Armadillo): 
    1. Download at http://arma.sourceforge.net;
    2. Extract files;
    3. In a terminal window, change into the directory that was created by unpacking Armadillo, and type `cmake .` (the full stop separated from "cmake" by a space is important);
    4. To generate the run-time armadillo library, type `make`;
    5. Type `sudo make install`

- I have run into a problem when Anaconda is installed. Armadillo could not found **libhdf5.so.101**, so I had to do one of the following:
    1. Add this line to **~/.bashrc**: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/lib` (in my case, that was `/home/user/anaconda3/lib`). Then, I needed to reload the .bashrc by running `source ~\.bashrc`. After that, METiS compiled and linked just fine. However, this solution may affect other application, so I don't recommend it.

    OR

    2. The other option is to include the following lines in the compiling and running scripts, in such a way that it does not affect other applications:
```    
    LD_LIBRARY_PATH=/path/to/lib
    export LD_LIBRARY_PATH
```    
        

    

### **Compiling and running**
After installing all the previous dependencies/third party tools listed above, METiS can be compiled, linked and run using the following commands (see the previous section in case you run in problems concerning **libhdf5.so.101**): 
```
    g++ -c src/*.cpp                        # Compile the source files
    g++ -o METiS -O2 *.o -larmadillo        # Link the resulting object files and Armadillo
    ./METiS <filePath>                      # Run METiS
```    

To make this process easier, two bash scripts named **build_metis.sh** and **run_metis.sh** are included. They can be run as
```
    .\build_metis             # Compile the source files and link with Armadillo + move METiS to bin folder
    .\run_metis fileName      # Run METiS with the input file specified by fileName
```    

This process will be conducted with CMake in the future (I don't know when).




## **Windows**
### **Installing dependencies/third party tools:**
- Download Visual Studio Community or other IDE of your choice.

- The Armadillo files that need to be included are provided in **/METiS/src/armadillo-8.400.0_include**. Alternatively, you can download Armadillo at http://arma.sourceforge.net and copy the entire **include** folder to a convenient location.

- Download and install Intel (R) Math Kernel Library (MKL) at https://software.intel.com/en-us/mkl (or other libraries like OpenBLAS).

- Modify **include/armadillo_bits/config.hpp** to indicate which libraries are currently available on your system (see the readme file provided with Armadillo for details).

- Configure a project in your IDE (the procedure below is given for Visual Studio):
    1. Create a new Empty Project or a new Windows Console Application;
    2. If you are using Visual Studio 2017 or above, downgrade your project to Visual Studio 2015. This is done by right clicking on the project name and selecting **Properties**. In **Configuration Properties**, find **Platform Toolset** and change it to **Visual Studio 2015 (v140)**. I honestly don't know why this is necessary;
    3. Add all the files included in **/src/** (except for Armadillo's **include** folder) to the project;       
    4. Locate the MKL folder. It is usually in **C:\Program Files (x86)\IntelSWTools\compiler_and_libraries_xxxx/mkl**;
    5. Tell the compiler to look for header files in Armadillo's and MKL's **include** folders. In Visual Studio 2017, right click on the project name and select properties. Unfold the **C/C++** list, pick **General** -> **Additional Include Directories**, and select **Edit**. Use the yellow button to insert a new entry and paste the path to **mkl\include** directory. Do the same thing for Armadillo's **include** folder;
    6. Next, go to **Linker** -> **General** -> **Additional Library Directories** and select **mkl\lib\intel64_win** (select intel32_win if your computer is 32bits);
    7. Finally, you need to go to **Linker** -> **Input** -> **Additional Dependencies** and add *mkl_core.lib*, *mkl_sequential.lib*, and *mkl_intel_lp64.lib*.


# METiS
METiS-USP (**M**orison's **E**quation **Ti**me domain **S**imulator - University of São Paulo) is a time-domain solver for the analysis of floating offshore wind turbines written in C++. By disabling some effects, it can also be applied to other floating structures or fixed offshore/onshore wind turbines.

The following effects are modeled in METiS-USP:
- Hydrodynamics: wave loads are computed using a modified version of Rainey's equation, which can be seen as an extended version of Morison's equation to include second-order wave loads acting on the floater. Drag on the structure is considered using the traditional quadratic term from Morison's equation. Details can be found in [Carmo, 2021](https://www.teses.usp.br/teses/disponiveis/3/3135/tde-03022022-120253/publico/LucasHenriqueSouzadoCarmoCorr21.pdf) and [Carmo and Simos, 2022](https://www.sciencedirect.com/science/article/pii/S0029801822012446).
- Aerodynamics: aerodynamic loads on the rotor are included using Blade Element Momentum Theory (BEMT). Details can be found in [Pegoraro, 2018](https://www.teses.usp.br/teses/disponiveis/3/3135/tde-31012019-075149/publico/BrunoPegoraroCorr18.pdf) (in Portuguese) and [Carmo, 2021](https://www.teses.usp.br/teses/disponiveis/3/3135/tde-03022022-120253/publico/LucasHenriqueSouzadoCarmoCorr21.pdf).
- For now, body dynamics is modeled using the equations of motion of a rigid body. Moorings are not modeled yet, but their effects can be partially included by specifying a stiffness matrix and a constant vertical force.  

# **Installing, compiling and running**

## **Linux**
### **Installing dependencies/third party tools:**
- Ensure a C++ compiler is installed. Otherwise, a suggestion is to install g++: `sudo apt-get install g++`

- Install CMake: `sudo apt-get install cmake`

- Instal Intel (R) MKL. By the time I wrote these Linux instructions, it is part of Intel (R) oneAPI Base Toolkit, available at https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.gdas9u

- Need to set environment variables: https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-dlfd-linux/top/before-you-begin.html


- Armadillo (details in the readme file provided with Armadillo):
    1. Download at http://arma.sourceforge.net;
    2. Extract files;
    3. As MKL is probably installed in a non standard location, I had to run this extra step that is not included in the Armadillo readme file `sudo ln -s /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_rt.so /usr/local/lib/libmkl_rt.so`. It creates a symbolic link so that `ld` can find MKL;
    4. In a terminal window, change into the directory that was created by unpacking Armadillo, and type `cmake .` (the full stop separated from "cmake" by a space is important);
    5. Activate MKL in your current terminal window source `/opt/intel/oneapi/setvars.sh`    
    7. Type `sudo make install`


### **Compiling and running**
After installing all the previous dependencies/third party tools listed above, METiS can be compiled and linked using `./build_metis`. It can then be run using `./METiS path/to/input/file.txt`. Note that to use Intel MKL you have to run `source /opt/intel/oneapi/setvars.sh` or equivalent depending on your MKL path.


## **Windows**
- Download Visual Studio Community or other IDE of your choice.

- The Armadillo files that need to be included are provided in **/METiS/src/armadillo-version_include**. Alternatively, you can download Armadillo at http://arma.sourceforge.net and copy the entire **include** folder to a convenient location.

- Download and install Intel (R) Math Kernel Library (MKL) at https://software.intel.com/en-us/mkl (Note: The following steps were written when MKL was not part of the OneAPI Base Toolkit, so they might need update).

- Modify **include/armadillo_bits/config.hpp** to indicate which libraries are currently available on your system (see the readme file provided with Armadillo for details).

- Configure a project in your IDE (the procedure below is given for Visual Studio):
    1. Create a new Empty Project;
    2. If you are using Visual Studio 2017 or above, downgrade your project to Visual Studio 2015. This is done by right clicking on the project name and selecting **Properties**. In **Configuration Properties**, find **Platform Toolset** and change it to **Visual Studio 2015 (v140)**. I honestly don't know why this is necessary;
    3. Add all the files included in **/src/** (except for Armadillo's **include** folder) to the project;       
    4. Locate the MKL folder. It is usually in **C:\Program Files (x86)\IntelSWTools\compiler_and_libraries_xxxx\windows\mkl**;
    5. Tell the compiler to look for header files in Armadillo's and MKL's **include** folders. In Visual Studio 2017, right click on the project name and select properties. Unfold the **C/C++** list, pick **General** -> **Additional Include Directories**, and select **Edit**. Use the yellow button to insert a new entry and paste the path to **mkl\include** directory. Do the same thing for Armadillo's **include** folder;
    6. Next, go to **Linker** -> **General** -> **Additional Library Directories** and select **mkl\lib\intel64_win** (select intel32_win if your computer is 32bits);
    7. Finally, you need to go to **Linker** -> **Input** -> **Additional Dependencies** and add *mkl_core.lib*, *mkl_sequential.lib*, and *mkl_intel_lp64.lib* (if 32bits, replace *mkl_intel_lp64.lib* by ** mkl_intel_c_dll.lib**).


### Running
METiS is a CLI software. To run it on Windows, open a PowerShell window and type:

`& "Path\to\METis\executable.exe" "Path\to\test_case.txt"`

with the paths replaced by the ones specific to your case.

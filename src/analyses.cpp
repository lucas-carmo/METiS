#include "analyses.h"

#include <iostream>
#include <armadillo>
#include <iomanip> // For input/output manipulators

// TEMPORARY: Needed usleep function to test output of progress bar
#ifdef __unix__
# include <unistd.h>
#elif defined _WIN32
# include <windows.h>
#define usleep(x) Sleep((x)/1000)
#endif

void timeDomainAnalysis(FOWT &fowt, ENVIR &envir)
{
    for ( ; envir.time() <= envir.timeTotal(); envir.stepTime() )    
    {  
        if (envir.time() == 0)
        {
            IO::print2outFile("TIME");

            IO::print2outFile("FLU_VEL_X");
			IO::print2outFile("FLU_VEL_Y");
			IO::print2outFile("FLU_VEL_Z");

			IO::print2outFile("FLU_ACC_X");
			IO::print2outFile("FLU_ACC_Y");
			IO::print2outFile("FLU_ACC_Z");

			IO::newLineOutFile();
			IO::print2outFile("(s)");

			IO::print2outFile("(m/s)");
			IO::print2outFile("(m/s)");
			IO::print2outFile("(m/s)");

			IO::print2outFile("(m/s^2)");
			IO::print2outFile("(m/s^2)");
			IO::print2outFile("(m/s^2)");
        }
        
        IO::newLineOutFile();
        IO::print2outFile(envir.time());

        IO::print2outFile(envir.fluidVel(0,0,0).at(0));
		IO::print2outFile(envir.fluidVel(0,0,0).at(1));
		IO::print2outFile(envir.fluidVel(0,0,0).at(2));

		IO::print2outFile(envir.fluidAcc(0, 0, 0).at(0));
		IO::print2outFile(envir.fluidAcc(0, 0, 0).at(1));
		IO::print2outFile(envir.fluidAcc(0, 0, 0).at(2));
		

        std::cout << round(100 * envir.time() / envir.timeTotal()) << "%" << '\r';
        std::fflush(stdout);
    }
}
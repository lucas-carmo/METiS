#include "analyses.h"

#include <iostream>
#include <armadillo>
#include <iomanip> // For input/output manipulators

void timeDomainAnalysis(FOWT &fowt, ENVIR &envir)
{
    for ( ; envir.time() <= envir.timeTotal(); envir.stepTime() )    
    {  
        if (envir.time() == 0)
        {
            IO::print2outFile("TIME");

            IO::print2outFile("FORCE_X");
			IO::print2outFile("FORCE_Y");
			IO::print2outFile("FORCE_Z");
			IO::print2outFile("MOMENT_X");
			IO::print2outFile("MOMENT_Y");
			IO::print2outFile("MOMENT_Z");

        }
        
        IO::newLineOutFile();
		IO::print2outFile(envir.time());
		arma::vec::fixed<6> hydroForce = fowt.hydrodynamicForce(envir);

		for (int jj = 0; jj < hydroForce.size(); ++jj)
		{
			IO::print2outFile(hydroForce.at(jj));
		}      

        std::cout << round(100 * envir.time() / envir.timeTotal()) << "%" << '\r';
        std::fflush(stdout);
    }
}
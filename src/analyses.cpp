#include "analyses.h"

#include <iostream>
#include <armadillo>
#include <iomanip> // For input/output manipulators

void timeDomainAnalysis(FOWT &fowt, ENVIR &envir)
{    
    IO::print2outLineHeader_turnOn(); // The header of the formatted output file is written during the first time step

    for ( ; envir.time() <= envir.timeTotal(); envir.stepTime() )    
    {          
        IO::print2outLine_turnOn();
		IO::print2outLine(envir.time());		

		arma::vec::fixed<6> hdForce = fowt.hydrodynamicForce(envir);        
		arma::vec::fixed<6> hsForce = fowt.hydrostaticForce(envir);

				

        // After the first time step, we do not need to print anything else to the header of the formatted output file
        if (envir.time() == 0)
        {
			IO::print2outLineHeader_turnOff();
			IO::printOutLineHeader2outFile();
        }

		IO::print2outLine_turnOff();
		IO::printOutLine2outFile();
        


		// Print progress to the screen
        std::cout << round(100 * envir.time() / envir.timeTotal()) << "%" << '\r';
        std::fflush(stdout);        
    }
}
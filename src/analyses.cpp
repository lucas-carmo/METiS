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

		arma::vec::fixed<6> hydroForce = fowt.hydrodynamicForce(envir);        
        IO::print2outLine(IO::OUTFLAG_HD_FORCE, hydroForce);        

        // After the first time step, we do not need to print anything else to the header of the formatted output file
        if (envir.time() == 0)
        {
            IO::print2outLineHeader_turnOff();            
        }

        

        std::cout << round(100 * envir.time() / envir.timeTotal()) << "%" << '\r';
        std::fflush(stdout);        
    }
}
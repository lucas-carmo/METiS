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
        arma::mat teste(1, 1, arma::fill::randu);        

        if (envir.time() == 0)
        {
            IO::print2outFile("TIME");
            IO::print2outFile("TesteClm1");
            IO::print2outFile("TesteClm2");
        }
        
        IO::newLineOutFile();
        IO::print2outFile(envir.time());
        IO::print2outFile(teste(0));
        IO::print2outFile(static_cast<int>(round(100 * envir.time() / envir.timeTotal())));


        std::cout << round(100 * envir.time() / envir.timeTotal()) << "%" << '\r';
        std::fflush(stdout);
        usleep(200);
    }
}
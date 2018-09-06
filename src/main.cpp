#include <iostream>
#include <chrono> // For std::chrono functions. It is useful to time the code and verify whether one method or another will be more performant
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include <string>
#include <stdexcept> // For std::exception

#include "IO.h" // IO is a pure static class that manages input/output
#include "FOWT.h"
#include "ENVIR.h"


// TEMPORARY: Needed usleep function to test output of progress bar
#ifdef __unix__
# include <unistd.h>
#elif defined _WIN32
# include <windows.h>
#define usleep(x) Sleep((x)/1000)
#endif

// METiS Version
extern const std::string g_METIS_VERSION{ "0.0.1" };

int main(int argc, char *argv[])
{
	std::cout << IO::METiS_Header() << '\n';
   
	try{
		if (argc != 2)
		{
			throw std::runtime_error("Please provide one input file");
			return 0;
		}

		FOWT fowt;
		ENVIR envir;

		std::cout << "Running METiS with file '" << argv[1] << "'\n";
		IO::setFiles(argv[1]); // Set paths to input file and output files
		IO::readInputFile(fowt, envir); // Read data from input file to fowt and envir
		IO::printSumFile(fowt, envir);

		std::cout << '\n';
		for (int ii = 1; ii <= 100; ++ii)
		{		
			std::cout << ii << "%" << '\r';		
			std::fflush(stdout);
			usleep(10000);
		}		
	}

	catch(std::exception &exception)
	{
		IO::writeErrorMessage( std::string(exception.what()) );
	}

	catch(...)
    {
        IO::writeErrorMessage( "Unknown exception thrown." );
    }


	std::cin.sync();
	std::cout << "\nPress enter to exit.\n";
	std::cin.get();

	return 0;
}




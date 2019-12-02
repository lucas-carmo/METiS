#include <iostream>
#include <string>
#include <stdexcept> // For std::exception

#include "IO.h" // IO is a pure static class that manages input/output
#include "FOWT.h"
#include "ENVIR.h"
#include "analyses.h"

// METiS Version
extern const std::string g_METIS_VERSION{ "0.0.1" };

int main(int argc, char *argv[])
{

	// Timer for measuring the elapsed time
	arma::wall_clock timer;
	timer.tic();

	std::cout << IO::METiS_Header() << '\n';
   
	// All the errors in the code are dealt by exceptions thrown in the different components
	try
	{ 
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
		IO::printSumFile(fowt, envir);	// Print the summary file for later verification of the input by the user

		// All the time domain calculation is done inside this function
		timeDomainAnalysis(fowt, envir);	
	}

	catch(std::exception &exception)
	{
		IO::print2log( "ERROR: " + std::string(exception.what()) );
	}

	catch(...)
    {
        IO::print2log( "ERROR: Unknown exception thrown." );
    }
	

	IO::print2log("Elapsed time: " + std::to_string(timer.toc()) + " seconds.");
	return 0;
}




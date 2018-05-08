#include <iostream>
#include <chrono> // For std::chrono functions. It is useful to time the code and verify whether one method or another will be more performant
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include <string>

#include "FOWT.h"
#include "IO.h" // IO is a pure static class that manages input/output
#include "ENVIR.h"



int main()
{
	FOWT fowt;
	ENVIR envir;

    std::string inFlNm = "Test/mainInputExample.txt";
    //std::string inFlNm = "../Test/mainInputExample.txt";

    IO::readInputFile(inFlNm, fowt, envir);

	int ii;
	std::cin >> ii;

	return 0;
}

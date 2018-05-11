#include "ENVIR.h"
#include "IO.h"

#include <iostream>
#include <vector>
//#include <stdlib.h>     // For 'atof' function, used for converting strings to double

void ENVIR::readTimeStep(const std::string &data)
{
    std::vector<std::string> input = stringTokenize(getData( data ), " \t");

    // 1) Verify if timeInput contains exactly 1 element (corresponding to timeStep)
    if ( input.size() != 1 )
        std::cout << "Deu ruim \n";

    // 2) Convert input from string to its corresponding numeric format (double, float, ...)
	if ( !string2num(input.at(0), m_timeStep) )
	{
        // Throw an exception if the conversion fails
		std::cout << "Deu ruim \n";
	}

    // 3) Verificar se os valores fazem sentido (se são positivos e maiores que zero, que timestep é menor que Tmax, etc)

    std::cout << "Time Step: " << m_timeStep << "\n";
}


void ENVIR::readTimeTotal(const std::string &data)
{
    std::vector<std::string> input = stringTokenize(getData( data ), " \t");

    // 1) Verify if timeInput contains exactly 1 element (corresponding to timeStep)
    if ( input.size() != 1 )
        std::cout << "Deu ruim \n";

    // 2) Convert input from string to its corresponding numeric format (double, float, ...)
	if ( !string2num(input.at(0), m_timeTotal) )
	{
        // Throw an exception if the conversion fails
		std::cout << "Deu ruim \n";
	}

    std::cout << "Total Time: " << m_timeTotal << "\n";
}



void ENVIR::readTimeRamp(const std::string &data)
{
    std::vector<std::string> timeInput = stringTokenize(getData( data ), " \t");

    // 1) Verify if timeInput contains exactly 1 element
    if (timeInput.size() != 1)
        std::cout << "Deu ruim \n";

    // 2) Convert input from string to its corresponding numeric format (double, float, ...)
	if ( !string2num(timeInput.at(0), m_timeRamp) )
	{
        // Throw an exception if the conversion fails
		std::cout << "Deu ruim \n";
	}

    // 3) Verificar se os valores fazem sentido (se são positivos e maiores que zero, que timestep é menor que Tmax, etc)
    std::cout << "Time Ramp: " << m_timeRamp << "\n";
}



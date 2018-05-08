#include "ENVIR.h"
#include "auxiliarFunctions.h"

#include <iostream>
#include <vector>
#include <stdlib.h>     // For 'atof' function, used for converting strings to double

void ENVIR::readTimeStepAndMaxTime(const std::string &data)
{
    std::vector<std::string> timeInput = stringTokenize(data, " \t");

    // 1) Verify if timeInput contains exactly 2 elements (corresponding to timeStep and maxTime)
    if (timeInput.size() != 2)
        std::cout << "Deu ruim \n";

    // 2) Convert inputs from string to their corresponding numeric format (double, float, ...)
	if ( !FromString(timeInput.at(0), m_timeStep) )
	{
		// Throw an exception if the conver
	}

	if ( !FromString(timeInput.at(1), m_maxTime) )
	{
		// Throw an exception if the conver
	}
	
    // 3) Verificar se os valores fazem sentido (se são positivos e maiores que zero, que timestep é menor que Tmax, etc)

    std::cout << "Time Step: " << m_timeStep << "\n";
    std::cout << "Max Time: " << m_maxTime << "\n";
}


#include "Floater.h"

#include "IO.h"
#include <iostream>
#include <vector>


void Floater::readMass(const std::string &data)
{
	readDataFromString(data, m_mass);
	std::cout << "Floater Mass: " << m_mass << "\n";
}

void Floater::readCoG(const std::string &data)
{
	// The coordinates of the center of gravity are separated by commas in the input string
	std::vector<std::string> input = stringTokenize(data, ",");

	if (input.size() != 3)
	{
		std::cout << "Deu ruim \n";
	}

	for (int ii = 0; ii < 3; ++ii)
	{
		readDataFromString(input.at(ii), m_CoG(ii));
	}

	std::cout << "CoG: \n" << m_CoG << "\n";
}



//mat Floater::rotatMat(const vec &FOWTpos)
//{}
//
//
//vec Floater::hydrodynamicForce(const ENVIR &envir, const vec &FOWTpos, const vec &FOWTvel, const vec &FOWTacc)
//{}


